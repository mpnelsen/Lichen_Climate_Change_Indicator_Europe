#Check for correlations among raster layers:

library(maptools)
library(mapproj)
library(rgeos)
library(rgdal)
library(ggplot2)
library(jsonlite)
library(RCurl)
library(raster)
library(RColorBrewer)
library(MASS)
library(magrittr)

#define spatial extent
ext<-c(5,25,45,55)

#read raster layers, crop and stack them.
for(x in 1:19){
	lsb1<-raster(paste("/Users/matthewnelsen/Documents/papers_reviews/papers/ants_ecogeography/new_enm/ascii_layers_for_maxent/bio_",x,".asc",sep=""))
	lsb1<-crop(lsb1,ext)
	if(x==1){
		lsb.stack<-lsb1
	}
	if(x>1){
		lsb.stack<-stack(lsb.stack,lsb1)	
	}
}

#saves cropped bioclim files
my.filenames<-paste("/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/raster_layers/croppedfor_study_area_bio_",1:19,".asc",sep="")
writeRaster(x=lsb.stack,filename=my.filenames,format="ascii",bylayer=TRUE)


#read and stack cropped and masked bioclim files
for(x in 1:19){
	lsb1<-raster(paste("/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/raster_layers/croppedfor_study_area_bio_",x,".asc",sep=""))
	if(x==1){
		lsb.stack<-lsb1
	}
	if(x>1){
		lsb.stack<-stack(lsb.stack,lsb1)	
	}	
}

#test for correlations among layers
require(ENMTools)
rc<-raster.cor.matrix(lsb.stack,method="spearman")
corrs.clean<-rc
corrs.clean[upper.tri(corrs.clean,diag=TRUE)]<-NA
write.csv(corrs.clean,file="/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/spearman_correlations_on_cropped_masked.csv",row.names=TRUE)


rc<-raster.cor.matrix(lsb.stack,method="pearson")
corrs.clean<-rc
corrs.clean[upper.tri(corrs.clean,diag=TRUE)]<-NA
write.csv(corrs.clean,file="/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/pearson_correlations_on_cropped_masked.csv",row.names=TRUE)

keeps<-c(1,2,3,4,12)








################




###Reduce files to include unique records up to 1969 and then from 1970-present and make a table
require(rgeos)
require(maptools)
require(raster)

#get list of taxa
tax<-read.csv(file="~/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_upd2_4jan2019.csv",stringsAsFactors=FALSE)
taxlist<-tax$TaxonSearched[!is.na(tax$Occs)]
taxlist
#Degelia plumbea and Parmotrema pseudoreticulatum had records in Europe, but not in focus area, so remove them
taxlist<-taxlist[!taxlist %in% c("Degelia plumbea","Parmotrema pseudoreticulatum")]
taxlist
taxlistnosp<-gsub(" ","_",taxlist)
file.path<-"~/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020_w_colls_added_13april2020/"

#make empty

cols<-c("Taxon","UniquePre1970","Unique1970Plus")
df<-data.frame(matrix(nrow=length(taxlist),ncol=length(cols)))
colnames(df)<-cols
df$Taxon<-taxlist

ext<-c(5,25,45,55)
lsb1<-raster(paste("/Users/matthewnelsen/Documents/papers_reviews/papers/ants_ecogeography/new_enm/ascii_layers_for_maxent/bio_",1,".asc",sep=""),values=TRUE)
lsb1<-crop(lsb1,ext)
projection(lsb1)<-CRS("+proj=longlat +datum=WGS84")
lsb1<-setMinMax(lsb1)




for(x in 1:length(taxlist)){
	fil<-read.csv(file=paste(file.path,taxlist[x],"_gbif_records_10jan2020.csv",sep=""),stringsAsFactors=FALSE)
	fil$b1test<-NA
	
	#add name
	fil$name<-taxlist[x]

	#add in BBarcode column if not present
	if(!"BBarcode" %in% colnames(fil)){
		fil$BBarcode<-NA
	}
	
	#only retain those inside restricted study area
	fil<-fil[fil$decimalLongitude>=5,]
	fil<-fil[fil$decimalLongitude<=25,]
	fil<-fil[fil$decimalLatitude>=45,]
	fil<-fil[fil$decimalLatitude<=55,]
	
	#if there are still records left...
	if(nrow(fil)==0){
		df[x,2]<-0
		df[x,3]<-0
	}
	
	if(nrow(fil)>0){
	
		#only retain those over bioclim raster layer
		xy<-SpatialPoints(cbind(fil$decimalLongitude,fil$decimalLatitude))
		fil$b1test<-extract(lsb1,xy)
		fil<-fil[!is.na(fil$b1test),]

		#partition into time periods
		fil.old<-fil[fil$year<=1969,]
	
		#remove duplicated latlong
		fil.old<-fil.old[!duplicated(paste(fil.old$decimalLongitude,fil.old$decimalLatitude,sep="_")),]
		write.csv(fil.old,file=paste(file.path,taxlist[x],"_pre1970_uniques_13april2020.csv",sep=""),row.names=FALSE)
	
		#number of records	
		df[x,2]<-nrow(fil.old)
	
		#partition into time periods
		fil.new<-fil[fil$year>1969,]
	
		#remove duplicated latlong
		fil.new<-fil.new[!duplicated(paste(fil.new$decimalLongitude,fil.new$decimalLatitude,sep="_")),]
		write.csv(fil.new,file=paste(file.path,taxlist[x],"_1970plus_uniques_13april2020.csv",sep=""),row.names=FALSE)
	
		#number of records	
		df[x,3]<-nrow(fil.new)	
	}
}


#save df with number of unique records pre/post
write.csv(df,file=paste(file.path,"uniques_prepost_cleaned_13april2020.csv",sep=""),row.names=FALSE)

#how many with 10 or more records?
df[df$UniquePre1970>9,]



######Make bias file

#essentially taken from here: http://stackoverflow.com/questions/27338512/color-countries-on-world-map-based-on-iso3-codes-in-r-using-ggplot
library(maptools)
library(mapproj)
library(rgeos)
library(rgdal)
library(ggplot2)
library(jsonlite)
library(RCurl)
library(raster)
library(RColorBrewer)
library(MASS)
library(magrittr)
#https://scottrinnan.wordpress.com/2015/08/31/how-to-construct-a-bias-file-with-r-for-use-in-maxent-modeling/

#define spatial extent
ext<-c(5,25,45,55)
lsb1<-raster(paste("/Users/matthewnelsen/Documents/papers_reviews/papers/ants_ecogeography/new_enm/ascii_layers_for_maxent/bio_",1,".asc",sep=""))
lsb1<-crop(lsb1,ext)

tax<-read.csv(file="~/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_upd2_4jan2019.csv",stringsAsFactors=FALSE)
taxlist<-tax$TaxonSearched[!is.na(tax$Occs)]
taxlist
#Degelia plumbea and Parmotrema pseudoreticulatum had records in Europe, but not in focus area, so remove them
taxlist<-taxlist[!taxlist %in% c("Degelia plumbea","Parmotrema pseudoreticulatum")]
taxlist
taxlistnosp<-gsub(" ","_",taxlist)
file.path<-"~/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020_w_colls_added_13april2020/"


for(x in 1:length(taxlist)){
	if(x==1){
		dat<-read.csv(file=paste(file.path,taxlist[x],"_1970plus_uniques_13april2020.csv",sep=""),stringsAsFactors=FALSE)
	}
	if(x>1){
		jj<-read.csv(file=paste(file.path,taxlist[x],"_1970plus_uniques_13april2020.csv",sep=""),stringsAsFactors=FALSE)
		dat<-rbind(dat,jj)
	}
}

#save uniques new
write.csv(dat,file=paste(file.path,"combined_1970plus_uniques_13april2020.csv",sep=""),row.names=FALSE)



for(x in 1:length(taxlist)){
	if(x==1){
		dat<-read.csv(file=paste(file.path,taxlist[x],"_pre1970_uniques_13april2020.csv",sep=""),stringsAsFactors=FALSE)
	}
	if(x>1){
		jj<-read.csv(file=paste(file.path,taxlist[x],"_pre1970_uniques_13april2020.csv",sep=""),stringsAsFactors=FALSE)
		dat<-rbind(dat,jj)
	}
}

#save uniques old
write.csv(dat,file=paste(file.path,"combined_pre1970_uniques_13april2020.csv",sep=""),row.names=FALSE)


dat<-dat[,c("decimalLongitude","decimalLatitude")]
dat.clean<-dat[!duplicated(paste(dat$decimalLongitude,dat$decimalLatitude,sep="_")),]
occur.ras<-rasterize(dat.clean,lsb1,1)

presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(ncol(lsb1),nrow(lsb1)), lims=c(extent(lsb1)[1:4]))
dens.ras <- raster(x=dens)
dens.ras
plot(dens.ras)
dens.ras.ext<-dens.ras
extent(dens.ras.ext)<-extent(lsb1)
dens.ras.ext
plot(dens.ras.ext)

dens.ras.ext.mask<-mask(dens.ras.ext,lsb1)
plot(dens.ras.ext.mask)

#writeRaster(dens.ras.ext.mask, "~/Desktop/ant_bias_file.tif")
writeRaster(x=dens.ras.ext.mask,filename="/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/bias_layer_13april2020.asc",format="ascii",bylayer=TRUE,overwrite=TRUE)



library(maptools) 
library(maps)
library(rgeos)
library(rgdal)
library(raster)
library(rnaturalearth)
junk<-map('world', fill = TRUE, col = "grey",plot=FALSE)
countries <- ne_countries(scale=110)
ext<-c(5,25,45,55)
countries<-crop(countries,ext)
bias<-raster("/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/bias_layer_13april2020.asc")
names(bias)<-"bias"
projection(bias)<-CRS("+proj=longlat +datum=WGS84")
bias<-setMinMax(bias)
plot(bias)
plot(countries, col=NA, add=TRUE)
points(dat.clean$decimalLongitude,dat.clean$decimalLatitude,col="red",cex=0.4,pch=16)

##########################
##########################
require(raster)
library(mapproj)
library(rgeos)
library(maptools)
library(rgdal)
library(ggplot2)
library(jsonlite)
library(RCurl)
library(raster)
library(RColorBrewer)
library(MASS)
library(magrittr)
require(dismo)
require(parallel)
require(doParallel)
require(foreach)
library(PresenceAbsence)
library(raster)
library(dismo)
library(zoo)
require(SDMTools)

filepath<-"/home/mpnelsen/lichens_global_change/"
ext<-c(5,25,45,55)

#load and stack rasters for present scenario
stack.rasters.regular<-function(rasterpathprefix){
	for(x in c(1,2,3,4,12)){
		lsb1<-raster(paste(filepath,"rasters/croppedfor_study_area_bio_",x,".asc",sep=""))
		names(lsb1)<-paste("bio_",x,sep="")
		projection(lsb1)<-CRS("+proj=longlat +datum=WGS84")
		lsb1<-setMinMax(lsb1)
		if(x==1){
			rast.stk<-lsb1
		}
		if(x>1){
			rast.stk<-stack(rast.stk,lsb1)	
		}	
	}
return(rast.stk)
}

present<-stack.rasters.regular(filepath)

bias<-raster(paste(filepath,"rasters/bias_layer_13april2020.asc",sep=""))
names(bias)<-"bias"
projection(bias)<-CRS("+proj=longlat +datum=WGS84")
bias<-setMinMax(bias)


dat<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_pre1970_uniques_13april2020.csv",stringsAsFactors=FALSE)
dat$name<-gsub(" ","_",dat$name)
occs<-dat[,c("name","decimalLongitude","decimalLatitude")]


#reduce to those with 10 records or more
colnames(occs)<-c("SPECIES","LONGITUDE","LATITUDE")
occs_s<-occs[occs$SPECIES %in% names(table(occs$SPECIES))[(table(occs$SPECIES)>9)],]
species<-unique(occs_s$SPECIES)




#this is a function to use parallel processing to predict suitability
#http://rstudio-pubs-static.s3.amazonaws.com/3911_c549d81ac2814caf9919dd0ec1b5829a.html
parallel.predictor<-function(no.cores,sets,model,climate.scenario,filenameclimscen,taxon){
	#make empties
	mod_allF_bias<-NULL
	lichenmaxStack_allF_bias<-NULL
	# set up parallel computing
	cl <- makeCluster(no.cores) # number of cores to use
	registerDoParallel(cl)
	#predict
	mod_allF_bias<-foreach(j=1:sets, .verbose=T, .packages=c("dismo","rJava","rgdal")) %dopar% predict(model[[j]], climate.scenario)
	lichenmaxStack_allF_bias<-stack(mod_allF_bias)
	names(lichenmaxStack_allF_bias)<-1:sets
	#write stacked predictions as R objects
	save(lichenmaxStack_allF_bias,file=paste("/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/",taxon,"_", filenameclimscen,".RData", sep=""))
  	#write write mean probability and sd of predictions as rasters
  	#maxent - ALL FEATURES - ***bias corrected***
 	writeRaster(sum(lichenmaxStack_allF_bias)/dim(lichenmaxStack_allF_bias)[3],paste("/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/",taxon,"_",filenameclimscen, "_maxent_mean_allF_biasCorrected_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
  	writeRaster(calc(lichenmaxStack_allF_bias, sd), paste("/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/",taxon,"_", filenameclimscen,"_maxent_sd_allF_biasCorrected_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	stopCluster(cl)	
	registerDoSEQ()
}

#now, run a for loop through species
for(sp in 1:length(species)){
	
	#subset
	occsb<-NULL
	occsb<-occs_s[occs_s$SPECIES %in% species[sp],]
	occsb<-occsb[,2:3]

	########### make sets of training/test (t/t) data
	ttSplit = 0.20 # ttSplit = test/train split percentage
	# function to partition data into t/t
	fold <- function(ttSplit){ 
		k = round(1/ttSplit, 0)
		fold <- kfold(occsb, k=k)
	}

	# make sets of t/t data
	sets = 50 # number of t/t sets
	folds <- replicate(sets, fold(ttSplit))
  
	# now loop through to build lists of t/t data
	lichenTrain <- list()
	lichenTest <- list()
	for(h in 1:sets){
		lichenTrain[[h]] <- occsb[folds[,h]!=1,]
		lichenTest[[h]] <- occsb[folds[,h]==1,]
	}

	# set up parallel computing
	cl <- makeCluster(12) # number of cores to use
	registerDoParallel(cl)

	#### MAXENT - DEFAULT FEATURES - BIAS CORRECTED BACKGROUND
	# run & predict (in parallel) maxent models for k randomizations
	# 10,000 locations selected using probabilistic target-group sampling from bias surface file
	lichenmaxMods_allF_bias <- foreach(k=1:sets, .verbose=T, .packages=c("dismo", "rJava", "rgdal")) %dopar% maxent(present, lichenTrain[[k]], a=SpatialPoints(randomPoints(bias, n=10000, prob=TRUE)), args=c("-P",c("-m", 10000)))
	save(lichenmaxMods_allF_bias,lichenTrain, lichenTest,file=paste("/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/",species[sp],"_","maxent_model_13april2020",".RData", sep=""))
	stopCluster(cl)
	registerDoSEQ()
	parallel.predictor(no.cores=10,sets=sets,model=lichenmaxMods_allF_bias,climate.scenario=present,filenameclimscen="present_13april2020",taxon=species[sp])
}
	
	








#require(SSDM)
library(mapproj)
library(rgeos)
library(maptools)
library(rgdal)
library(ggplot2)
library(jsonlite)
library(RCurl)
library(raster)
library(RColorBrewer)
library(MASS)
library(magrittr)
require(dismo)
require(parallel)
require(doParallel)
require(foreach)
library(PresenceAbsence)
library(zoo)
require(SDMTools)

#load and stack rasters for present scenario
stack.rasters.regular<-function(rasterpathprefix){
	for(x in c(1,2,3,4,12)){
		lsb1<-raster(paste(filepath,"rasters/croppedfor_study_area_bio_",x,".asc",sep=""))
		names(lsb1)<-paste("bio_",x,sep="")
		projection(lsb1)<-CRS("+proj=longlat +datum=WGS84")
		lsb1<-setMinMax(lsb1)
		if(x==1){
			rast.stk<-lsb1
		}
		if(x>1){
			rast.stk<-stack(rast.stk,lsb1)	
		}	
	}
return(rast.stk)
}

filepath<-"/home/mpnelsen/lichens_global_change/"
present<-stack.rasters.regular(filepath)

dat<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_pre1970_uniques_13april2020.csv",stringsAsFactors=FALSE)
dat$name<-gsub(" ","_",dat$name)
occs<-dat[,c("name","decimalLongitude","decimalLatitude")]


#reduce to those with 10 records or more
colnames(occs)<-c("SPECIES","LONGITUDE","LATITUDE")
occs_s<-occs[occs$SPECIES %in% names(table(occs$SPECIES))[(table(occs$SPECIES)>9)],]
species<-unique(occs_s$SPECIES)

sets = 50

#make a summary data frame for all of this
colnames_bigsum<-c("Taxon","AUC_mean","AUC_sd","meanProb_mean","meanProb_sd","meanBG_mean","meanBG_sd","Thresh95_mean","Thresh95_sd")
sum.df<-as.data.frame(matrix(nrow=length(species),ncol=length(colnames_bigsum)))
colnames(sum.df)<-colnames_bigsum

#now, run a for loop through species
for(sp in 1:length(species)){
	#loading this gives us the bg and train test info
	load(paste("/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/",species[sp],"_","maxent_model_13april2020",".RData", sep=""))
	#loading this gives us the stacked environment fit	
	load(paste("/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/",species[sp],"_","present_13april2020",".RData", sep=""))
	#separately calculate AUC etc. for present, and then do cross-comparisons
	#make a summary data frame
	colnames_ind<-c("Taxon","Replicate","AUC","meanProb","meanBG","thresh95","predArea95")
	ind.sum<-as.data.frame(matrix(nrow=sets,ncol=length(colnames_ind)))
	colnames(ind.sum)<-colnames_ind
	ind.sum$Taxon<-species[sp]
	ind.sum$Replicate<-1:sets
	#now calculate stats on fit model
	AUCmod <- thresh95 <- meanProb <- meanBG <- NULL  	
  	area <- cellStats(lichenmaxStack_allF_bias[[1]]<=1, sum)
  	#predicted probability at random background points 	
  	probBG <- extract(lichenmaxStack_allF_bias, SpatialPoints(randomPoints(present, 10000)))
  	for(ii in 1:dim(lichenmaxStack_allF_bias)[3]){    
    	# predicted probability at test points
    	probTest <- as.numeric(na.omit(extract(lichenmaxStack_allF_bias[[ii]], SpatialPoints(lichenTest[[ii]]))))
    	evalDismo <- evaluate(p=probTest, a=probBG[,ii])
    	AUCmod[[ii]] <- evalDismo@auc
		ind.sum$AUC[ii]<-AUCmod[[ii]]
		meanProb[[ii]] <- mean(probTest)
		ind.sum$meanProb[ii]<-meanProb[[ii]]
		meanBG[[ii]] <- mean(probBG[,ii])
    	ind.sum$meanBG[ii]<-meanBG[[ii]]
    	#made thresh95[ii] into thresh95[[ii]]
    	thresh95[[ii]] <- sort(probTest, decreasing=T)[round(length(probTest)*0.95,0)]
    	#made thresh95[ii] into thresh95[[ii]]
    	ind.sum$thresh95[ii]<-thresh95[[ii]]
	}	
	write.csv(ind.sum,file=paste("/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/",species[sp],"_present_evaluation_13april2020.csv", sep=""),row.names=FALSE)
	sum.df$Taxon[sp]<-species[sp]
	sum.df$AUC_mean[sp]<-mean(ind.sum$AUC)
	sum.df$AUC_sd[sp]<-sd(ind.sum$AUC)
	sum.df$meanProb_mean[sp]<-mean(ind.sum$meanProb)
	sum.df$meanProb_sd[sp]<-sd(ind.sum$meanProb)
	sum.df$meanBG_mean[sp]<-mean(ind.sum$meanBG)
	sum.df$meanBG_sd[sp]<-sd(ind.sum$meanBG)
	#added Thresh95
	sum.df$Thresh95_mean[sp]<-mean(ind.sum$thresh95)
	sum.df$Thresh95_sd[sp]<-sd(ind.sum$thresh95)
	write.csv(sum.df,file=paste("/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/summary_mean_model_evaluation_13april2020.csv", sep=""),row.names=FALSE)
}






















#USED THIS


require(raster)
require(rasterVis)
require(viridis)
require(gridExtra)
require(SDMTools)
require(lattice)
require(ggplot2)
require(grid)
require(Cairo)
library(maptools) 
library(maps)
library(rnaturalearth)
require(latticeExtra)
require(cowplot)
filepath<-"/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/"
sum.dat<-read.csv(file=paste(filepath,"summary_mean_model_evaluation_13april2020.csv",sep=""),stringsAsFactors=FALSE)
species<-sum.dat$Taxon

species<-species[c(3,5:7,11:12,14:15,1:2,4,9,8,10,13,16:17)]
sum.dat<-sum.dat[match(species, sum.dat$Taxon),]

historic.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_pre1970_uniques_13april2020.csv",stringsAsFactors=FALSE)
modern.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_1970plus_uniques_13april2020.csv",stringsAsFactors=FALSE)

junk<-map('world', fill = TRUE, col = "grey",plot=FALSE)
countries <- ne_countries(scale=110)
ext<-c(5,25,45,55)
countries<-crop(countries,ext)

library(calecopal)

# all palettes
names(cal_palettes)
mytheme <- rasterTheme(region = cal_palette("chaparral3"))



	historica<-historic.all[historic.all$name %in% gsub("_"," ",species[1]),]
	historica<-SpatialPoints(cbind(historica$decimalLongitude,historica$decimalLatitude))
	#modern points
	moderna<-modern.all[modern.all$name %in% gsub("_"," ",species[1]),]
	moderna<-SpatialPoints(cbind(moderna$decimalLongitude,moderna$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[1],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[1]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[1]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[1],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot1<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[1]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historica,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderna,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot1<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[1]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historica,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderna,pch=19,col="plum4"),columns=1)
	

	
	
	historicb<-historic.all[historic.all$name %in% gsub("_"," ",species[2]),]
	historicb<-SpatialPoints(cbind(historicb$decimalLongitude,historicb$decimalLatitude))
	#modern points
	modernb<-modern.all[modern.all$name %in% gsub("_"," ",species[2]),]
	modernb<-SpatialPoints(cbind(modernb$decimalLongitude,modernb$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[2],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[2]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[2]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[2],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot2<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[2]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicb,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernb,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot2<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[2]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicb,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernb,pch=19,col="plum4"),columns=1)

	
	
	
	historicc<-historic.all[historic.all$name %in% gsub("_"," ",species[3]),]
	historicc<-SpatialPoints(cbind(historicc$decimalLongitude,historicc$decimalLatitude))
	#modern points
	modernc<-modern.all[modern.all$name %in% gsub("_"," ",species[3]),]
	modernc<-SpatialPoints(cbind(modernc$decimalLongitude,modernc$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[3],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[3]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[3]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[3],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot3<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[3]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicc,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernc,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot3<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[3]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicc,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernc,pch=19,col="plum4"),columns=1)







	historicd<-historic.all[historic.all$name %in% gsub("_"," ",species[4]),]
	historicd<-SpatialPoints(cbind(historicd$decimalLongitude,historicd$decimalLatitude))
	#modern points
	modernd<-modern.all[modern.all$name %in% gsub("_"," ",species[4]),]
	modernd<-SpatialPoints(cbind(modernd$decimalLongitude,modernd$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[4],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[4]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[4]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[4],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot4<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[4]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicd,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernd,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot4<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[4]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicd,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernd,pch=19,col="plum4"),columns=1)




#this works well'ish...not for legend2...might have to just use colors and make new thing like before
historic.all$Period<-as.factor("Pre-1970's")
historic.all$Habitat<-as.factor("Suitable")
historic.all$HabCol<-1
modern.all$Period<-as.factor("1970-Present")
modern.all$Habitat<-as.factor("Unsuitable")
modern.all$HabCol<-0
tots<-rbind(modern.all,historic.all)

library(paletteer)
paletteer_d("calecopal::chaparral3",5)

p <- ggplot(tots, aes(decimalLatitude, decimalLongitude))+geom_point(aes(shape=Period,color=Period))+scale_color_manual(values=c("plum4", "darkgoldenrod"))+guides(colour = guide_legend(override.aes = list(size=3)))+ theme(legend.key = element_rect(fill = NA),legend.title.align=0.5, legend.title = element_text(face="bold"),legend.position="bottom") 
legend <- get_legend(p)


p2 <- ggplot(tots, aes(decimalLatitude, decimalLongitude))+geom_point(aes(color=Habitat))+scale_color_manual(values=c("#D3E3CAFF","#2F3525FF"))+guides(colour = guide_legend(override.aes = list(size=5,shape=15)))+ theme(legend.key = element_rect(fill = NA),legend.title.align=0.5, legend.title = element_text(face="bold"),legend.position="bottom") 
legend2 <- get_legend(p2)

Cairo(file=paste(filepath,"combo1thresh95_combined_13april2020_4.png",sep=""), type="png",bg="white",width = 4800, height = 4400, res=800, units = "px")
print(grid.arrange(present.threshplot1,present.threshplot2,present.threshplot3,present.threshplot4, legend2, legend, ncol=2, nrow = 3, layout_matrix = rbind(c(1,2), c(3,4), c(5,6)), widths = c(2.5, 2.5), heights = c(2, 2, 0.25)))
dev.off()


	Cairo(file=paste(filepath,"combo1_combined_13april2020_4.png",sep=""), type="png",bg="white",width = 4800, height = 4800, res=800, units = "px")
	print(grid.arrange(present.plot1,present.plot2,present.plot3,present.plot4,legend,ncol=2,nrow=3,layout_matrix = rbind(c(1,2), c(3,4), 5), widths = c(2.5, 2.5), heights = c(2, 2, 0.25)))
	dev.off()




























require(raster)
require(rasterVis)
require(viridis)
require(gridExtra)
require(SDMTools)
require(lattice)
require(ggplot2)
require(grid)
require(Cairo)
library(maptools) 
library(maps)
library(rnaturalearth)
require(latticeExtra)
require(cowplot)
filepath<-"/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/"
sum.dat<-read.csv(file=paste(filepath,"summary_mean_model_evaluation_13april2020.csv",sep=""),stringsAsFactors=FALSE)
species<-sum.dat$Taxon

species<-species[c(3,5:7,11:12,14:15,1:2,4,9,8,10,13,16:17)]
sum.dat<-sum.dat[match(species, sum.dat$Taxon),]


historic.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_pre1970_uniques_13april2020.csv",stringsAsFactors=FALSE)
modern.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_1970plus_uniques_13april2020.csv",stringsAsFactors=FALSE)

junk<-map('world', fill = TRUE, col = "grey",plot=FALSE)
countries <- ne_countries(scale=110)
ext<-c(5,25,45,55)
countries<-crop(countries,ext)

library(calecopal)

# all palettes
names(cal_palettes)
mytheme <- rasterTheme(region = cal_palette("chaparral3"))




	historice<-historic.all[historic.all$name %in% gsub("_"," ",species[5]),]
	historice<-SpatialPoints(cbind(historice$decimalLongitude,historice$decimalLatitude))
	#modern points
	moderne<-modern.all[modern.all$name %in% gsub("_"," ",species[5]),]
	moderne<-SpatialPoints(cbind(moderne$decimalLongitude,moderne$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[5],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[5]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[5]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[5],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot5<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[5]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historice,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderne,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot5<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[5]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historice,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderne,pch=19,col="plum4"),columns=1)







	historicf<-historic.all[historic.all$name %in% gsub("_"," ",species[6]),]
	historicf<-SpatialPoints(cbind(historicf$decimalLongitude,historicf$decimalLatitude))
	#modern points
	modernf<-modern.all[modern.all$name %in% gsub("_"," ",species[6]),]
	modernf<-SpatialPoints(cbind(modernf$decimalLongitude,modernf$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[6],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[6]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[6]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[6],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot6<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[6]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicf,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernf,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot6<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[6]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicf,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernf,pch=19,col="plum4"),columns=1)





	historicg<-historic.all[historic.all$name %in% gsub("_"," ",species[7]),]
	historicg<-SpatialPoints(cbind(historicg$decimalLongitude,historicg$decimalLatitude))
	#modern points
	moderng<-modern.all[modern.all$name %in% gsub("_"," ",species[7]),]
	moderng<-SpatialPoints(cbind(moderng$decimalLongitude,moderng$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[7],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[7]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[7]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[7],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot7<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[7]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicg,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderng,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot7<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[7]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicg,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderng,pch=19,col="plum4"),columns=1)






	historich<-historic.all[historic.all$name %in% gsub("_"," ",species[8]),]
	historich<-SpatialPoints(cbind(historich$decimalLongitude,historich$decimalLatitude))
	#modern points
	modernh<-modern.all[modern.all$name %in% gsub("_"," ",species[8]),]
	modernh<-SpatialPoints(cbind(modernh$decimalLongitude,modernh$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[8],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[8]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[8]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[8],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot8<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[8]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historich,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernh,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot8<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[8]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historich,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernh,pch=19,col="plum4"),columns=1)



#this works well'ish...not for legend2...might have to just use colors and make new thing like before
historic.all$Period<-as.factor("Pre-1970's")
historic.all$Habitat<-as.factor("Suitable")
historic.all$HabCol<-1
modern.all$Period<-as.factor("1970-Present")
modern.all$Habitat<-as.factor("Unsuitable")
modern.all$HabCol<-0
tots<-rbind(modern.all,historic.all)

library(paletteer)
paletteer_d("calecopal::chaparral3",5)

p <- ggplot(tots, aes(decimalLatitude, decimalLongitude))+geom_point(aes(shape=Period,color=Period))+scale_color_manual(values=c("plum4", "darkgoldenrod"))+guides(colour = guide_legend(override.aes = list(size=3)))+ theme(legend.key = element_rect(fill = NA),legend.title.align=0.5, legend.title = element_text(face="bold"),legend.position="bottom") 
legend <- get_legend(p)


p2 <- ggplot(tots, aes(decimalLatitude, decimalLongitude))+geom_point(aes(color=Habitat))+scale_color_manual(values=c("#D3E3CAFF","#2F3525FF"))+guides(colour = guide_legend(override.aes = list(size=5,shape=15)))+ theme(legend.key = element_rect(fill = NA),legend.title.align=0.5, legend.title = element_text(face="bold"),legend.position="bottom") 
legend2 <- get_legend(p2)

Cairo(file=paste(filepath,"combo2thresh95_combined_13april2020_4.png",sep=""), type="png",bg="white",width = 4800, height = 4400, res=800, units = "px")
print(grid.arrange(present.threshplot5,present.threshplot6,present.threshplot7,present.threshplot8, legend2, legend, ncol=2, nrow = 3, layout_matrix = rbind(c(1,2), c(3,4), c(5,6)), widths = c(2.5, 2.5), heights = c(2, 2, 0.25)))
dev.off()


	Cairo(file=paste(filepath,"combo2_combined_13april2020_4.png",sep=""), type="png",bg="white",width = 4800, height = 4800, res=800, units = "px")
	print(grid.arrange(present.plot5,present.plot6,present.plot7,present.plot8,legend,ncol=2,nrow=3,layout_matrix = rbind(c(1,2), c(3,4), 5), widths = c(2.5, 2.5), heights = c(2, 2, 0.25)))
	dev.off()

































require(raster)
require(rasterVis)
require(viridis)
require(gridExtra)
require(SDMTools)
require(lattice)
require(ggplot2)
require(grid)
require(Cairo)
library(maptools) 
library(maps)
library(rnaturalearth)
require(latticeExtra)
require(cowplot)
filepath<-"/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/"
sum.dat<-read.csv(file=paste(filepath,"summary_mean_model_evaluation_13april2020.csv",sep=""),stringsAsFactors=FALSE)
species<-sum.dat$Taxon

species<-species[c(3,5:7,11:12,14:15,1:2,4,9,8,10,13,16:17)]
sum.dat<-sum.dat[match(species, sum.dat$Taxon),]


historic.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_pre1970_uniques_13april2020.csv",stringsAsFactors=FALSE)
modern.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_1970plus_uniques_13april2020.csv",stringsAsFactors=FALSE)

junk<-map('world', fill = TRUE, col = "grey",plot=FALSE)
countries <- ne_countries(scale=110)
ext<-c(5,25,45,55)
countries<-crop(countries,ext)

library(calecopal)

# all palettes
names(cal_palettes)
mytheme <- rasterTheme(region = cal_palette("chaparral3"))





	historici<-historic.all[historic.all$name %in% gsub("_"," ",species[9]),]
	historici<-SpatialPoints(cbind(historici$decimalLongitude,historici$decimalLatitude))
	#modern points
	moderni<-modern.all[modern.all$name %in% gsub("_"," ",species[9]),]
	moderni<-SpatialPoints(cbind(moderni$decimalLongitude,moderni$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[9],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[9]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[9]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[9],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot9<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[9]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historici,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderni,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot9<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[9]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historici,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderni,pch=19,col="plum4"),columns=1)









	historicj<-historic.all[historic.all$name %in% gsub("_"," ",species[10]),]
	historicj<-SpatialPoints(cbind(historicj$decimalLongitude,historicj$decimalLatitude))
	#modern points
	modernj<-modern.all[modern.all$name %in% gsub("_"," ",species[10]),]
	modernj<-SpatialPoints(cbind(modernj$decimalLongitude,modernj$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[10],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[10]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[10]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[10],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot10<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[10]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicj,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernj,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot10<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[10]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicj,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernj,pch=19,col="plum4"),columns=1)









	historick<-historic.all[historic.all$name %in% gsub("_"," ",species[11]),]
	historick<-SpatialPoints(cbind(historick$decimalLongitude,historick$decimalLatitude))
	#modern points
	modernk<-modern.all[modern.all$name %in% gsub("_"," ",species[11]),]
	modernk<-SpatialPoints(cbind(modernk$decimalLongitude,modernk$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[11],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[11]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[11]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[11],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot11<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[11]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historick,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernk,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot11<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[11]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historick,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernk,pch=19,col="plum4"),columns=1)








	historicl<-historic.all[historic.all$name %in% gsub("_"," ",species[12]),]
	historicl<-SpatialPoints(cbind(historicl$decimalLongitude,historicl$decimalLatitude))
	#modern points
	modernl<-modern.all[modern.all$name %in% gsub("_"," ",species[12]),]
	modernl<-SpatialPoints(cbind(modernl$decimalLongitude,modernl$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[12],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[12]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[12]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[12],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot12<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[12]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicl,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernl,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot12<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[12]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicl,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernl,pch=19,col="plum4"),columns=1)



#this works well'ish...not for legend2...might have to just use colors and make new thing like before
historic.all$Period<-as.factor("Pre-1970's")
historic.all$Habitat<-as.factor("Suitable")
historic.all$HabCol<-1
modern.all$Period<-as.factor("1970-Present")
modern.all$Habitat<-as.factor("Unsuitable")
modern.all$HabCol<-0
tots<-rbind(modern.all,historic.all)

library(paletteer)
paletteer_d("calecopal::chaparral3",5)

p <- ggplot(tots, aes(decimalLatitude, decimalLongitude))+geom_point(aes(shape=Period,color=Period))+scale_color_manual(values=c("plum4", "darkgoldenrod"))+guides(colour = guide_legend(override.aes = list(size=3)))+ theme(legend.key = element_rect(fill = NA),legend.title.align=0.5, legend.title = element_text(face="bold"),legend.position="bottom") 
legend <- get_legend(p)


p2 <- ggplot(tots, aes(decimalLatitude, decimalLongitude))+geom_point(aes(color=Habitat))+scale_color_manual(values=c("#D3E3CAFF","#2F3525FF"))+guides(colour = guide_legend(override.aes = list(size=5,shape=15)))+ theme(legend.key = element_rect(fill = NA),legend.title.align=0.5, legend.title = element_text(face="bold"),legend.position="bottom") 
legend2 <- get_legend(p2)

Cairo(file=paste(filepath,"combo3thresh95_combined_13april2020_4.png",sep=""), type="png",bg="white",width = 4800, height = 4400, res=800, units = "px")
print(grid.arrange(present.threshplot9,present.threshplot10,present.threshplot11,present.threshplot12, legend2, legend, ncol=2, nrow = 3, layout_matrix = rbind(c(1,2), c(3,4), c(5,6)), widths = c(2.5, 2.5), heights = c(2, 2, 0.25)))
dev.off()


	Cairo(file=paste(filepath,"combo3_combined_13april2020_4.png",sep=""), type="png",bg="white",width = 4800, height = 4800, res=800, units = "px")
	print(grid.arrange(present.plot9,present.plot10,present.plot11,present.plot12,legend,ncol=2,nrow=3,layout_matrix = rbind(c(1,2), c(3,4), 5), widths = c(2.5, 2.5), heights = c(2, 2, 0.25)))
	dev.off()






























require(raster)
require(rasterVis)
require(viridis)
require(gridExtra)
require(SDMTools)
require(lattice)
require(ggplot2)
require(grid)
require(Cairo)
library(maptools) 
library(maps)
library(rnaturalearth)
require(latticeExtra)
require(cowplot)
filepath<-"/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/"
sum.dat<-read.csv(file=paste(filepath,"summary_mean_model_evaluation_13april2020.csv",sep=""),stringsAsFactors=FALSE)
species<-sum.dat$Taxon

species<-species[c(3,5:7,11:12,14:15,1:2,4,9,8,10,13,16:17)]
sum.dat<-sum.dat[match(species, sum.dat$Taxon),]


historic.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_pre1970_uniques_13april2020.csv",stringsAsFactors=FALSE)
modern.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_1970plus_uniques_13april2020.csv",stringsAsFactors=FALSE)

junk<-map('world', fill = TRUE, col = "grey",plot=FALSE)
countries <- ne_countries(scale=110)
ext<-c(5,25,45,55)
countries<-crop(countries,ext)

library(calecopal)

# all palettes
names(cal_palettes)
mytheme <- rasterTheme(region = cal_palette("chaparral3"))


	historicm<-historic.all[historic.all$name %in% gsub("_"," ",species[13]),]
	historicm<-SpatialPoints(cbind(historicm$decimalLongitude,historicm$decimalLatitude))
	#modern points
	modernm<-modern.all[modern.all$name %in% gsub("_"," ",species[13]),]
	modernm<-SpatialPoints(cbind(modernm$decimalLongitude,modernm$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[13],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[13]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[13]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[13],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot13<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[13]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicm,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernm,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot13<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[13]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicm,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernm,pch=19,col="plum4"),columns=1)









	historicn<-historic.all[historic.all$name %in% gsub("_"," ",species[14]),]
	historicn<-SpatialPoints(cbind(historicn$decimalLongitude,historicn$decimalLatitude))
	#modern points
	modernn<-modern.all[modern.all$name %in% gsub("_"," ",species[14]),]
	modernn<-SpatialPoints(cbind(modernn$decimalLongitude,modernn$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[14],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[14]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[14]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[14],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot14<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[14]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicn,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernn,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot14<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[14]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicn,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernn,pch=19,col="plum4"),columns=1)



	historico<-historic.all[historic.all$name %in% gsub("_"," ",species[15]),]
	historico<-SpatialPoints(cbind(historico$decimalLongitude,historico$decimalLatitude))
	#modern points
	moderno<-modern.all[modern.all$name %in% gsub("_"," ",species[15]),]
	moderno<-SpatialPoints(cbind(moderno$decimalLongitude,moderno$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[15],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[15]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[15]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[15],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot15<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[15]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historico,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderno,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot15<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[15]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historico,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(moderno,pch=19,col="plum4"),columns=1)






	historicp<-historic.all[historic.all$name %in% gsub("_"," ",species[16]),]
	historicp<-SpatialPoints(cbind(historicp$decimalLongitude,historicp$decimalLatitude))
	#modern points
	modernp<-modern.all[modern.all$name %in% gsub("_"," ",species[16]),]
	modernp<-SpatialPoints(cbind(modernp$decimalLongitude,modernp$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[16],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[16]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[16]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[16],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot16<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[16]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicp,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernp,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot16<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[16]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicp,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernp,pch=19,col="plum4"),columns=1)








	historicq<-historic.all[historic.all$name %in% gsub("_"," ",species[17]),]
	historicq<-SpatialPoints(cbind(historicq$decimalLongitude,historicq$decimalLatitude))
	#modern points
	modernq<-modern.all[modern.all$name %in% gsub("_"," ",species[17]),]
	modernq<-SpatialPoints(cbind(modernq$decimalLongitude,modernq$decimalLatitude))
	#read in prediction rasters
	present<-raster(paste(filepath,species[17],"_present_13april2020_maxent_mean_allF_biasCorrected_13april2020.tif",sep=""))
	#convert rasters to presence/absence based on threshold and write them
	presentThresh95<-present 	
	presentThresh95[presentThresh95>=sum.dat$Thresh95_mean[17]]<-1
	presentThresh95[presentThresh95<sum.dat$Thresh95_mean[17]]<-0
 	writeRaster(presentThresh95, paste(filepath,species[17],"_presentThresh95_13april2020.tiff", sep=""),type="GTiff",overwrite=T)
	#create plots
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.plot17<-levelplot(present,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[17]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=list(space="right",width=0.8),par.settings=viridisTheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicq,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernq,pch=19,col="plum4"),columns=1)
	#add this if have viridisTheme properly installed viridisTheme(region=viridis(10))
	present.threshplot17<-levelplot(presentThresh95,maxpixels=1e7,margin=FALSE,xlab=NULL,ylab=NULL,main=list(label=gsub("_"," ",species[17]),x=0.55,just="center",fontsize=12,fontface=4),colorkey=NULL,par.settings=mytheme)+latticeExtra::layer(sp.polygons(countries,col="white",lwd=0.15),columns=1)+latticeExtra::layer(sp.points(historicq,pch=17,col="darkgoldenrod"),columns=1)+latticeExtra::layer(sp.points(modernq,pch=19,col="plum4"),columns=1)




#this works well'ish...not for legend2...might have to just use colors and make new thing like before
historic.all$Period<-as.factor("Pre-1970's")
historic.all$Habitat<-as.factor("Suitable")
historic.all$HabCol<-1
modern.all$Period<-as.factor("1970-Present")
modern.all$Habitat<-as.factor("Unsuitable")
modern.all$HabCol<-0
tots<-rbind(modern.all,historic.all)

library(paletteer)
paletteer_d("calecopal::chaparral3",5)

p <- ggplot(tots, aes(decimalLatitude, decimalLongitude))+geom_point(aes(shape=Period,color=Period))+scale_color_manual(values=c("plum4", "darkgoldenrod"))+guides(colour = guide_legend(override.aes = list(size=3)))+ theme(legend.key = element_rect(fill = NA),legend.title.align=0.5, legend.title = element_text(face="bold"),legend.position="bottom") 
legend <- get_legend(p)


p2 <- ggplot(tots, aes(decimalLatitude, decimalLongitude))+geom_point(aes(color=Habitat))+scale_color_manual(values=c("#D3E3CAFF","#2F3525FF"))+guides(colour = guide_legend(override.aes = list(size=5,shape=15)))+ theme(legend.key = element_rect(fill = NA),legend.title.align=0.5, legend.title = element_text(face="bold"),legend.position="bottom") 
legend2 <- get_legend(p2)


Cairo(file=paste(filepath,"combo4thresh95_combined_13april2020_4.png",sep=""), type="png",bg="white",width = 4800, height = 6000, res=800, units = "px")
print(grid.arrange(present.threshplot13,present.threshplot14,present.threshplot15,present.threshplot16,present.threshplot17,legend2, legend, ncol=2, nrow = 4, layout_matrix = rbind(c(1,2), c(3,4), c(5,NULL), c(7,8)), widths = c(2.5, 2.5), heights = c(2, 2, 2, 0.25)))
dev.off()


	Cairo(file=paste(filepath,"combo4_combined_13april2020_4.png",sep=""), type="png",bg="white",width = 4800, height = 6000, res=800, units = "px")
	print(grid.arrange(present.plot13,present.plot14,present.plot15,present.plot16,present.plot17,legend,ncol=2,nrow=3,layout_matrix = rbind(c(1,2), c(3,4), c(5,6)), widths = c(2.5, 2.5), heights = c(2, 2, 2)))
	dev.off()

















#add proportion of modern records that are in historically unsuitable habitat

require(raster)
filepath<-"/home/mpnelsen/lichens_global_change/lichen_models_fit_13april2020/"

historic.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_pre1970_uniques_13april2020.csv",stringsAsFactors=FALSE)
modern.all<-read.csv(file="/home/mpnelsen/lichens_global_change/combined_1970plus_uniques_13april2020.csv",stringsAsFactors=FALSE)

dat<-read.csv(file="/home/mpnelsen/lichens_global_change/uniques_prepost_cleaned_13april2020.csv",stringsAsFactors=FALSE)
dat<-dat[dat$UniquePre1970>9,]
dat[,"Proportion Modern In Historically Unsuitable"]<-NA

for(x in 1:nrow(dat)){
	present<-raster(paste(filepath,gsub(" ","_",dat$Taxon[x]),"_presentThresh95_13april2020.tif",sep=""))
	recents<-modern.all[modern.all$name %in% dat$Taxon[x],c("decimalLongitude","decimalLatitude")]
	vals<-extract(present, SpatialPoints(recents))
	dat[x,"Proportion Modern In Historically Unsuitable"]<-1-round(sum(vals)/length(vals),2)
}

write.csv(dat,file="/home/mpnelsen/lichens_global_change/uniques_prepost_cleaned_wprop_13april2020.csv",row.names=FALSE)





