
read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/indicator_taxa_for_dist.csv",header=FALSE,stringsAsFactors=FALSE)->taxlist
taxlist$V1<-gsub("_"," ",taxlist$V1)

get.det.dist<-function(listoftaxa){
	library(rgbif)
	colnames(listoftaxa)<-"TaxonSearched"
	listoftaxa[c("taxonKey","Occs","Pre1900","Pre1930","Pre1950","Pre1970","Pre1980","1800to1809","1810to1819","1820to1829","1830to1839","1840to1849","1850to1859","1860to1869","1870to1879","1880to1889","1890to1899","1900to1909","1910to1919","1920to1929","1930to1939","1940to1949","1950to1959","1960to1969","1970to1979","1980to1989","1990to1999","2000to2009","2010to2019")]<-NA
	for(i in 1:nrow(listoftaxa)){
		print(listoftaxa$TaxonSearched[i])
		pkey<-try(name_lookup(query=as.character(listoftaxa$TaxonSearched[i]),rank="species",return="data"),silent=TRUE)
		if(any(class(pkey)=="try-error")){
			}
		else{
			if(all(!names(pkey) %in% "nubKey")){
				}
			else{
				pkey<-subset(pkey,!is.na(pkey$nubKey))
				if(pkey$canonicalName[1]!=as.vector(listoftaxa$TaxonSearched[i])){			
					}
				else{		
					listoftaxa$taxonKey[i]<-pkey$nubKey[1]
					z<-occ_search(taxonKey=pkey$nubKey[1],hasCoordinate=TRUE,hasGeospatialIssue=FALSE,limit=200000,basisOfRecord="PRESERVED_SPECIMEN",geometry='POLYGON((-20 35, -20 90, 0 90, 40 90, 40 35, 0 35, -20 35))',eventDate="1200,2019",fields=c('name','key','year','month','day','recordedBy','decimalLatitude','decimalLongitude','country','catalogNumber','collectionCode','institutionCode','datasetKey','publishingOrg','stateProvince','locality'),return="data")
					if(is.null(nrow(z))){
						print(paste("........NO occurrences found for ", listoftaxa$TaxonSearched[i],sep=""))
					}
					if(!is.null(nrow(z))){
						write.csv(z,paste(file="/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020/",listoftaxa$TaxonSearched[i],"_gbif_records_10jan2020.csv",sep=""),row.names=FALSE)
						listoftaxa$Occs[i]<-nrow(z)
						listoftaxa$Pre1980[i]<-sum(z$year<1980)
						listoftaxa$Pre1970[i]<-sum(z$year<1970)
						listoftaxa$Pre1950[i]<-sum(z$year<1950)
						listoftaxa$Pre1930[i]<-sum(z$year<1930)
						listoftaxa$Pre1900[i]<-sum(z$year<1900)
						listoftaxa[i,"1800to1809"]<-nrow(z[z$year>=1800 & z$year<=1809,])
						listoftaxa[i,"1810to1819"]<-nrow(z[z$year>=1810 & z$year<=1819,])
						listoftaxa[i,"1820to1829"]<-nrow(z[z$year>=1820 & z$year<=1829,])
						listoftaxa[i,"1830to1839"]<-nrow(z[z$year>=1830 & z$year<=1839,])
						listoftaxa[i,"1840to1849"]<-nrow(z[z$year>=1840 & z$year<=1849,])
						listoftaxa[i,"1850to1859"]<-nrow(z[z$year>=1850 & z$year<=1859,])
						listoftaxa[i,"1860to1869"]<-nrow(z[z$year>=1860 & z$year<=1869,])
						listoftaxa[i,"1870to1879"]<-nrow(z[z$year>=1870 & z$year<=1879,])
						listoftaxa[i,"1880to1889"]<-nrow(z[z$year>=1880 & z$year<=1889,])
						listoftaxa[i,"1890to1899"]<-nrow(z[z$year>=1890 & z$year<=1899,])
						listoftaxa[i,"1900to1909"]<-nrow(z[z$year>=1900 & z$year<=1909,])
						listoftaxa[i,"1910to1919"]<-nrow(z[z$year>=1910 & z$year<=1919,])
						listoftaxa[i,"1920to1929"]<-nrow(z[z$year>=1920 & z$year<=1929,])
						listoftaxa[i,"1930to1939"]<-nrow(z[z$year>=1930 & z$year<=1939,])
						listoftaxa[i,"1940to1949"]<-nrow(z[z$year>=1940 & z$year<=1949,])
						listoftaxa[i,"1950to1959"]<-nrow(z[z$year>=1950 & z$year<=1959,])
						listoftaxa[i,"1960to1969"]<-nrow(z[z$year>=1960 & z$year<=1969,])
						listoftaxa[i,"1970to1979"]<-nrow(z[z$year>=1970 & z$year<=1979,])
						listoftaxa[i,"1980to1989"]<-nrow(z[z$year>=1980 & z$year<=1989,])
						listoftaxa[i,"1990to1999"]<-nrow(z[z$year>=1990 & z$year<=1999,])
						listoftaxa[i,"2000to2009"]<-nrow(z[z$year>=2000 & z$year<=2009,])
						listoftaxa[i,"2010to2019"]<-nrow(z[z$year>=2010 & z$year<=2019,])
					}
				}
			}
		}
	}		
write.csv(listoftaxa,file="/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_10jan2020.csv",row.names=FALSE)
}


#Now get records from GBIF
get.det.dist(taxlist)->taxlistdist1



#Subsequently focused study area

#Remove records outside focus area from GBIF files and save in 
#get the list of indicator taxa with data
tax<-read.csv(file="~/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_upd2_4jan2019.csv",stringsAsFactors=FALSE)
taxlist<-tax$TaxonSearched[!is.na(tax$Occs)]
taxlistnosp<-gsub(" ","_",taxlist)

#compared to 10jan2020 to make sure same taxa...yes, they are the same
#tax2<-read.csv(file="~/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_10jan2020.csv",stringsAsFactors=FALSE)
#taxlist2<-tax2$TaxonSearched[!is.na(tax2$Occs)]
#taxlistnosp2<-gsub(" ","_",taxlist2)
#taxlist==taxlist2
#taxlistnosp==taxlistnosp2

#define spatial extent
#ext<-c(5,25,45,55)

file.path1<-"~/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020/"
file.path<-"~/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020_w_colls_added_13april2020/"
for(x in 1:length(taxlist)){
	gb<-read.csv(file=paste(file.path1,taxlist[x],"_gbif_records_10jan2020.csv",sep=""),stringsAsFactors=FALSE)
	gb$DB<-NA
	nrow(gb)
	gb<-gb[gb$decimalLatitude>=45,]
	nrow(gb)
	gb<-gb[gb$decimalLatitude<=55,]
	nrow(gb)
	gb<-gb[gb$decimalLongitude>=5,]
	nrow(gb)
	gb<-gb[gb$decimalLongitude<=25,]
	nrow(gb)
	if(nrow(gb)>0){
		gb$DB[1:nrow(gb)]<-"GBIF"
	}
	write.csv(gb,file=paste(file.path,taxlist[x],"_gbif_records_10jan2020.csv",sep=""),row.names=FALSE)
}
#files have now been reduced and saved in a new folder







#get the list of indicator taxa with data
tax<-read.csv(file="~/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_upd2_4jan2019.csv",stringsAsFactors=FALSE)
taxlist<-tax$TaxonSearched[!is.na(tax$Occs)]
taxlistnosp<-gsub(" ","_",taxlist)
file.path<-"~/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020_w_colls_added_13april2020/"

#Read in Berlin data from Robert that was georeferenced
dd<-read.csv("/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/Berlin_collections_from_HTL/from_robert1_sorted_added4.csv",stringsAsFactors=FALSE)
nrow(dd)
#206

#drop those with dubious/unknown coordinates or without year or are duplicates with some in GBIF
dd<-dd[dd$Status %in% "Keep",]
nrow(dd)
#178

#Make sure restricted to area of interest
dd<-dd[dd$decimalLatitude>=45,]
nrow(dd)
dd<-dd[dd$decimalLatitude<=55,]
nrow(dd)
dd<-dd[dd$decimalLongitude>=5,]
nrow(dd)
dd<-dd[dd$decimalLongitude<=25,]
nrow(dd)




#Taxon
dd$Taxon<-paste(dd$Gen,dd$Spec,sep=" ")

#any taxa that are not in other list? Nope...all good
missing<-unique(dd$Taxon)[!unique(dd$Taxon) %in% taxlist]

#get list of taxa in munich
ddtt<-sort(unique(dd$Taxon))
ddtt


#add these records into files now and remove duplicates already in there (running number = GBIF CatalogNumber)
for(x in 1:length(ddtt)){
	gb<-read.csv(file=paste(file.path,ddtt[x],"_gbif_records_10jan2020.csv",sep=""),stringsAsFactors=FALSE)
	#add empty column for Berlin Barcode	
	gb$BBarcode<-NA
	
	#reduce dd to records of taxon of interest
	ddsp<-dd[dd$Taxon %in% ddtt[x],]
	nrow(ddsp)
		
	#remove ddsp records already in GBIF (ddsp running number = GBIF CatalogNumber)
	#reduce GBIF records to those held in Berlin
	cn<-gb$catalogNumber[gb$institutionCode=="B"]
	cn<-cn[!is.na(cn)]

	#also check for Barcode Number in Berlin (B 60 xxxxx)
	#trim down GBIF B 60... to be in format of Berlin database
	cn<-gsub(" /.*","",cn)
	
	#reduce ddsp to those with catalog numbers NOT in GBIF and from Berlin
	ddsp<-ddsp[!ddsp[,"Running.nr"] %in% cn,]
	nrow(ddsp)
	
	#reduce ddsp to those with catalog numbers NOT in GBIF and from Berlin
	ddsp<-ddsp[!ddsp[,"Barcode"] %in% cn,]
	nrow(ddsp)
	
	origrow<-nrow(gb)
	newrow<-nrow(ddsp)

	if(newrow>0){	
		#add extra rows
		gb[origrow+newrow,]<-NA
	
		#add in info
		for(z in 1:newrow){
			#gb$name[origrow+z]<-ddsp$Taxon[z]
			gb$decimalLongitude[origrow+z]<-ddsp$decimalLongitude[z]
			gb$decimalLatitude[origrow+z]<-ddsp$decimalLatitude[z]
			gb$year[origrow+z]<-ddsp$Collecting.year[z]
			gb$locality[origrow+z]<-ddsp$Loc....Hab.[z]
			gb$institutionCode[origrow+z]<-"BERLIN_DB"
			gb$catalogNumber[origrow+z]<-ddsp$Running.nr[z]
			gb$BBarcode[origrow+z]<-ddsp$Barcode[z]
			gb$recordedBy[origrow+z]<-ddsp$Collector.s.[z]
			gb$country[origrow+z]<-ddsp$Country[z]
			gb$month[origrow+z]<-ddsp$Datum[z]
			gb$DB[origrow+z]<-"BERLIN_DB"
		}
	}	
	gb_sorted<-gb[with(gb,order(-year)),]
	write.csv(gb_sorted,file=paste(file.path,ddtt[x],"_gbif_records_10jan2020.csv",sep=""),row.names=FALSE)
}










#manually double-checked whether any of the B collections Thorsten found were already in Berlin DB or GBIF...definitely some overlap
#add in HTL's from Berlin

#get the list of indicator taxa with data
tax<-read.csv(file="~/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_upd2_4jan2019.csv",stringsAsFactors=FALSE)
taxlist<-tax$TaxonSearched[!is.na(tax$Occs)]
taxlistnosp<-gsub(" ","_",taxlist)
file.path<-"~/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020_w_colls_added_13april2020/"

#Read in Berlin data from HTL that was georeferenced
dd<-read.csv("/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/Berlin_collections_from_HTL/Berlin_Collections_Climate_indicators_mpn_htl_updated2.csv",stringsAsFactors=FALSE)
nrow(dd)
#95

#drop those without year
dd<-dd[!is.na(dd$Year),]
nrow(dd)
#78

#drop those without coordinates or to drop for other reasons like already in one of the other db's
dd<-dd[dd$Status %in% "Keep",]
nrow(dd)
#49

#Make sure restricted to area of interest
dd<-dd[dd$decimalLatitude>=45,]
nrow(dd)
dd<-dd[dd$decimalLatitude<=55,]
nrow(dd)
dd<-dd[dd$decimalLongitude>=5,]
nrow(dd)
dd<-dd[dd$decimalLongitude<=25,]
nrow(dd)

#Taxon
dd$Taxon<-paste(dd$Gen,dd$Spec,sep=" ")

#any taxa that are not in other list? Nope...all good
missing<-unique(dd$Taxon)[!unique(dd$Taxon) %in% taxlist]

#get list of taxa in munich
ddtt<-sort(unique(dd$Taxon))
ddtt

#add these records into files now and remove duplicates already in there (running number = GBIF CatalogNumber)
for(x in 1:length(ddtt)){
	gb<-read.csv(file=paste(file.path,ddtt[x],"_gbif_records_10jan2020.csv",sep=""),stringsAsFactors=FALSE)
	
	#reduce dd to records of taxon of interest
	ddsp<-dd[dd$Taxon %in% ddtt[x],]
	nrow(ddsp)
		
	origrow<-nrow(gb)
	newrow<-nrow(ddsp)

	if(newrow>0){	
		#add extra rows
		gb[origrow+newrow,]<-NA
	
		#add in info
		for(z in 1:newrow){
			#gb$name[origrow+z]<-ddsp$Taxon[z]
			gb$decimalLongitude[origrow+z]<-ddsp$decimalLongitude[z]
			gb$decimalLatitude[origrow+z]<-ddsp$decimalLatitude[z]
			gb$year[origrow+z]<-ddsp$Year[z]
			gb$locality[origrow+z]<-ddsp$Locality[z]
			gb$institutionCode[origrow+z]<-"BERLIN_HTL"
			gb$recordedBy[origrow+z]<-ddsp$Collector[z]
			gb$country[origrow+z]<-ddsp$Country[z]
			gb$month[origrow+z]<-ddsp$Datum[z]
			gb$DB[origrow+z]<-"BERLIN_HTL"
		}
	}	
	gb_sorted<-gb[with(gb,order(-year)),]
	write.csv(gb_sorted,file=paste(file.path,ddtt[x],"_gbif_records_10jan2020.csv",sep=""),row.names=FALSE)
}


###Add in Munich collections from HTL


tax<-read.csv(file="~/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_upd2_4jan2019.csv",stringsAsFactors=FALSE)
taxlist<-tax$TaxonSearched[!is.na(tax$Occs)]
taxlistnosp<-gsub(" ","_",taxlist)
file.path<-"~/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020_w_colls_added_13april2020/"


#combine Thorsten's data from Munich w GBIF
munich<-read.csv("/Users/matthewnelsen/Documents/papers_reviews/papers/lichen_global_change_indicators/munich/M_Klimazeiger_Lumbsch_190326_mpn4_htl_updated2.csv",stringsAsFactors=FALSE)
nrow(munich)

#drop those with dubious/unknown coordinates
munich<-munich[munich$Status %in% "Keep",]
nrow(munich)

#drop those lacking a collection date
munich<-munich[!is.na(munich$Year),]
nrow(munich)
sort(unique(munich$Year))

#add colnames for lat/long
colnames(munich)[13:14]<-c("decimalLatitude","decimalLongitude")

#Taxon
munich$Taxon<-paste(munich$Gattung,munich$Epitheton,sep=" ")

#any taxa in Munich list that are not in other list? Yes...Parmelina_carporrhizans
missing<-unique(munich$Taxon)[!unique(munich$Taxon) %in% taxlist]

#drop this
nrow(munich)
munich<-munich[!munich$Taxon %in% missing,]
nrow(munich)

#get list of taxa in munich
mt<-sort(unique(munich$Taxon))
mt

#add these records into files now
for(x in 1:length(mt)){
	gb<-read.csv(file=paste(file.path,mt[x],"_gbif_records_10jan2020.csv",sep=""),stringsAsFactors=FALSE)
	msp<-munich[munich$Taxon %in% mt[x],]
	origrow<-nrow(gb)
	newrow<-nrow(msp)
	
	#add extra rows
	gb[origrow+newrow,]<-NA
	
	#add in info
	for(z in 1:newrow){
		#gb$name[origrow+z]<-msp$Taxon[z]
		gb$decimalLongitude[origrow+z]<-msp$decimalLongitude[z]
		gb$decimalLatitude[origrow+z]<-msp$decimalLatitude[z]
		gb$year[origrow+z]<-msp$Year[z]
		gb$locality[origrow+z]<-msp$Loc[z]
		gb$month[origrow+z]<-msp$Datum[z]
		gb$recordedBy[origrow+z]<-msp$Collector.s.[z]
		gb$institutionCode[origrow+z]<-"MUNICH_HTL"
		gb$catalogNumber[origrow+z]<-msp$Barcode[z]
		gb$DB[origrow+z]<-"MUNICH_HTL"
	}
	gb_sorted<-gb[with(gb,order(-year)),]
	write.csv(gb_sorted,file=paste(file.path,mt[x],"_gbif_records_10jan2020.csv",sep=""),row.names=FALSE)
}


#Hypotrachyna revoluta, Parmelia submontana, Parmelina quercina, and Punctelia jeckeri had Munich specimens already in GBIF - took care of



#also are some duplicates within GBIF.
########

#Lecanographa amylacea...retained the bottom two as they are separate collections in same herbarium.  Others seem like they are likely duplicates of these in different herbaria and with different coordinates
#1030438533	84797ee6-f762-11e1-a439-00145eb45e9a	11.52483	54.68247	NA	1952	4	12	Denmark	129490	Christiansen, Mogens Skytte	MSC	Lolland, Holeby, Bremersvold	Lichen	GBIF
#1042857695	aab0cf80-0c64-11dd-84d1-b8a03c50a862	11.533333	54.683333	Lolland	1952	4	12	Denmark	1687873	Mogens Skytte Christiansen	LD	Holleby	General	GBIF
#1148856044	8f7e3c45-4d76-4982-9478-651439e0cd4b	11.53333	54.68333	NA	1952	4	12	Denmark	TU63768	Mogens Skytte Christiansen	TU(M)	Lolland, Holeby, Bremersvold	NA	GBIF
#1212157073	85739778-f762-11e1-a439-00145eb45e9a	11.5333	54.6833	NA	1952	4	12	Denmark	186396	Mogens Skytte Christiansen	B	DENMARK, Lolland: Holeby, Bremersvold.	Lichen Herbarium Berlin	GBIF
#1324747120	fb1ecd28-f09e-4747-8bde-0b3d7a6f78d1	11.5333	54.6833	Lolland	1952	4	12	Denmark	95579	Christiansen, M.S.	BG	Holeby, Bremersvold	L	GBIF
#2237876536	84d26682-f762-11e1-a439-00145eb45e9a	11.52719	54.6828	NA	1952	4	12	Denmark	DMS-1013595	Mogens Skytte Christiansen	NA	Bremersvold	NA	GBIF
#2237878078	84d26682-f762-11e1-a439-00145eb45e9a	11.52719	54.6828	NA	1952	4	12	Denmark	DMS-1013596	Mogens Skytte Christiansen	NA	Bremersvold	NA	GBIF
#dropped "1030438533","1042857695","1148856044","1212157073","1324747120"

########

#Opegrapha vermicellifera
#2237904769	84d26682-f762-11e1-a439-00145eb45e9a	11.81666	54.78333	1943	5	28	Denmark	DMS-1039589	Mogens Skytte Christiansen	NA	Fuglsang Storskov	NA	NA	GBIF	NA
#914765	aab0cf80-0c64-11dd-84d1-b8a03c50a862	11.798333	54.722778	1943	NA	NA	Denmark	1052032	M. Skytte Christiansen	LD	Fuglsang	General	Lolland	GBIF	NA
#dropped "914765"

########

#Parmotrema perlatum - unclear where these are stored, but one will be removed later on when reduce to those with unique lat/long, so leave for now.
#2237883549	84d26682-f762-11e1-a439-00145eb45e9a	9.451133	54.84274	1939	8	17	Denmark	Kollund	DMS-1019762	Ove Almborn	NA	NA	NA	GBIF	NA
#2237885979	84d26682-f762-11e1-a439-00145eb45e9a	9.451133	54.84274	1939	8	17	Denmark	Kollund	DMS-1019763	Ove Almborn	NA	NA	NA	GBIF	NA

########

#Pertusaria hymenea - these seem like duplicates
#2237883187	84d26682-f762-11e1-a439-00145eb45e9a	9.591564	54.91952	NA	1939	8	15	Denmark	DMS-1019872	Ove Almborn	NA	Gr√•sten	NA	GBIF	NA
#788853441	aab0cf80-0c64-11dd-84d1-b8a03c50a862	9.594444	54.921944	Jylland	1939	8	15	Denmark	1637496	Ove Almborn	LD	Gr√•sten	General	GBIF	NA
#dropped "2237883187"

########

#Thelotrema lepadinum dropped "31326775","2270773269","2237905685","2237885415"
#2270704957	98bceab6-3c3e-4163-945d-6ca10576ebf1	10.9	47.538889	2004	9	5	Germany	Ammergauer Alpen, 15 km ESE of F√ºssen, uppermost valley of the rivulet Linder, Neualmgrie√ü	563562	W. Obermayer	ASU	Bavaria	NA	GBIF	NA
#31326775	835f8992-f762-11e1-a439-00145eb45e9a	10.9	47.540001	2004	9	5	Germany	Ammergauer Alpen, 15 km ESE of F√ºssen	5192-1	W. Obermayer	MAF	By	MAF-Lich	GBIF	NA

#1212249008	85739778-f762-11e1-a439-00145eb45e9a	14.6833	47.6917	1995	5	3	Austria	AUSTRIA, Steiermark, Ennstaler Alpen, 3.5 km E of St. Gallen, Wolfsbachgraben. Alt. 500 m.	138294	H. Mayrhofer & G. B√∂ttger	B	NA	Lichen Herbarium Berlin	GBIF	NA
#2270773269	98bceab6-3c3e-4163-945d-6ca10576ebf1	14.683333	47.691667	1995	5	3	Austria	Ennstaler Alpen, 3.5 km E of St. Gallen, Wolfsbachgraben, MTB 8354/1	503420	H. Mayrhofer; G. B√∂ttger	ASU	Styria	NA	GBIF	NA

#1068342	aab0cf80-0c64-11dd-84d1-b8a03c50a862	11.798333	54.722778	1943	5	28	Denmark	Fuglsang	1186402	M. Skytte-Christiansen	LD	Lolland	General	GBIF	NA
#2237905685	84d26682-f762-11e1-a439-00145eb45e9a	11.81666	54.78333	1943	5	28	Denmark	Fuglsang Storskov	DMS-1039628	Mogens Skytte Christiansen	NA	NA	NA	GBIF	NA

#1068265	aab0cf80-0c64-11dd-84d1-b8a03c50a862	9.835	54.902778	1939	8	15	Denmark	Als: S√∏nderskov	1109622	Ove Almborn	LD	Jylland	General	GBIF	NA
#2237885415	84d26682-f762-11e1-a439-00145eb45e9a	9.834309	54.903362	1939	8	15	Denmark	Als S√∏nderskov	DMS-1020562	Ove Almborn	NA	NA	NA	GBIF	NA

#these are potentially duplicates based on lat/long and date, but collector missing and locality info seems a bit different - retained
#918998577	828a3d8c-f762-11e1-a439-00145eb45e9a	10.9	47.533333	2004	9	5	Germany	Allg√§u, Lindertal	16580	NA	FR	NA	Herbarium Senckenbergianum (FR) - Fungi	GBIF	NA
#NA	NA	10.902493	47.540557	2004	5 Sep.	NA	GERMANY	GERMANY, Bayern (=Bavaria), Ammergauer Alpen (Ammergebirge), 15 km ESE of F√ºssen, uppermost valley of the rivulet Linder, Neualmgrie√ü. Alt. 1050-1100 m.	138522	W. Obermayer	BERLIN_DB	NA	NA	BERLIN_DB	
#dropped "31326775","2270773269","2237905685","2237885415"

########

#Usnea florida dropped "1273866841","1273868074"
#1099104857	7ba21c92-f762-11e1-a439-00145eb45e9a	11.448376	48.025852	NA	1979	3	18	Germany	M-0103162 / 538100 / 221203	Hertel, R.	SNSB-M	BSMlichenscoll	Deutschland, Bayern, Lkr. M√ºnchen, M√ºnchner Schotterebene: Forstenrieder Park: Eichenallee am s√ºdlichen Ende des Karolinenger√§umts (ca. 3 km NNE Baierbrunn) (MTB 7934/4); 600 m.	GBIF	NA
#1273866841	aab0cf80-0c64-11dd-84d1-b8a03c50a862	11.483333	48.020278	Bayern	1979	3	18	Germany	1884300	Rainer & H. Hertel	LD	General	Baierbrunn	GBIF	NA

#1099106006	7ba21c92-f762-11e1-a439-00145eb45e9a	14.377327	46.496559	NA	1979	7	NA	Austria	M-0196139 / 580844 / 252694	Follmann, B.	SNSB-M	BSMlichenscoll	S√ºd√∂sterreich, Karawanken: Hangwald unter den Matzent√ºrmen unweit Waidisch, 1500 m.	GBIF	NA
#1273868074	aab0cf80-0c64-11dd-84d1-b8a03c50a862	14.205278	46.439722	K√§rnten	1979	7	NA	Austria	1896456	G. Follmann & B. A. Follmann	LD	General	Karawanken	GBIF	NA
#dropped "1273866841","1273868074"

########

#Diploica canescens
#1212098084	85739778-f762-11e1-a439-00145eb45e9a	16.3833	47.3167	1991	4	11	Austria	138241	J. Hafellner & W. Maurer	B	AUSTRIA, Burgenland, G√ºnser Gebirge, 5 km WNW Rechnitz, vineyard S of Althodis. Alt. 440 m.	Lichen Herbarium Berlin	NA	GBIF	NA
#2270773538	98bceab6-3c3e-4163-945d-6ca10576ebf1	16.383333	47.316667	1991	4	11	Austria	503521	J. Hafellner; W. Maurer	ASU	G√ºnser Gebirge, 5 km WNW Rechnitz, vineyard S of Althodis	NA	Burgenland	GBIF	NA

#2237881942	84d26682-f762-11e1-a439-00145eb45e9a	10.71376	54.77292	1967	7	14	Denmark	DMS-1019465	Gunnar Degelius	NA	Magleby	NA	NA	GBIF	NA
#741254	aab0cf80-0c64-11dd-84d1-b8a03c50a862	10.713611	54.772778	1967	7	14	Denmark	1087658	Gunnar Degelius	LD	Magleby	General	Langeland	GBIF	NA
#dropped "2270773538", "2237881942"

########

#Melanohalea laciniatula
#2237865844	84d26682-f762-11e1-a439-00145eb45e9a	12.53299	54.97022	1939	8	6	Denmark	Store Klinteskov	DMS-1004630	Ove Almborn	NA	NA	NA	GBIF	NA
#885744	aab0cf80-0c64-11dd-84d1-b8a03c50a862	12.547778	54.970833	1939	8	6	Denmark	Store Klint	1021153	Ove Almborn	LD	General	M√∏n	GBIF	NA
#remove 885744

########

#Micarea adnata
#1212166962	85739778-f762-11e1-a439-00145eb45e9a	14.4333	47.8333	NA	1999	8	3	Austria	138376	R. T√ºrk	B	AUSTRIA, Ober√∂sterreich (=Upper Austria), Kalkalpen National Park, Reichraminger Hintergebirge, 5.8 km SSW of Reichraming, Z√∂belboden. Alt. 890 m.	Lichen Herbarium Berlin	GBIF
#2270772972	98bceab6-3c3e-4163-945d-6ca10576ebf1	14.433333	47.833333	Upper Austria	1999	8	3	Austria	ASUL000627	R. T√ºrk	ASU	Kalkalpen National Park, Reichraminger Hintergebirge, 5.8 km SSW of Reichraming, Z√∂belboden, MTB 8152	NA	GBIF
#removed 2270772972

########

#Parmelia submontana
#1056515003	7e380070-f762-11e1-a439-00145eb45e9a	9.55427	47.151846	1988	6	19	Liechtenstein	NA	BM001105095	Philippe Clerc	NHMUK	BOT	NA	GBIF	NA
#1632918387	aab0cf80-0c64-11dd-84d1-b8a03c50a862	9.557222	47.138611	1988	6	19	Liechtenstein	No locality information available	1895482	P. Clerc	LD	General	NA	GBIF	NA

#these might be the same, but are tough to tell because collector not listed in one and locations are different.  Leave both in.#######################
#1632918385	aab0cf80-0c64-11dd-84d1-b8a03c50a862	16.38088	51.07989	1951	9	23	Poland	No locality information available	1895391	Z. Tobolewski	LD	General	Dolnoslaskie	GBIF	NA
#912724377	828a3d8c-f762-11e1-a439-00145eb45e9a	16.336667	50.4625	1951	9	23	Poland	Sudeten, Gory Stolowe (pow. klodzki), ad viam Karlow - Kudowa	3884	NA	FR	Herbarium Senckenbergianum (FR) - Fungi	NA	GBIF	NA
#dropped 1632918387

########

#Parmelia quercina
#drop - same day, same location, but collector missing
#788525433	aab0cf80-0c64-11dd-84d1-b8a03c50a862	10.275833	47.408611	1950	8	2	Germany	Oberstdorf	1281961	E. Putzler	LD	General	Bayern	GBIF	NA
#912724433	828a3d8c-f762-11e1-a439-00145eb45e9a	10.264717	47.388883	1950	8	2	Germany	Allg√§u, Oberstdorf	3902	NA	FR	Herbarium Senckenbergianum (FR) - Fungi	NA	GBIF	NA
#dropped 912724433

########

#Punctelia jeckeri
#1135331489	c1a13bf0-0c71-11dd-84d4-b8a03c50a862	10.64732	47.78104	Bayern	2014	8	9	Germany	L-693328	Walter Obermayer	UPS	Schwaben (=Swabia), Alpenvorland, 2.2 km east of the centre of Marktoberdorf	BOT	GBIF	NA
#1262183055	85739778-f762-11e1-a439-00145eb45e9a	10.6473	47.781	NA	2014	8	9	Germany	193561	Walter Obermayer	B	GERMANY, Bayern (=Bavaria), Schwaben (=Swabia), Alpenvorland, 2.2 km east of the centre of Marktoberdorf. Alt. 770 m.	Lichen Herbarium Berlin	GBIF	NA
#2270735866	98bceab6-3c3e-4163-945d-6ca10576ebf1	10.643056	47.781042	Bavaria	2014	8	9	Germany	604027	Walter Obermayer	ASU	Schwaben (= Swabia), Alpenvorland, 2.2 km east of the centre of Marktoberdorf	NA	GBIF	NA
#drop 1135331489, 2270735866

#these look similar but are very different dates - leave both in:
#914864552	828a3d8c-f762-11e1-a439-00145eb45e9a	9.345283	48.794717	NA	1949	1	27	Germany	3964	NA	FR	Stetten/Remstal	Herbarium Senckenbergianum (FR) - Fungi	GBIF	NA
#788548945	aab0cf80-0c64-11dd-84d1-b8a03c50a862	9.340556	48.791667	Baden-W√ºrttemberg	1949	2	10	Germany	1306083	E. Putzler	LD	Stetten im Remstal	General	GBIF	NA

########

#Punctelia subrudecta
#possibly the same, but tough to know without collector?  seems like probably...barcode not same because other specimen not in Berlin. Dropped Berlin
#914864706	828a3d8c-f762-11e1-a439-00145eb45e9a	6.457217	49.92555	1957	4	13	Germany	3966	NA	FR	Dockendorf	Herbarium Senckenbergianum (FR) - Fungi	NA	GBIF	NA
#NA	NA	6.45652	49.925605	1957	17. April	NA	GERMANY	97687	Th. M√ºller	BERLIN_DB	DEUTSCHLAND, Rheinland-Pfalz: Kr. Bitburg, Dockendorf, an Pyrus malus.	NA	NA	BERLIN_DB	B 60 0141433

#likely the same - drop Berlin.
#788567629	aab0cf80-0c64-11dd-84d1-b8a03c50a862	22.22597	49.96282	1958	9	20	Poland	1325265	K. Glanc	LD	No locality information available	General	Podkarpackie	GBIF	NA
#914864469	828a3d8c-f762-11e1-a439-00145eb45e9a	22.666667	49.1375	1958	9	20	Poland	3961	NA	FR	Bieszczady Zachodnie, prope rivulum Wolosaty in vico Berezki	Herbarium Senckenbergianum (FR) - Fungi	NA	GBIF	NA
#NA	NA	22.681819	49.19664	1958	20.9.1958	NA	Polen	NA	K. Glanc	BERLIN_HTL	Wo≈Çosatym	NA	NA	BERLIN_HTL	NA
#Drop 788567629

#seem potentially identical but collector not given.  Treated as identical.
#2449353352	85739778-f762-11e1-a439-00145eb45e9a	7.6717	47.6333	1929	7	4	Germany	118851	V. J. Grummann	B	DEUTSCHLAND, Baden-W√ºrttemberg: Haagen bei L√∂rrach, R√∂ttler Wald.	Lichen Herbarium Berlin	NA	GBIF	NA
#914864614	828a3d8c-f762-11e1-a439-00145eb45e9a	7.68305	47.64055	1929	7	4	Germany	3965	NA	FR	Haagen (L√∂rrach)	Herbarium Senckenbergianum (FR) - Fungi	NA	GBIF	NA
#Drop 914864614

########

#Pyrenula nitida
#1212213392	85739778-f762-11e1-a439-00145eb45e9a	17.3667	50.3	1922	4	11	Poland	115742	V. J. Grummann	B	Altvatergebirge [Hrub√Ω Jesen√≠k], Ziegenhals [G≈Çucho≈Çazy], Hohenzollernwarte s√ºdlich.	Lichen Herbarium Berlin	NA	GBIF	NA
#914898702	828a3d8c-f762-11e1-a439-00145eb45e9a	17.375	50.302783	1922	4	11	Poland	540	NA	FR	Altvatergebirge, Bad Ziegenhals (Glucholazy)	Herbarium Senckenbergianum (FR) - Fungi	NA	GBIF	NA
#914898702	828a3d8c-f762-11e1-a439-00145eb45e9a	17.375	50.302783	1922	4	11	Poland	537	NA	FR	Altvatergebirge, Bad Ziegenhals (Glucholazy)	Herbarium Senckenbergianum (FR) - Fungi	NA	GBIF	NA
#dropped 914898702, 914898702

########


#now go through files and remove accessions with these keys...these are what appear to be duplicated collections within GBIF
drops<-c("1030438533","1042857695","1148856044","1212157073","1324747120","914765","2237883187","31326775","2270773269","2237905685","2237885415","1273866841","1273868074","2270773538","2237881942","885744","2270772972","1632918387","912724433","1135331489","2270735866","788567629","914864614","914898702","914898702")

tax<-read.csv(file="~/Documents/papers_reviews/papers/lichen_global_change_indicators/summary_w_latlong_upd2_4jan2019.csv",stringsAsFactors=FALSE)
taxlist<-tax$TaxonSearched[!is.na(tax$Occs)]
taxlistnosp<-gsub(" ","_",taxlist)
file.path<-"~/Documents/papers_reviews/papers/lichen_global_change_indicators/distribution_files_w_latlong_10jan2020_w_colls_added_13april2020/"

#add these records into files now
for(x in 1:length(taxlist)){
	gb<-read.csv(file=paste(file.path,taxlist[x],"_gbif_records_10jan2020.csv",sep=""),stringsAsFactors=FALSE)
	gb<-gb[!gb$key %in% drops,]
	gb_sorted<-gb[with(gb,order(-year)),]
	write.csv(gb_sorted,file=paste(file.path,taxlist[x],"_gbif_records_10jan2020.csv",sep=""),row.names=FALSE)
}

