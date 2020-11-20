require(rgbif)
require(stringr)
filepath<-"...distribution_files_w_latlong_10jan2020/"
filepath.new<-"...distribution_files_w_latlong_10jan2020_DOI_GBIF_CITATION_ADDED_20NOV2020/"
files<-list.files(filepath)
species<-gsub("_.*","",files)
unique.dois<-NULL

for(z in 1:length(species)){
	xx<-NULL
	xx<-read.csv(file=paste(filepath,species[z],"_gbif_records_10jan2020.csv",sep=""),stringsAsFactors=FALSE)
	xx$notes<-"proper GBIF citations were not initially included here and were added subsequently on 20 November 2020"
	xx$citation<-NA
	xx$doi<-NA
	for(q in 1:nrow(xx)){
		xx$citation[q]<-gbif_citation(x=xx$datasetKey[q])$citation$citation
		xx$doi[q]<-gsub(pattern=" a",replacement="",str_extract(xx$citation[q],"https://doi.org/.* a"))
	}
	unique.dois<-sort(unique(c(xx$doi,unique.dois)))
	write.csv(xx,file=paste(filepath.new,species[z],"_gbif_records_10jan2020_DOI_GBIF_CITATION_ADDED_20NOV2020.csv",sep=""),row.names=FALSE)
}


length(unique.dois)
write.csv(unique.dois,"...dois_for_gbif_data_used.csv",row.names=FALSE)