library(openxlsx)
library(tidyverse)
###############################################################################################
## Identification of calls specific to bng (not detected in cma) from comparison reports ######
###############################################################################################
overlap_filelist <- list.files("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Validation_comparison_output/Samples_Discussed/")
overlap_filelist
for (x in overlap_filelist){
## read in sample report

xl.input=paste0("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Validation_comparison_output/Samples_Discussed/",x,"/report_summary-",x,".xlsx",sep="")
cma=read.xlsx(xl.input, sheet=4)
bng.pass=read.xlsx(xl.input,sheet=5)
bng.fail=read.xlsx(xl.input,sheet=6)
bngXcma=read.xlsx(xl.input, sheet=2)
#need to adjust sheets to possibly combine results of passed and failed and in roder to get around filtering for control %
bng.all<- rbind(bng.fail,bng.pass)
## pull only reportable cma calls
cma.reportable=cma[which(cma$Size>=200000),]


## replace default bng size with the actual difference between coordinates
#bng.pass$Size=bng.pass$End-bng.pass$Start
#Need to figure out if should be adjusted

## pull only reportable bng calls
idx.cnv=which((bng.pass$Type_general=="<DEL>"|bng.pass$Type_general=="<DUP>"|bng.pass$Type_general=="<INS>") & bng.pass$Size>=200000)
idx.sv=which(bng.pass$Type_general=="<INV>" & bng.pass$Size>=5000000)
bng.pass.reportable=bng.pass[c(idx.cnv, idx.sv),]

#Testing for pulling on all
## pull only reportable bng calls
idx.cnv_all=which((bng.all$Type_general=="<DEL>"|bng.all$Type_general=="<DUP>"|bng.all$Type_general=="<INS>") & bng.all$Size>=200000)
idx.sv_all=which(bng.all$Type_general=="<INV>" & bng.all$Size>=5000000)
bng.all.reportable=bng.all[c(idx.cnv_all, idx.sv_all),]

#filter from all
idx.all.reportable <- which((bng.all.reportable$Present_in_BNG_control_samples<=1 & bng.all.reportable$Present_in_BNG_control_samples_w_same_enzyme<=1)|is.na(bng.all.reportable$Present_in_BNG_control_samples))
bng.all.reportable <- bng.all.reportable[idx.all.reportable,]
## do we want to only consider non-masked bionano calls?
bng.all.reportable <- bng.all.reportable[which(bng.all.reportable$Maskedoverlap<.30),]


## determine if any remaining bng reportable calls are detected in cma
bng.all.reportable$Id           ##reportable BNG calls
bngXcma$Corresponding_BNG.ID     ##reportable BNG calls detected by CMA
#calls_specificToBionano=setdiff(bng.all.reportable$Id, bngXcma$Corresponding_BNG.ID[which(!is.na(bngXcma$Corresponding_BNG.ID))])
if(is.character(bngXcma$Corresponding_BNG.ID[which(!is.na(bngXcma$Corresponding_BNG.ID))])){
  calls_specificToBionano=setdiff(bng.all.reportable$Id, unlist(strsplit(bngXcma$Corresponding_BNG.ID[which(!is.na(bngXcma$Corresponding_BNG.ID))],";")))
} else{
  calls_specificToBionano=setdiff(bng.all.reportable$Id, bngXcma$Corresponding_BNG.ID[which(!is.na(bngXcma$Corresponding_BNG.ID))])
}

calls_specificToBionano          ##reportable BNG calls NOT detected by CMA


## pull information for unique bionano calls
bng_unique <- bng.all.reportable[which(bng.all.reportable$Id%in%calls_specificToBionano),]

#attach sample_id to unique calls
sample_name <- sub("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Validation_comparison_output/Samples_Discussed/","",xl.input)
sample_name <- sub(".report.*","",sample_name)
sample_name
dim(bng_unique)
if (dim(bng_unique)[1]!= 0){
  bng_unique$Sample_name <- sample_name

  bionanoUnique <- read.xlsx("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Bionano_unique_calls/bionano_unique.xlsx")
  bionanoUnique <- rbind(bionanoUnique,bng_unique)
  write.xlsx(bionanoUnique,"W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Bionano_unique_calls/bionano_unique.xlsx")
  }
}

#filtering the excel
bionanoUnique <- read.xlsx("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Bionano_unique_calls/bionano_unique.xlsx")
bionanoUnique
bionanoUnique <- bionanoUnique[which(bionanoUnique$Confidence!=0),]

bionanoUnique <- bionanoUnique %>% filter(!(Type=="gain" & Confidence<0.99 ))
bionanoUnique <- bionanoUnique %>% filter(!(Type=="loss" & Confidence<0.99 ))
bionanoUnique <- bionanoUnique %>% distinct(Sample_name,Start,End,Size, .keep_all = TRUE)
write.xlsx(bionanoUnique,"W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Bionano_unique_calls/bionano_unique.xlsx")

bionanoUnique_filtered <- read.xlsx("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Bionano_unique_calls/bionano_unique_filtered.xlsx")
