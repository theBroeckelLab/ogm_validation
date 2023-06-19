library(tidyverse)
library(archive)


fileList <- list.files("D:/Keyhole",full.names= TRUE)

for (f in fileList){
  #f <- "W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//BN_output_test/MCW_SVI_0001_-_De_novo_pipeline_results.zip"
  bioname<- sub("_-_De_novo_rerun_pipeline_results.zip.*","",f)
  bioname<- sub("_-_De_novo_pipeline_results.zip.*","",bioname)
  bioname <- sub(".*D:/Keyhole/","",bioname)
  
  
  smap <- read_tsv(archive_read(f,file = "output/contigs/annotation/variants_combine_filters_inMoleRefine1.smap"))
  headstart <- smap %>% pull(var = 1) %>% str_which("SmapEntryID")
  
  
  smap <- read_tsv(archive_read(f,file = "output/contigs/annotation/variants_combine_filters_inMoleRefine1.smap"),skip = headstart)
  smap <- dplyr::rename(smap,"Id"="#h SmapEntryID")
  headers<-colnames(smap)
  smap <- read_tsv(archive_read(f,file = "output/contigs/annotation/variants_combine_filters_inMoleRefine1.smap"),col_names = headers, skip = (headstart+2))
  #cnv <- read_tsv(archive_read(f,file = "output/contigs/annotation/CNVs_annotation/cnv_calls_exp_annotation_results.txt"), skip =3)
  smap$Id=paste0("SMAP", smap$Id)
  smap<- filter(smap,!(grepl('trans',Type)))
  smap<- filter(smap,!(grepl('inversion',Type)))
  smap <- smap %>% filter(Confidence >= 0)
  smap <- smap %>% filter(Self_molecule_count >= 5)
  smap <- smap %>% mutate(Sample_ID=bioname)
  smap <- smap %>% filter(Size >= 1500)
  
  smap <- dplyr::rename(smap,"Start"="RefStartPos")
  smap <- dplyr::rename(smap,"End"="RefEndPos")
  smap <- dplyr::rename(smap,"Chromosome"="RefcontigID1")
  smap <- smap %>% mutate(diff_size = (smap$Size/(smap$End-smap$Start)))
  smap <- smap %>% mutate(diff_labels = (smap$End-smap$Start))
  
  smap <- smap %>% filter()
  
  label_Mastersheet <- read_csv("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//label_Mastersheet.csv",col_types = "ccidddcdd")
  
  smap <- smap %>% select(Sample_ID,Id,Chromosome,Start,End,Size,Type,diff_size,diff_labels)
  
  label_Mastersheet <- label_Mastersheet %>% full_join(smap)
  
  label_Mastersheet <- label_Mastersheet[!duplicated(label_Mastersheet), ]
  
  write_csv(label_Mastersheet,"W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//label_Mastersheet.csv")
  
}

ggplot2
