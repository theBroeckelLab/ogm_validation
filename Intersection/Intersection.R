#Intersect Pipeline
#Array output files from ChAS 4.2
#7/7/2023
library(tidyverse)
library(archive)
library(GenomicRanges)
#Paths to folders not direct files
fileListbn <- list.files("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//BN_output_test",full.names= TRUE)
fileListArray <- list.files("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//Array_output_test",full.names= TRUE)

##Check for files found
for (f in fileListbn){
  print(f)
}
for (f in fileListArray){
  print(f)
}


#Intersection of Array and Bionano results
for (x in fileListArray){
  array_input <- read_tsv(x)
  #Strips to name of file for mathcing with corropsonding array file
  #Need to change based on used folder locations
  arrayname <- sub(".*W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//Array_output_test/","",x)
  arrayname <- sub(".*?_","",arrayname)
  arrayname <- sub(".cyhd.cychp.segments.txt*","",arrayname)
  #rename omim gene in array to remove special characters
  array_input <- array_input %>% dplyr::rename("Array_OMIM_genes"=`OMIM ® Genes`)
  array_input <- array_input %>% dplyr::rename("Array_OMIM_genes_count"=`OMIM ® Genes Count`)
  dummyname <- arrayname
  dummysub <- substring(dummyname,1,2)
  if (dummysub == "NA"){
    dummyname <- paste("HM",substring(dummyname,3,),sep = "")
  }
  #Arrange columns for better readability
  array_input <- array_input %>%
    relocate("Chromosome", .before = "File")
  array_input <- array_input %>%
    relocate("Min", .after = "Chromosome")
  array_input <- array_input %>%
    relocate("Max", .after = "Min")
  #Remove any AOH/LOH
  array_input <- array_input %>%
    drop_na("Chromosome")
  
  #Rename size(kb) to size
  array_input <- array_input %>%
    dplyr::rename(`Size` = `Size (kbp)`)
  #Change size from kbp to bp
  array_input <- array_input %>%
    mutate(`Size` = `Size`*1000)
  #add IDs to array
  array_input$Call_ID=paste0("CMA", 1:nrow(array_input))
  
  #change X to 23 and Y to 24
  array_input <- array_input %>% mutate(Chromosome = replace(Chromosome, Chromosome == "X","23"))
  array_input <- array_input %>% mutate(Chromosome = replace(Chromosome, Chromosome == "Y","24"))
  #Add Masking to the Arrays before intersection with full combined
  #Make grange object for Array
  gra=GRanges(seqnames=array_input$Chromosome, ranges=paste0(array_input$Min,"-",array_input$Max), strand=NULL) 
  ##INTERSECT WITH MASK FILE
  mask=read.csv("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/vcf/Corriell_vcf/bedfiles/bionano_annotation_files/hg38_CombinedMask_cnv-nbase-segdup-common.txt", sep="\t", header=F)
  gr3=GRanges(seqnames=mask$V1, ranges=paste0(mask$V2,"-",mask$V3), strand=NULL) 
  gr3.reduced=IRanges::reduce(gr3)
  ##overlap BED with masks
  hits=findOverlaps(gra, gr3.reduced)
  ##compile into one data frame
  comb=cbind(as.data.frame(gra[queryHits(hits)]), variantID=array_input$Call_ID[as.data.frame(hits)[,1]], as.data.frame(gr3.reduced[subjectHits(hits)]))
  colnames(comb)=c("sv_chr","sv_start","sv_end","sv_size","sv_strand","sv_ID","mask_chr","mask_start","mask_end","mask_size","mask_strand")
  #filter based on percent overlap
  comb$percentOverlap <- width(pintersect(gra[queryHits(hits)], gr3.reduced[subjectHits(hits)])) / width(gra[queryHits(hits)])
  ## map genes to svs
  svs_map2mask=unique(comb$sv_ID)
  array_input$Maskedoverlap=0
  for (i in 1:length(svs_map2mask)) {
    array_input$Maskedoverlap[which(array_input$Call_ID==svs_map2mask[i])]=sum(comb$percentOverlap[which(comb$sv_ID==svs_map2mask[i])])
  }
  
  ##############################################################################################################
  #INTERSECT CMA WITH MERGED MORBID TRACK FROM BIONANO
  ##read in merged morbid
  mergedMorbid=read.csv("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/vcf/Corriell_vcf/bedfiles/bionano_annotation_files/MERGED_morbid_GRCh38_GTF.gtf", sep="\t", header=F)
  mergedMorbid$gene=gsub("gene_id ", "", sapply(strsplit(mergedMorbid$V9, "; "), `[`, 1))
  mergedMorbid$V1=gsub("chr", "", mergedMorbid$V1)
  mergedMorbid.filtered=mergedMorbid[,c(1,4,5,10)]
  ##add 3kbp buffer to merged morbid coordinates as per bionano analysis protocal
  mergedMorbid.filtered$V4=mergedMorbid.filtered$V4-3000
  mergedMorbid.filtered$V5=mergedMorbid.filtered$V5+3000
  ##get files into GRanges objects
  gra=GRanges(seqnames=array_input$Chromosome, ranges=paste0(array_input$Min,"-",array_input$Max), strand=NULL) 
  gr2=GRanges(seqnames=mergedMorbid.filtered$V1, ranges=paste0(mergedMorbid.filtered$V4,"-",mergedMorbid.filtered$V5), strand=NULL)        
  ##overlap BED with mergedMorbid
  hits=findOverlaps(gra, gr2)
  ##compile into one data frame
  comb=cbind(as.data.frame(gra[queryHits(hits)]), variantID=array_input$Call_ID[as.data.frame(hits)[,1]], as.data.frame(gr2[subjectHits(hits)]), gene=mergedMorbid.filtered$gene[as.data.frame(hits)[,2]])
  colnames(comb)=c("sv_chr","sv_start","sv_end","sv_size","sv_strand","sv_ID","gene_chr","gene_start","gene_end","gene_size","gene_strand","gene_name")
  ## map genes to svs
  svs_map2genes=unique(comb$sv_ID)
  array_input$OMIMgenes=NA
  for (i in 1:length(svs_map2genes)) {
    gene.pull=paste0(unique(comb$gene_name[which(comb$sv_ID==svs_map2genes[i])]), collapse=",")
    array_input$OMIMgenes[which(array_input$Call_ID==svs_map2genes[i])]=gene.pull
  }

  #Loop to match Array input to corrosponding bionano input
  for (f in fileListbn){
    #f <- "W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//BN_output_test/MCW_SVI_0001_-_De_novo_pipeline_results.zip"
    bioname<- sub("_-_De_novo_rerun_pipeline_results.zip.*","",f)
    bioname<- sub("_-_De_novo_pipeline_results.zip.*","",bioname)
    bioname <- sub(".*W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//BN_output_test/","",bioname)
    
    #Checks against dummyname to account for difference in naming for corielle sample run on array vs bionano(HM)
    if (dummyname == bioname){
      smap <- read_tsv(archive_read(f,file = "output/contigs/annotation/variants_combine_filters_inMoleRefine1.smap"))
      headstart <- smap %>% pull(var = 1) %>% str_which("SmapEntryID")
      
      
      smap <- read_tsv(archive_read(f,file = "output/contigs/annotation/variants_combine_filters_inMoleRefine1.smap"),skip = headstart)
      smap <- dplyr::rename(smap,"ID"="#h SmapEntryID")
      headers<-colnames(smap)
      smap <- read_tsv(archive_read(f,file = "output/contigs/annotation/variants_combine_filters_inMoleRefine1.smap"),col_names = headers, skip = (headstart+2))
      cnv <- read_tsv(archive_read(f,file = "output/contigs/annotation/CNVs_annotation/cnv_calls_exp_annotation_results.txt"), skip =3)
      report <- read_tsv(archive_read(f,file = "output/exp_informaticsReport.txt"))
      qcresults <-list()
      #What we care about from the quality report
      qclist <- list("Solve Version:","Total number of molecules:","Total length \\(Mbp\\)", "Average length \\(kbp\\)", "Molecule N50 \\(kbp\\)","Label density \\(/100kb\\)","NLV","PLV","bpp","res:","sd:","sf:","sr", "Total number of aligned molecules","Fraction of aligned molecules","Total molecule align length \\(Mbp\\) ","Total reference align length \\(Mbp\\)",
                     "Effective coverage of reference \\(X\\)", "Average aligned length \\(kbp\\) ", "Fraction aligned length", "Average confidence"            
      )
      #create list to get index of info from quality report
      occur <- list()
      for (y in qclist){
        #could create dummy list and than take the 1st value from that append it and have it reset the list each loop
        dummy_list <-list()
        dummy_list<- append(dummy_list, report %>% pull(var = 1) %>% str_which(y))
        #Takes 1st occurance of item from dummy list and appends it to a list so we dont retain locations we do not want
        occur <- append(occur,dummy_list[1])
        # report$`/program/bionano-1.5.1/1.5.1/Pipeline/1.0/pipelineCL.py -t /program/bionano-1.5.1/1.5.1/RefAligner/1.0 -b RawMolecules.bnx -a optArguments_haplotype_DLE1_saphyr_human.xml -d -U -y -i 5 -F 1 -W 1 -c 1 -f 0.1 -J 64 -jp 64 -N 4 -l output -r hg38_DLE1_0kb_0labels.cmap -G hg38_DLE1_gap_common_segdup_min2_com10kb_seg50kb.bed --vapini annotation.ini -json {"toolsVer":"1.5.1","solveVer":"Solve3.5.1_01142020"} -C /program/bionano-1.5.1/1.5.1/Pipeline/1.0/clusterArguments_xeon.xml -T 144`[9] 
      }
      #Creates List to get the actual values based on locations gotten from loop
      qc_values <- list()
      for (z in occur){
        qc_values <- append(qc_values, report %>% pluck(1,z))
      }
      
      #Gets values out of list of list
      sv <- unlist(qc_values)
      #Loop to string the information before : from the list
      i=1
      for (z in sv){
        sv[i] <- sub(".*:","",z)
        sv[i] = str_trim(sv[i])
        print(sv[i])
        i=i+1
      }
      sv
      class(as.numeric(sv[3]))
      sv[3] = as.numeric(sv[3])
      str(sv[3])
      sv[1]
      class(sv[3])
      
     
      
      #Filter and make BED
    
      colnames(cnv)[1]="Id"
      
      
      #Manipulate SMAP file
      smap$ID=paste0("SMAP", smap$ID)
      smap$OverlapGenes[which(smap$OverlapGenes=="-")]=FALSE
      smap$OverlapGenes[which(smap$OverlapGenes!=FALSE)]=TRUE
      colnames(smap)
      plot(smap$SVsize, smap$Size) ##Size vs SVsize
      
      if(length(smap) > 48){
        
      smap=smap[,c(1,3,7,8,9,10,30,25,32,33,44,46,42,48,49)]    ##SMAP format for training/coriell samples
      } else {
      smap=smap[,c(1,3,7,8,9,10,30,25,32,33,38,40,36,42,43)]     ##SMAP format for 'newer' samples
      }
      
      #Manipulate CNV file
      cnv$Id=paste0("CNV", cnv$Id)
      cnv$OverlapGenes[which(cnv$OverlapGenes=="-")]=FALSE
      cnv$OverlapGenes[which(cnv$OverlapGenes!=FALSE)]=TRUE
      cnv$Present_in_BNG_control_samples=NA
      cnv$Present_in_BNG_control_samples_w_same_enzyme=NA
      cnv$Fail_chimeric_score=NA
      cnv$Self_molecules=NA
      cnv$Self_molecules_count=NA
      colnames(cnv)
      cnv=cnv[,c(1:4,9,6,12,5,19,20,13,15,21:23)]
      
      #Combine SMAP and CNVs
      colnames(smap)=colnames(cnv)
      comb=as.data.frame(rbind(smap, cnv))
      comb=comb[order(comb$Chromosome, comb$Start),]
      
      #Add column for masked variants
      ##**## masking only tracked for SVs in CNV file, not SMAP file ##**##
      comb$Masked=FALSE
      comb$Masked[grep("masked", comb$Type)]=TRUE
      
      #Recode labels for SV types
      comb$Type=as.character(comb$Type)
      table(comb$Type)
      comb$Type_general=comb$Type
      comb$Type_general[c(grep("deletion", comb$Type), grep("loss", comb$Type))]="<DEL>"
      comb$Type_general[c(grep("duplication", comb$Type), grep("gain", comb$Type))]="<DUP>"
      comb$Type_general[grep("insertion", comb$Type)]="<INS>"
      comb$Type_general[grep("inversion", comb$Type)]="<INV>"
      comb$Type_general[grep("trans", comb$Type)]="<BND>"
      table(comb$Type_general)
      
      
      
      
      #Remove translocations for testing
      translocations<- filter(comb,(grepl('<BND>',Type_general)))
      translocations <- translocations %>%
        filter(Present_in_BNG_control_samples<=1) %>%
        filter(Present_in_BNG_control_samples_w_same_enzyme<=1) %>%
        filter(Confidence>=.05) %>%
        filter(Self_molecules_count >= 5)
        
      
      
      #Filters out types not feasable to intersection
      #Code might be redundant, check out and fix
      bed1 <- filter(comb,!(grepl('<BND>',Type_general)))
      bed1 <- filter(bed1,!(grepl('inversion_partial',Type)))
      #change formating to bed  
      
      ##INTERSECT WITH MERGED MORBID
      ##read in merged morbid
      mergedMorbid=read.csv("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/vcf/Corriell_vcf/bedfiles/bionano_annotation_files/MERGED_morbid_GRCh38_GTF.gtf", sep="\t", header=F)
      mergedMorbid$gene=gsub("gene_id ", "", sapply(strsplit(mergedMorbid$V9, "; "), `[`, 1))
      mergedMorbid$V1=gsub("chr", "", mergedMorbid$V1)
      mergedMorbid.filtered=mergedMorbid[,c(1,4,5,10)]
      ##add 3kbp buffer to merged morbid coords
      mergedMorbid.filtered$V4=mergedMorbid.filtered$V4-3000
      mergedMorbid.filtered$V5=mergedMorbid.filtered$V5+3000
      ##get files into GRanges objects
      gr1=GRanges(seqnames=bed1$Chromosome, ranges=paste0(bed1$Start,"-",bed1$End), strand=NULL) 
      gr2=GRanges(seqnames=mergedMorbid.filtered$V1, ranges=paste0(mergedMorbid.filtered$V4,"-",mergedMorbid.filtered$V5), strand=NULL)        
      ##overlap BED with mergedMorbid
      hits=findOverlaps(gr1, gr2)
      hits
      ##compile into one data frame
      comb=cbind(as.data.frame(gr1[queryHits(hits)]), variantID=bed1$Id[as.data.frame(hits)[,1]], as.data.frame(gr2[subjectHits(hits)]), gene=mergedMorbid.filtered$gene[as.data.frame(hits)[,2]])
      colnames(comb)=c("sv_chr","sv_start","sv_end","sv_size","sv_strand","sv_ID","gene_chr","gene_start","gene_end","gene_size","gene_strand","gene_name")
      ## map genes to svs
      svs_map2genes=unique(comb$sv_ID)
      bed1$OMIMgenes=NA
      for (i in 1:length(svs_map2genes)) {
        gene.pull=paste0(unique(comb$gene_name[which(comb$sv_ID==svs_map2genes[i])]), collapse=",")
        bed1$OMIMgenes[which(bed1$Id==svs_map2genes[i])]=gene.pull
      }
      
      
      ##INTERSECT WITH MASK FILE
      mask=read.csv("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/vcf/Corriell_vcf/bedfiles/bionano_annotation_files/hg38_CombinedMask_cnv-nbase-segdup-common.txt", sep="\t", header=F)
      gr3=GRanges(seqnames=mask$V1, ranges=paste0(mask$V2,"-",mask$V3), strand=NULL) 
      gr3.reduced=IRanges::reduce(gr3)
      ##overlap BED with masks
      hits=findOverlaps(gr1, gr3.reduced)
      ##compile into one data frame
      comb=cbind(as.data.frame(gr1[queryHits(hits)]), variantID=bed1$Id[as.data.frame(hits)[,1]], as.data.frame(gr3.reduced[subjectHits(hits)]))
      colnames(comb)=c("sv_chr","sv_start","sv_end","sv_size","sv_strand","sv_ID","mask_chr","mask_start","mask_end","mask_size","mask_strand")
      #filter based on percent overlap
      comb$percentOverlap <- width(pintersect(gr1[queryHits(hits)], gr3.reduced[subjectHits(hits)])) / width(gr1[queryHits(hits)])
      ## map genes to svs
      svs_map2mask=unique(comb$sv_ID)
      bed1$Maskedoverlap=0
      for (i in 1:length(svs_map2mask)) {
        bed1$Maskedoverlap[which(bed1$Id==svs_map2mask[i])]=sum(comb$percentOverlap[which(comb$sv_ID==svs_map2mask[i])])
      }
      
      #View(bed1)
      #Mark SV's with any masked overlap to true
      bed1$Masked[which(bed1$Maskedoverlap != 0)] = TRUE
      
      
      
      ##**## ADD TO JARED SCRIPT ##**################################################################################
      hits=findOverlaps(gra, gr1, select = "all")
      ##**##array "ID" changed to "Call ID"##**##
      comb=cbind(arrayID=array_input$Call_ID[as.data.frame(hits)[,1]], as.data.frame(gra[queryHits(hits)]), 
                 bngID=bed1$Id[as.data.frame(hits)[,2]], as.data.frame(gr1[subjectHits(hits)]))
      colnames(comb)=c("array_ID", "array_chr","array_start","array_end","array_size","array_strand","bng_ID","bng_chr","bng_start","bng_end","bng_size","bng_strand")
      setdiff(array_input$Call_ID, unique(comb$array_ID))
      #filter based on percent overlap
      comb$array_percentOverlapping.bng <- width(pintersect(gra[queryHits(hits)], gr1[subjectHits(hits)])) / width(gra[queryHits(hits)])
      comb$bng_percentOverlapping.array <- width(pintersect(gra[queryHits(hits)], gr1[subjectHits(hits)])) / width(gr1[subjectHits(hits)])
      comb$perc.PASS=F
      comb$perc.PASS[which(comb$array_percentOverlapping.bng>=0.50&comb$bng_percentOverlapping.array>=0.50)]=TRUE
      
      ##**##If there are CMA calls with no overlaps, add to comb dataframe##**##
      cma.unique=setdiff(array_input$Call_ID, unique(comb$array_ID))
      if(length(cma.unique!=0)) {
        cma.unique_df=cbind(as.data.frame(array_input[match(cma.unique, array_input$Call_ID),c(19,1,2,3,9)]), "*")
        colnames(cma.unique_df)=colnames(comb)[1:ncol(cma.unique_df)]
        add.cols=colnames(comb)[7:ncol(comb)]
        cma.unique_df[add.cols]=NA
        comb=as.data.frame(rbind(comb, cma.unique_df))
        comb=comb[order(as.numeric(gsub("CMA","",comb$array_ID))),]
      }
      
      ########################
      ## GENERATE FILE #1 ####
      #########################
      ##output array results with or without corresponding bionano calls
      f1.out=data.frame(Sample_ID=rep(dummyname, nrow(comb)))
      f1.out$CMA_ID=comb$array_ID
      f1.out$CMA_chr=comb$array_chr
      f1.out$CMA_start=comb$array_start
      f1.out$CMA_end=comb$array_end
      f1.out$CMA_size=comb$array_size
      f1.out$CMA_type=array_input$Type[match(f1.out$CMA_ID, array_input$Call_ID)]
      f1.out$CMA_ReportableRange=F; 
      f1.out$CMA_ReportableRange[which(f1.out$CMA_end-f1.out$CMA_start>=200000)]=T
      f1.out$CMA_percent.masked=array_input$Maskedoverlap[match(f1.out$CMA_ID, array_input$Call_ID)]
      f1.out$CMA_genes=array_input$OMIMgenes[match(f1.out$CMA_ID, array_input$Call_ID)]
      for (i in 1:nrow(f1.out)) {
        if(f1.out$CMA_ID[i]%in%cma.unique) {f1.out$CMA_number.overlappingBNG[i]=0; next}
        f1.out$CMA_number.overlappingBNG[i]=length(which(f1.out$CMA_ID==f1.out$CMA_ID[i]))
      }
      f1.out$CMA_percent.overlappingBNG=comb$array_percentOverlapping.bng
      f1.out$BNG_ID=comb$bng_ID
      f1.out$BNG_chr=comb$bng_chr
      f1.out$BNG_start=comb$bng_start
      f1.out$BNG_end=comb$bng_end
      f1.out$BNG_size=comb$bng_size
      f1.out$BNG_type=bed1$Type_general[match(f1.out$BNG_ID, bed1$Id)]
      f1.out$BNG_percent.masked=bed1$Maskedoverlap[match(f1.out$BNG_ID, bed1$Id)]
      f1.out$BNG_genes=bed1$OMIMgenes[match(f1.out$BNG_ID, bed1$Id)]
      f1.out$BNG_percent.overlappingCMA=comb$bng_percentOverlapping.array
      f1.out$cond1_overlap.50percent=FALSE
      f1.out$cond1_overlap.50percent[which(is.na(f1.out$BNG_ID))]=NA
      f1.out$cond1_overlap.50percent[which(f1.out$CMA_percent.overlappingBNG>=0.50&f1.out$BNG_percent.overlappingCMA>=0.50)]=TRUE
      f1.out$cond2_mapped.sameGenes=NA
      for (i in 1:nrow(f1.out)) {
        if(is.na(f1.out$CMA_genes[i])&is.na(f1.out$BNG_genes[i])) {next}
        genes.array=strsplit(f1.out$CMA_genes[i], ",")[[1]]
        genes.bng=strsplit(f1.out$BNG_genes[i], ",")[[1]]
        if(identical(genes.array, genes.bng)) {f1.out$cond2_mapped.sameGenes[i]="YES"; next}
        if(length(intersect(genes.array, genes.bng))!=0) {f1.out$cond2_mapped.sameGenes[i]="PARTIAL"}
        if(length(intersect(genes.array, genes.bng))==0) {f1.out$cond2_mapped.sameGenes[i]="NO"}
      }
      
      
      ## RESULT SUMMARY
      report.summary=data.frame(CMA.calls=rep(c("Reportable Calls", "All Calls"), each=8), 
                                Masked.regions=rep(rep(c("Masked and Unmasked regions","Unmasked Only"), each=4),2),
                                Results=rep(c("Variants","Variants with 50% overlap","Variants with 50% overlap and same OMIM genes","Discordent Variants"),4),
                                Count=NA, CMA_IDs=NA)
      ##indexes for conditions
      idx.cond1_50overlap=which(f1.out$cond1_overlap.50percent==T)
      idx.cond2_genes=which(f1.out$cond2_mapped.sameGenes!="NO"|is.na(f1.out$cond2_mapped.sameGenes))
      ##reportable masked/unmasked 
      idx.cma_reportable=which(f1.out$CMA_ReportableRange==T)
      report.summary$Count[1]=length(unique(f1.out$CMA_ID[idx.cma_reportable]))
      report.summary$CMA_IDs[1]=paste0(unique(f1.out$CMA_ID[idx.cma_reportable]), collapse=",")
      report.summary$Count[2]=length(unique(f1.out$CMA_ID[intersect(idx.cma_reportable, idx.cond1_50overlap)]))
      report.summary$CMA_IDs[2]=paste0(unique(f1.out$CMA_ID[intersect(idx.cma_reportable, idx.cond1_50overlap)]), collapse=",")
      report.summary$Count[3]=length(unique(f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable, idx.cond1_50overlap))]))
      report.summary$CMA_IDs[3]=paste0(unique(f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable, idx.cond1_50overlap))]), collapse=",")
      report.summary$Count[4]=length(setdiff(unique(f1.out$CMA_ID[idx.cma_reportable]), f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable, idx.cond1_50overlap))]))
      report.summary$CMA_IDs[4]=paste0(setdiff(unique(f1.out$CMA_ID[idx.cma_reportable]), f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable, idx.cond1_50overlap))]), collapse=",")
      ##reportable unmasked only
      idx.cma_reportable.unmasked=which(f1.out$CMA_ReportableRange==T&f1.out$CMA_percent.masked==0)
      report.summary$Count[5]=length(unique(f1.out$CMA_ID[idx.cma_reportable.unmasked]))
      report.summary$CMA_IDs[5]=paste0(unique(f1.out$CMA_ID[idx.cma_reportable.unmasked]), collapse=",")
      report.summary$Count[6]=length(unique(f1.out$CMA_ID[intersect(idx.cma_reportable.unmasked, idx.cond1_50overlap)]))
      report.summary$CMA_IDs[6]=paste0(unique(f1.out$CMA_ID[intersect(idx.cma_reportable.unmasked, idx.cond1_50overlap)]), collapse=",")
      report.summary$Count[7]=length(unique(f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable.unmasked, idx.cond1_50overlap))]))
      report.summary$CMA_IDs[7]=paste0(unique(f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable.unmasked, idx.cond1_50overlap))]), collapse=",")
      report.summary$Count[8]=length(setdiff(unique(f1.out$CMA_ID[idx.cma_reportable.unmasked]), f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable.unmasked, idx.cond1_50overlap))]))
      report.summary$CMA_IDs[8]=paste0(setdiff(unique(f1.out$CMA_ID[idx.cma_reportable.unmasked]), f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable.unmasked, idx.cond1_50overlap))]), collapse=",")
      ##all masked/unmasked
      report.summary$Count[9]=length(unique(f1.out$CMA_ID))
      report.summary$CMA_IDs[9]=paste0(unique(f1.out$CMA_ID), collapse=",")
      report.summary$Count[10]=length(unique(f1.out$CMA_ID[idx.cond1_50overlap]))
      report.summary$CMA_IDs[10]=paste0(unique(f1.out$CMA_ID[idx.cond1_50overlap]), collapse=",")
      report.summary$Count[11]=length(unique(f1.out$CMA_ID[intersect(idx.cond2_genes, idx.cond1_50overlap)]))
      report.summary$CMA_IDs[11]=paste0(unique(f1.out$CMA_ID[intersect(idx.cond2_genes,  idx.cond1_50overlap)]), collapse=",")
      report.summary$Count[12]=length(setdiff(unique(f1.out$CMA_ID), f1.out$CMA_ID[intersect(idx.cond2_genes,  idx.cond1_50overlap)]))
      report.summary$CMA_IDs[12]=paste0(setdiff(unique(f1.out$CMA_ID), f1.out$CMA_ID[intersect(idx.cond2_genes, idx.cond1_50overlap)]), collapse=",")
      ##all unmasked only
      idx.cma_unmasked=which(f1.out$CMA_percent.masked==0)
      report.summary$Count[13]=length(unique(f1.out$CMA_ID[idx.cma_unmasked]))
      report.summary$CMA_IDs[13]=paste0(unique(f1.out$CMA_ID[idx.cma_unmasked]), collapse=",")
      report.summary$Count[14]=length(unique(f1.out$CMA_ID[intersect(idx.cma_unmasked, idx.cond1_50overlap)]))
      report.summary$CMA_IDs[14]=paste0(unique(f1.out$CMA_ID[intersect(idx.cma_unmasked, idx.cond1_50overlap)]), collapse=",")
      report.summary$Count[15]=length(unique(f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_unmasked, idx.cond1_50overlap))]))
      report.summary$CMA_IDs[15]=paste0(unique(f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_unmasked, idx.cond1_50overlap))]), collapse=",")
      report.summary$Count[16]=length(setdiff(unique(f1.out$CMA_ID[idx.cma_unmasked]), f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_unmasked, idx.cond1_50overlap))]))
      report.summary$CMA_IDs[16]=paste0(setdiff(unique(f1.out$CMA_ID[idx.cma_unmasked]), f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_unmasked, idx.cond1_50overlap))]), collapse=",")
      
      #Create folders for output
      
      if (!dir.exists(paste0("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//Validation_comparison_output//",arrayname, sep =""))){
        dir.create(paste0("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//Validation_comparison_output//",arrayname, sep =""))
      }
      
      
      
      
      ###################################################################################################################
      #Get Array IDs from intersection
      Array_con <- as_tibble(comb$array_ID)
      #Get all ID's from the array input
      Array_idall <- as_tibble(array_input$Call_ID)
      #Testing code to add guaranteed discordant ID
      #Array_idall <-add_row(Array_idall, value = 12)
      #Array_idall <-add_row(Array_idall, value = 13)
      
      #Get all the ID that are not concordant
      db <- setdiff(Array_idall,Array_con)
      #Test for discordant value that exists
      #db <-add_row(db, value = 5)
      # if (dim(db)[1] != 0){
      #   array_input %>% filter(Call_ID==db$value)
      #   array_unique <-subset(array_input, Call_ID %in% db$value)
      # }
      
      
      
      #Apply filters for bionano
      
      
      #general filters - not type specific
      idx.general=which((bed1$Present_in_BNG_control_samples<=1 | 
                           is.na(bed1$Present_in_BNG_control_samples)) &
                          (bed1$Present_in_BNG_control_samples_w_same_enzyme<=1 | 
                             is.na(bed1$Present_in_BNG_control_samples_w_same_enzyme)) &
                          (bed1$Self_molecules_count>5 | is.na(bed1$Self_molecules)) &
                          (bed1$OverlapGenes==TRUE | 
                             bed1$OverlapGenes==FALSE & bed1$NearestNonOverlapGeneDistance<=3000) &
                          (bed1$Fail_chimeric_score!="fail" | is.na(bed1$Fail_chimeric_score)) &
                          (bed1$Self_molecules!="no" | is.na(bed1$Self_molecules)))
      idx.general
      
      #type-specific filters - length and qual
      idx.spec=which((bed1$Type_general=="<INS>" & bed1$Size>1500 & bed1$Confidence>0) |
                       (bed1$Type_general=="<INV>" & bed1$Size>0 & bed1$Confidence>0.70) |
                       (bed1$Type_general=="<DEL>" & bed1$Algorithm=="assembly_comparison" & 
                          bed1$Size>1500 & bed1$Confidence>0) |
                       (bed1$Type_general=="<DEL>" & bed1$Algorithm=="Region-based" & 
                          bed1$Size>200000 & bed1$Confidence>0.99) | 
                       (bed1$Type_general=="<DUP>" & bed1$Algorithm=="assembly_comparison" & 
                          bed1$Size>0 & bed1$Confidence==(-1)) |
                       (bed1$Type_general=="<DUP>" & bed1$Algorithm=="Region-based" & 
                          bed1$Size>200000 & bed1$Confidence>0.99) |
                       (bed1$Type_general=="<BND>" & bed1$Confidence>=0.05))
      idx.spec
      
      ## Filter for final set of SVs
      filtered.svs=bed1[intersect(idx.general, idx.spec),]
      #View(filtered.svs)
      #Grab IDs form filtered bionano samples
      bionanoIds <- as_tibble(filtered.svs$Id)
      #Grab ID's from all bionano samples in the concordant combined tabled
      concordantIDs <- as_tibble(f1.out$BNG_ID)
      #Any samples from the concordant that is not found in the filtered bionano means it failed a QC metric, gotten by finding difference of sets
      confail <- setdiff(concordantIDs,bionanoIds)
      #Create a column for bionano qc and set it to true
      f1.out$bnqc = TRUE
      #Change true to false if the sample ID was in the tibble of concordant samples that failed a QC metric
      f1.out$bnqc[which(f1.out$BNG_ID %in% confail$value)] = FALSE
      ## length differences between start-end coordinates vs SV length column, why? compare to VCF coordinates?
      #View(data.frame(cbind(length.coords=(filtered.svs$End-filtered.svs$Start), length.sv=filtered.svs$Size)))
      
      ##Check if bedout is used, and delete code
      bedout <- filtered.svs
      #bedout <- rename(bedout,"#Chromosome"=Chromosome)
      bedout
  
      
      failed_bncall <- bed1 %>% filter(Id %in% confail$value)
      #Relocate Call_Id in array
      array_input <- array_input %>% dplyr::relocate(Call_ID, .before = Chromosome)
      
      
      ## CMA CONCORDENT/DISCORDENT CALLS
      calls.concord=unique(f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable, idx.cond1_50overlap))])
      df.concordance=data.frame()
      if (is_empty(calls.concord) != TRUE){
        df.concordance=data.frame(Concordant="PASS", as.data.frame(array_input[match(calls.concord, array_input$Call_ID),c(1, 2,3,4,6,7,8,9,10)]), Corresponding_BNG.ID="NA")
        for (i in 1:nrow(df.concordance)) {
          f1.out_cma2bng=f1.out[which(f1.out$CMA_ID%in%df.concordance$Call_ID[i]),]  
          idx.cma2bng=which(f1.out_cma2bng$cond1_overlap.50percent==T&(f1.out_cma2bng$cond2_mapped.sameGenes!="NO"|is.na(f1.out_cma2bng$cond2_mapped.sameGenes)))
          df.concordance$Corresponding_BNG.ID[i]=paste0(f1.out_cma2bng$BNG_ID[idx.cma2bng], collapse=";")
      }
        calls.discord=setdiff(unique(f1.out$CMA_ID[idx.cma_reportable]), f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable, idx.cond1_50overlap))])
        if (is_empty(calls.discord) != TRUE){
          df.concordance=as.data.frame(rbind(df.concordance, data.frame(Concordant="FAIL", as.data.frame(array_input[match(calls.discord, array_input$Call_ID),c(1, 2,3,4,6,7,8,9,10)]), Corresponding_BNG.ID="NA")))
        }
          df.concordance$Notes_from_array=NA
        df.concordance$Notes_from_bionano=NA
      } else{
        calls.discord=setdiff(unique(f1.out$CMA_ID[idx.cma_reportable]), f1.out$CMA_ID[intersect(idx.cond2_genes, intersect(idx.cma_reportable, idx.cond1_50overlap))])
        if (is_empty(calls.discord) != TRUE){
        df.concordance=as.data.frame(rbind(df.concordance, data.frame(Concordant="FAIL", as.data.frame(array_input[match(calls.discord, array_input$Call_ID),c(1, 2,3,4,6,7,8,9,10)]), Corresponding_BNG.ID="NA")))
        }
        df.concordance$Notes_from_array=NA
        df.concordance$Notes_from_bionano=NA
      }
      
      ##write output to excel
      library(openxlsx)
      fname.out=paste0("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Bionano output//Validation_comparison_output//",arrayname,"/report_summary-", arrayname, ".xlsx")
      list_toWrite=list(report.summary,df.concordance, f1.out, array_input,filtered.svs,failed_bncall,translocations)
      write.xlsx(list_toWrite, fname.out, sheetName=c("summary_report","reportable_overlap_results","full_overlap_results","array_input","bionano_filtered","bionano_failedqc","translocations"))
    
      
      
      
      
     #########Fix Mastersheet to be adujsted to new formating######
      ##Add array calls to mastersheet
      array_input
      array_input <-array_input %>%
        mutate(`Segment Interpretation`= as.character(`Segment Interpretation`)) %>%
        mutate(Sample_ID=arrayname)
      array_input <-array_input %>%
        relocate(Sample_ID,.before = Call_ID)
      Arraycall_Mastersheet <- read_csv("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//Arraycall_Mastersheet.csv",col_types = "cccddcdcccdddcdcdcccdcdccdc")
      Arraycall_Mastersheet <- Arraycall_Mastersheet %>% full_join(array_input)
      Arraycall_Mastersheet <- Arraycall_Mastersheet[!duplicated(Arraycall_Mastersheet), ]
      write_csv(Arraycall_Mastersheet,"W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//Arraycall_Mastersheet.csv")
      
      
      #Reads in Mastersheet and adds values from the samples QC to the QC Mastersheet
      QC_Mastersheet <- read_csv("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//QC_Mastersheet.csv", col_types ="ccddddddddddddddddddddd" )
      
      QC_Mastersheet <- add_row(QC_Mastersheet,Sample_name = bioname, Solve_Version = sv[1],	Total_number_of_molecules = as.numeric(sv[2]),	Total_length_Mbp = as.numeric(sv[3]),	Average_length_kbp =  as.numeric(sv[4]),	Molecule_N50_kbp =  as.numeric(sv[5]),	"Label_density/100kb" =  as.numeric(sv[6])	,NLV =  as.numeric(sv[7]),	PLV=  as.numeric(sv[8]),	bpp =  as.numeric(sv[9]),	res =  as.numeric(sv[10]),	sd =  as.numeric(sv[11]),	sf=  as.numeric(sv[12]),	sr=  as.numeric(sv[13]),	Total_number_of_aligned_molecules =  as.numeric(sv[14]),	Fraction_of_aligned_molecules= as.numeric(sv[15]),	Total_molecule_align_length_Mbp=  as.numeric(sv[16]),	Total_reference_alinged_length_Mbp=  as.numeric(sv[17]),	Effective_coverage_of_reference_X =  as.numeric(sv[18]),	Average_aligned_length_kbp=  as.numeric(sv[19]),	Fraction_aligned_length=  as.numeric(sv[20]),	Average_confidence= as.numeric(sv[21]))
      QC_Mastersheet <- QC_Mastersheet[!duplicated(QC_Mastersheet), ]
      write_csv(QC_Mastersheet,"W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//QC_Mastersheet.csv")
      
      
     ##Add Bionano filtered calls to mastersheet
      filtered.svs <- filtered.svs %>%
        mutate(Sample_ID=arrayname)
      filtered.svs <- filtered.svs %>% relocate(Sample_ID, .before = Id)
      #Read in mastersheet and add calls to it
      filtered.svs
      filtered.svs <- filtered.svs %>% mutate(OverlapGenes = as.logical(OverlapGenes))
      filtered.svs <- filtered.svs %>% mutate(Fail_chimeric_score = as.character(Fail_chimeric_score))
      
      Bionanocall_Mastersheet <- read_csv("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//Bionanocall_Mastersheet.csv",col_types = "ccidddccdddldccilccd" )
      Bionanocall_Mastersheet <- Bionanocall_Mastersheet %>% full_join(filtered.svs)
      Bionanocall_Mastersheet <- Bionanocall_Mastersheet[!duplicated(Bionanocall_Mastersheet), ]
      write_csv(Bionanocall_Mastersheet,"W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//Bionanocall_Mastersheet.csv")
    
      #Translocation mastersheet
      translocations <- translocations %>%
        mutate(Sample_ID=arrayname)
      translocations <- translocations %>% mutate(OverlapGenes = as.logical(OverlapGenes))
      translocations <- translocations %>% mutate(Fail_chimeric_score = as.character(Fail_chimeric_score))
      translocations <- translocations %>% relocate(Sample_ID, .before = Id)
      translocations_Mastersheet <- read_csv("W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//Translocation_Mastersheet.csv",col_types = "ccidddccdddldccilc" )
      translocations_Mastersheet <- translocations_Mastersheet %>% full_join(translocations)
      translocations_Mastersheet <- translocations_Mastersheet[!duplicated(translocations_Mastersheet), ]
      write_csv(translocations_Mastersheet,"W://Clinical Working//SVI Project//Bionano//Bionano Working Group//Mastersheets//Translocation_Mastersheet.csv")
      
      warnings()
      }
    
  }
    
  }



