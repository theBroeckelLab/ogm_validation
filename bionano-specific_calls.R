

###############################################################################################
## Identification of calls specific to bng (not detected in cma) from comparison reports ######
###############################################################################################


## read in sample report
xl.input="W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Validation_comparison_output/NA15186/report_summary-NA15186.xlsx"
cma=read.xlsx(xl.input, sheet=4)
bng.pass=read.xlsx(xl.input,sheet=5)
bngXcma=read.xlsx(xl.input, sheet=2)


## pull only reportable cma calls
cma.reportable=cma[which(cma$Size>=200000),]


## replace default bng size with the actual difference between coordinates
bng.pass$Size=bng.pass$End-bng.pass$Start


## pull only reportable bng calls
idx.cnv=which((bng.pass$Type_general=="<DEL>"|bng.pass$Type_general=="<DUP>"|bng.pass$Type_general=="<INS>") & bng.pass$Size>=200000)
idx.sv=which(bng.pass$Type_general=="<INV>" & bng.pass$Size>=5000000)
bng.pass.reportable=bng.pass[c(idx.cnv, idx.sv),]


## do we want to only consider non-masked bionano calls?
bng.pass.reportable=bng.pass.reportable[which(bng.pass.reportable$Masked==F),]


## determine if any remaining bng reportable calls are detected in cma
bng.pass.reportable$Id           ##reportable BNG calls
bngXcma$Corresponding_BNG.ID     ##reportable BNG calls detected by CMA
calls_specificToBionano=setdiff(bng.pass.reportable$Id, bngXcma$Corresponding_BNG.ID[which(!is.na(bngXcma$Corresponding_BNG.ID))])
calls_specificToBionano          ##reportable BNG calls NOT detected by CMA


## pull information for unique bionano calls
View(bng.pass.reportable[which(bng.pass.reportable$Id%in%calls_specificToBionano),])
