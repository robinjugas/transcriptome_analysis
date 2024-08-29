# source("http://bioconductor.org/biocLite.R")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.11")
# BiocManager::install("Rsubread")
# BiocManager::install("edgeR")
# BiocManager::install("M3C")
# BiocManager::install("umap")
# # biocLite() is the bioconductor installer function. Run it without any
# # arguments to install the core packages or update any installed packages. This
# # requires internet connectivity and will take some time!

rm(list = ls())  # vymazu Envrironment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(Rsubread)
library(edgeR)
library(ggplot2)
library(scatterpie)
library(RColorBrewer)
library(sqldf)
library(gplots)
############## export report
library(dplyr)
library(ggplot2)
library(magrittr)
library(knitr)
library(xtable)
library(flextable)
library(officer)
library(expss)
library(openxlsx)

################################################################################
RNAseqDATADIR <- "fastq/"
RSUBREAD_INDEX_PATH <- "/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/ref_data_GRCh37/" # referencni index GRCh37.p13.genome.fa
RSUBREAD_INDEX_BASE <-"GRCh37"#  "hg19"
outdir="output/"
sampath="SAMFILES"


file.names <- dir(RNAseqDATADIR, pattern =".fastq.gz$",recursive = FALSE)

samples <- data.table::fread("samples_RNAseq.tsv", sep = "\t")
cores = 16

j=1


for(j in 1:nrow(samples)){
  
  # sampleName <- strsplit(file.names[j], "[./]")[[1]][1]
  sampleName <- samples$Sample[j]
  
  #define folders and filenamaes
  outdirpac <- file.path(outdir,sampleName) #output folder per samples
  pacfile <- paste(sampleName,"_",RSUBREAD_INDEX_BASE,".bam",sep="") #name of BAM FILE
  inputfile1 <- file.path(RNAseqDATADIR,paste(sampleName,".fastq.gz",sep=""))
  if(!dir.exists(outdirpac)){dir.create(outdirpac,recursive=T)} #create the folder
  bamfile <- file.path(outdirpac,pacfile)
  
  # # ALIGNMENT
  Rsubread::align(index=file.path(RSUBREAD_INDEX_PATH,RSUBREAD_INDEX_BASE), readfile1=inputfile1, nthreads = cores, output_file=bamfile, output_format="BAM")
  # 
  # # FEATURE COUNTS AND ADD GENE NAMES
  detags.IDs=read.table(file="/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/suppdata/annomart.csv",sep=";",header=T)
  outputsamfile =file.path(outdirpac,pacfile)
  annotgtf=list.files(RSUBREAD_INDEX_PATH)[grep("gtf.gz",list.files(RSUBREAD_INDEX_PATH),fixed=T)]
  mycounts<-featureCounts(outputsamfile, annot.ext=file.path(RSUBREAD_INDEX_PATH,annotgtf), isGTFAnnotationFile=TRUE, isPairedEnd=F,nthreads = cores)

  tab=data.frame(mycounts$counts,substr(rownames(mycounts$counts),1,15))
  names(tab)[1]="counts"
  names(tab)[2]="ensembl"
  tab=merge(tab,detags.IDs,by.x="ensembl",by.y="ensembl_gene_id",all.x=T)

  # 
  # #ulozit pacientuv soubor
  filenam=paste(substr(pacfile,1,nchar(pacfile)-4),"_counts.csv",sep="")
  write.table(tab,file=file.path(outdirpac,filenam),sep=";",row.names=F)
  
  #nacteni souboru READ COUNTS
  stdfsubj=read.table(file=file.path(outdirpac,filenam),sep=";",header=T)
  
  # GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.txt
  #tabGTEx_v7_Annotations_SampleAttributesDS.txt
  
  #nacteni souboru GTEX
  stdfbcksel=read.table(file="/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/suppdata/sampleselbalanced498.csv",sep=";",header=T)
  stdfbcksel$ensembl=substr(stdfbcksel$Name,1,15)
  
  #MERGE SAMPLE AND GTEX, remove duplicated ENSG and GENENAMES
  m=merge(stdfsubj[stdfsubj$hgnc_symbol!="",c(1,3,2)],stdfbcksel[,3:ncol(stdfbcksel)],by.x="ensembl",by.y="ensembl",all.x=T)
  names(m)[names(m)=="counts"]=sampleName
  mx=m[,c("ensembl","hgnc_symbol")]
  dup=(sqldf("SELECT ensembl,hgnc_symbol,count(ensembl) as n FROM mx GROUP BY ensembl HAVING count(ensembl)>1"))[,1]
  m=m[!(m$ensembl %in% dup),]
  mx=m[,c("ensembl","hgnc_symbol")]
  dup=(sqldf("SELECT ensembl,hgnc_symbol,count(hgnc_symbol) as n FROM mx GROUP BY hgnc_symbol HAVING count(hgnc_symbol)>1"))[,2]
  m=m[!(m$hgnc_symbol %in% dup) & !is.na(m$hgnc_symbol),]
  rownames(m)=paste(m$hgnc_symbol)
  m=as.matrix(m[,3:ncol(m)])
  
  #filtrace
  mm=m[rowSums(m>10)>50 & rowSums(is.na(m))==0 & m[,1]>3,] # necham radky kde je vic jak 50 sloupcu s hodnotou > 10 AND zadne NA values AND SAMPLE ma READ COUNTS >3 
  dl=DGEList(counts=mm,group=1:ncol(mm))
  
  #normalizace na TMM
  dly=calcNormFactors(dl,method="TMM")
  
  cp=cpm(dly,log=T)
  
  # which(rownames(mm)=="ALK")
  
  dfcp=data.frame(cp)#,c(1,501:1000)
  colnames(dfcp)[1] <- sampleName
  dfcp$m=apply(dfcp[,2:ncol(dfcp)],1,mean)
  dfcp$sd=apply(dfcp[,2:ncol(dfcp)],1,sd)
  dfcp$fc=dfcp[,1]-dfcp$m
  
  # T-test
  t.result <- apply(dfcp, 1, function (x) t.test(x[c(2:(length(x)-3))],mu=x[1]))
  dfcp$p_value <- unlist(lapply(t.result, function(x) x$p.value))
  
  dfcp$abs=dfcp[,1]
  dfcp$gene=rownames(dfcp)
  
  #korekce nekonecna
  mp=min(dfcp$p_value[dfcp$p_value!=0])
  dfcp$p_value[dfcp$p_value==0]=mp
  dfcp$padj=p.adjust(dfcp$p_value,method="fdr")
  
  # Pathway aNotace ke genum
  misc=read.table("/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/suppdata/Misc.txt",sep="\t",header=T)
  pth=read.table("/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/suppdata/MyCancerGenomeNEW.txt",sep="\t",header=T)
  
  dcp=merge(dfcp,pth,by.x="gene",by.y="gene",all.x=T)
  
  numPathways <- length(unique(pth$pathway))
  
  write.table(dcp[,c(1:2,501:508)],file=file.path(outdirpac,paste(sampleName,"_bckgCTR_",ncol(dcp)-9,".csv",sep="")),sep=";",row.names=F)
  
  
  fileConn<-file(file.path(outdirpac,paste(sampleName,"_protokol_bckset.txt",sep="")))
  writeLines(names(dcp)[3:500],fileConn)
  close(fileConn)
  
  #export pro AIOMIC
  ai=data.frame(paste(dcp[,c("gene")]),paste(dcp[,sampleName]))
  names(ai)=c("Name","Expression")
  write.table(ai[!(ai$Name %in% c("","NA")) & ai$Expression!="NA",],file=file.path(outdirpac,paste("AIOMICtxt_",sampleName,".txt",sep="")),sep="\t",row.names=F)
  
  
  
  ####################################################################################################################
  # PDF VOLCANO PLOT
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  pdf(file.path(outdirpac,paste(sampleName,"_bckg",ncol(dfcp)-7,".pdf",sep="")),width = 11, height = 8 )
  print(ggplot(dcp, aes(x = fc, y = -log10(p_value), size = abs,label= gene ))+ geom_point(data= dcp[dcp$gene!="",],color="gray97") +
    geom_point(data= dcp[!is.na(dcp$pathway),],aes(color=pathway))+
    geom_point(data= dcp[ (dcp$gene %in% misc$gene) ,],color="blue",shape=1)+  scale_radius() +
    geom_text(data=dcp[(dcp$gene %in% misc$gene)| !is.na(dcp$pathway),],color="black",size=2)+
    geom_vline(xintercept=c(min(dcp$fc[(ecdf(dcp$fc)(dcp$fc))>0.95]),min(dcp$fc[(ecdf(dcp$fc)(dcp$fc))>0.975]),max(dcp$fc[(ecdf(dcp$fc)(dcp$fc))<0.05]),max(dcp$fc[(ecdf(dcp$fc)(dcp$fc))<0.025])))+
    scale_colour_manual(values=getPalette(numPathways))+theme(legend.text=element_text(size=5),legend.title=element_text(size=8))+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_line(color="gray",linewidth=0.1),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  )
  dev.off()
  Sys.sleep(10)
  
  # PERCENTILY
  pocn="499"
  dcp$perc=apply(dcp[,],1,function(x) (ecdf(as.numeric(x[2:500]))(as.numeric(x[2]))))
  dcp$fcperc=ecdf(dcp$fc)(dcp$fc)
  dcp$fcscale=ifelse(dcp$fcper>0.975,"++++++",ifelse(dcp$fcper>0.95,"+++++",ifelse(dcp$fcper>0.9,"++++",ifelse(dcp$fcper>0.85,"+++",ifelse(dcp$fcper>0.8,"++",ifelse(dcp$fcper>0.75,"+",""))))))
  
  dcp$percinv=ifelse(dcp$perc<0.5,1-dcp$perc,dcp$perc)
  
  write.table(dcp[,c(1:2,501:511)],file=file.path(outdirpac,paste(sampleName,"_bckgCTR_percentilyFC_",ncol(dcp)-13,".csv",sep="")),sep=";",row.names=F)
  
  fcref=min(dcp$fc[dcp$fcperc>0.95])
  pref=min(dcp$p_value[dcp$perc<0.05])
  
  
  #odlozeni souboru
  # file.rename(file.path(sampath,paste(sampleName,".fastq.gz",sep="")),file.path(sampath,paste("done_",sampleName,".fastq.gz",sep="")))
  # file.rename(bamfile,file.path("Output/BAM files",pacfile))
  
  ####################################################################################################################
  # WORD DOCX report
  
  dft=read.table(file=file.path(outdirpac,paste(sampleName,"_bckgCTR_percentilyFC_",pocn,".csv",sep="")),sep=";",header=T)
  
  dft=dft[!is.na(dft$pathway),]
  
  # unique(dft$pathway)
  
  wdft=dft[,c("pathway", "gene","fcscale","fc")]
  wdft=sqldf("SELECT DISTINCT * FROM wdft")
  wdft$fc=round(wdft$fc,2)
  
  names(wdft)=c("Pathway","Gene","Scale","FC")
  wdft$text=paste(wdft$Gene,wdft$FC,wdft$Scale)
  
  wdft$Pathway <- gsub("/","-", wdft$Pathway)
  
  ft=flextable(wdft[,1:4])
  ft <- bold(ft, i = ~ Scale%in%c("+++","++++","+++++","++++++") , 
             j = ~ Gene+Scale+FC)
  
  
  w=split(wdft[,2:4],wdft$Pathway)
  ft=vector(mode = "list", length = length(w))
  for(i in 1:length(w)){
    ft[[i]]=flextable(data.frame(w[[i]]))%>%bold(i = ~ Scale%in%c("+++","++++","+++++","++++++") , 
                                                 j = ~ Gene+Scale+FC)%>%align(j=2:3,part="all",align="center")%>%
      padding(padding=0.1)%>%fontsize(part="all",size=8)%>%fontsize(size=6)%>%
      autofit()
  }
  
  names(ft)=names(w)
  
  doc=read_docx()%>%
    body_add_fpar(value = fpar(ftext(paste("Report for ",substr(filenam,1,nchar(filenam)-11)," created on ",format(Sys.Date(),format="%d.%m.%Y"),sep=""),prop=fp_text(font.family="Arial"))))%>%
    body_add_par(value="")%>%
    body_end_section_continuous()
  
  big_b <- fp_border(color="black", width = 1)
  for(i in 1:length(w)){
    body_add_fpar(doc,value=fpar(ftext(names(ft[i]),prop=fp_text(font.size=8)),fp_p = fp_par(text.align = "center")))
    body_add_par(doc,value="")
    body_add_flextable(doc,value=hline_top(hline_bottom(border_remove(ft[[i]]),border = big_b, part = "body"),border = big_b, part = "body"))
    body_add_par(doc,value="")
  }
  doc=body_end_section_columns(doc,widths=c(2,2,2))
  
  print(doc, target=file.path(outdirpac,paste(sampleName,"_report.docx",sep="")))
  Sys.sleep(10)
  
  ####################################################################################################################
  ## EXCEL report
  dft <- dcp
  dft <- dft[!is.na(dft$pathway),]
  unik_pathways <- unique(dft$pathway)
  
  wb <- openxlsx::createWorkbook()
  
  for(pathway in unik_pathways){
    
    wdft <- dft[dft$pathway==pathway,c("pathway", "gene","fcscale","fc")]
    wdft$fc <- round(wdft$fc,2)
    colnames(wdft) <- c("Pathway","Gene","Scale","FC")
    
    sheetname <- pathway 
    if(pathway=="Receptor Tyrosine Kinase/Growth Factor Signaling"){sheetname <- "Receptor TK/GF Signaling"}
    if(pathway=="Mitogen Activated Protein (MAP) Kinase Signaling"){sheetname <- "MAP Kinase Signaling"}
    if(pathway=="Chromatin Remodeling/DNA Methylation"){sheetname <- "Chromatin/DNA Methylation"}
    if(pathway=="Janus Kinase (JAK)/ (STAT) Signaling"){sheetname <- "JAK/STAT Signaling"}
    if(pathway=="Cellular architecture and microenvironment"){sheetname <- "Cellular architecture"}
    if(pathway=="non-WNT/non-SHH medulloblastoma-related markers"){sheetname <- "Medulloblastoma markers"}
    
    sheetname <- gsub("/","-",sheetname)
    
    openxlsx::addWorksheet(wb, sheetName = sheetname)
    openxlsx::writeData(wb, sheet = sheetname, wdft)
  }
  
  openxlsx::saveWorkbook(wb, file=file.path(outdirpac,paste0(sampleName,"_BALANCED_report.xlsx")), overwrite = TRUE, returnValue = TRUE)
  
  Sys.sleep(10)
  
  
  ####################################################################################################################
  ####################################################################################################################
  ####################################################################################################################
  
  #specificka kontrola bckg
  # vyber skupiny pozadi
  
  zaldcp=dcp
  dest=read.table(file="/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/suppdata/tabGTEx_v7_Annotations_SampleAttributesDS.txt",sep="\t",quote="",header=T,fill=T)
  de=dest[,c(1,6:7)]
  nam=colnames(cp)[2:ncol(cp)]
  subn=data.frame(gsub(".","-",nam,fixed=T))
  names(subn)="gtex"
  de=de[de$SAMPID%in%subn$gtex,]
  
  namdes=sqldf("SELECT * FROM subn LEFT JOIN de ON subn.gtex=de.SAMPID") #správně seřazeno
  
  ##############################################################################
  dataset <- samples[[which(samples$Sample == sampleName, arr.ind = FALSE),2]]
  if(dataset=="no-brain"){
    #mimo mozky
    grnam <- "no_Brain"
    nogrnam <- "Brain"
    sels <- namdes$SMTS!=nogrnam
  }else{
    grnam <- dataset
    sels <- namdes$SMTS==grnam
  }
  
  #nebo jednotlive skupiny
  # grnam="Brain"
  # grnam="Kidney"
  # grnam="Liver"
  # grnam="Blood Vessel"
  # 
  
  
  
  ##############################################################################
  table(sels)
  
  ################ spec
  stdfsubj=read.table(file=file.path(outdirpac,filenam),sep=";",header=T)
  
  stdfbcksel=read.table(file="/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/suppdata/sampleselbalanced498.csv",sep=";",header=T)
  stdfbcksel$ensembl=substr(stdfbcksel$Name,1,15)
  
  m=merge(stdfsubj[stdfsubj$hgnc_symbol!="",c(1,3,2)],stdfbcksel[,3:ncol(stdfbcksel)],by.x="ensembl",by.y="ensembl",all.x=T)
  
  names(m)[names(m)=="counts"]=sampleName
  
  mx=m[,c("ensembl","hgnc_symbol")]
  
  dup=(sqldf("SELECT ensembl,hgnc_symbol,count(ensembl) as n FROM mx GROUP BY ensembl HAVING count(ensembl)>1"))[,1]
  m=m[!(m$ensembl %in% dup),]
  mx=m[,c("ensembl","hgnc_symbol")]
  dup=(sqldf("SELECT ensembl,hgnc_symbol,count(hgnc_symbol) as n FROM mx GROUP BY hgnc_symbol HAVING count(hgnc_symbol)>1"))[,2]
  m=m[!(m$hgnc_symbol %in% dup) & !is.na(m$hgnc_symbol),]
  
  rownames(m)=paste(m$hgnc_symbol)
  
  m=as.matrix(m[,3:ncol(m)]) #   [,c(T,sels)]
  #filtrace
  # pro vetsi komparatory
  mm=m[rowSums(m>10)>50 & rowSums(is.na(m))==0 & m[,1]>3,] 
  # ################
  # pro mensi komparatory
  # mm=m[rowSums(m>3)>10 & rowSums(is.na(m))==0 & m[,1]>3,]
  # ################
  
  #normalizace s celym souborem
  dl=DGEList(counts=mm,group=1:ncol(mm))
  dly=calcNormFactors(dl,method="TMM")
  
  cp=cpm(dly,log=T)
  
  #selekce pro srovnani specificke kontroly
  dfcp=data.frame(cp[,c(T,sels)])
  colnames(dfcp)[1] <- sampleName
  dfcp$m=apply(dfcp[,2:ncol(dfcp)],1,mean)
  dfcp$sd=apply(dfcp[,2:(ncol(dfcp)-1)],1,sd)
  dfcp$fc=dfcp[,1]-dfcp$m
  cnd=dfcp$sd>1e-10
  
  t.result <- apply(dfcp[cnd,], 1, function (x) t.test(x[c(2:(length(x)-3))],mu=x[1]))
  dfcp$p_value[cnd] <- unlist(lapply(t.result, function(x) x$p.value))
  
  dfcp$abs=dfcp[,1]
  dfcp$gene=rownames(dfcp)
  
  mp=min(dfcp$p_value[dfcp$p_value!=0],na.rm=T)
  dfcp$p_value[dfcp$p_value==0&!is.na(dfcp$p_value)]=mp
  dfcp$padj=p.adjust(dfcp$p_value,method="fdr")
  
  
  
  misc=read.table("/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/suppdata/Misc.txt",sep="\t",header=T)
  
  pth=read.table("/home/rj/4TB/CEITEC/TRANSKRIPTOM_DATA/MichalKDO/suppdata/MyCancerGenomeNEW.txt",sep="\t",header=T)
  
  dcp=merge(dfcp,pth,by.x="gene",by.y="gene",all.x=T)
  
  numPathways <- length(unique(pth$pathway))
  
  write.table(dcp[,c(1:2,(ncol(dcp)-7):ncol(dcp))],file=file.path(outdirpac,paste(sampleName,"_specCTR_",grnam,"_",ncol(dcp)-9,".csv",sep="")),sep=";",row.names=F)
  
  
  pdf(file.path(outdirpac,paste(sampleName,"_specBckg_",grnam,ncol(dfcp)-7,".pdf",sep="")),width = 11, height = 8 )
  print(ggplot(dcp, aes(x = fc, y = -log10(p_value), size = abs,label= gene ))+ geom_point(data= dcp[dcp$gene!="",],color="gray97") +
    geom_point(data= dcp[!is.na(dcp$pathway),],aes(color=pathway))+
    geom_point(data= dcp[ (dcp$gene %in% misc$gene) ,],color="blue",shape=1)+  scale_radius() +
    geom_text(data=dcp[(dcp$gene %in% misc$gene)| !is.na(dcp$pathway),],color="black",size=2)+
    scale_colour_manual(values=getPalette(numPathways))+theme(legend.text=element_text(size=5),legend.title=element_text(size=8))+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_line(color="gray",linewidth=0.1),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  )
  dev.off()
  Sys.sleep(10)
  ############ spec end
  
  
  pocn=table(sels)[2]+1
  dcp$perc=apply(dcp[,],1,function(x) (ecdf(as.numeric(x[2:(pocn+1)]))(as.numeric(x[2]))))
  dcp$fcperc=ecdf(dcp$fc)(dcp$fc)
  dcp$fcscale=ifelse(dcp$fcperc>0.975,"++++++",ifelse(dcp$fcperc>0.95,"+++++",ifelse(dcp$fcperc>0.9,"++++",ifelse(dcp$fcperc>0.85,"+++",ifelse(dcp$fcperc>0.8,"++",ifelse(dcp$fcperc>0.75,"+",""))))))
  
  dcp$percinv=ifelse(dcp$perc<0.5,1-dcp$perc,dcp$perc)
  
  write.table(dcp[,c(1:2,(pocn+2):(pocn+12))],file=file.path(outdirpac,paste(sampleName,"_specBckgCTR_",grnam,"_percentilyFC_",ncol(dcp)-13,".csv",sep="")),sep=";",row.names=F)
  
  fcref=min(dcp$fc[dcp$fcperc>0.95])
  pref=min(dcp$p_value[dcp$perc<0.05])
  
  ####################################################################################################################
  # WORD DOCX report
  # report pro specifickou kontrolu
  
  dft=read.table(file=file.path(outdirpac,paste(sampleName,"_specBckgCTR_",grnam,"_percentilyFC_",pocn,".csv",sep="")),sep=";",header=T)
  
  dft=dft[!is.na(dft$pathway),]
  
  #
  wdft=dft[,c("pathway", "gene","fcscale","fc")]
  wdft=sqldf("SELECT DISTINCT * FROM wdft")
  wdft$fc=round(wdft$fc,2)
  
  names(wdft)=c("Pathway","Gene","Scale","FC")
  wdft$text=paste(wdft$Gene,wdft$FC,wdft$Scale)
  
  wdft$Pathway <- gsub("/","-", wdft$Pathway)
  
  ft=flextable(wdft[,1:4])
  ft <- bold(ft, i = ~ Scale%in%c("+++","++++","+++++","++++++") , 
             j = ~ Gene+Scale+FC)
  
  
  w=split(wdft[,2:4],wdft$Pathway)
  ft=vector(mode = "list", length = length(w))
  for(i in 1:length(w)){
    ft[[i]]=flextable(data.frame(w[[i]]))%>%bold(i = ~ Scale%in%c("+++","++++","+++++","++++++") , 
                                                 j = ~ Gene+Scale+FC)%>%align(j=2:3,part="all",align="center")%>%
      padding(padding=0.1)%>%fontsize(part="all",size=8)%>%fontsize(size=6)%>%
      autofit()
  }
  #
  names(ft)=names(w)
  
  doc=read_docx()%>%
    body_add_fpar(value = fpar(ftext(paste("Report for ",substr(filenam,1,nchar(filenam)-11)," created on ",format(Sys.Date(),format="%d.%m.%Y"),sep=""),prop=fp_text(font.family="Arial"))))%>%
    body_add_par(value="")%>%
    body_end_section_continuous()
  
  big_b <- fp_border(color="black", width = 1)
  for(i in 1:length(w)){
    body_add_fpar(doc,value=fpar(ftext(names(ft[i]),prop=fp_text(font.size=8)),fp_p = fp_par(text.align = "center")))
    body_add_par(doc,value="")
    body_add_flextable(doc,value=hline_top(hline_bottom(border_remove(ft[[i]]),border = big_b, part = "body"),border = big_b, part = "body"))
    body_add_par(doc,value="")
  }
  doc=body_end_section_columns(doc,widths=c(2,2,2))
  
  print(doc, target=file.path(outdirpac,paste(sampleName,"_spec_",grnam,"_report.docx",sep="")))
  Sys.sleep(10)
  
  ####################################################################################################################
  ## EXCEL report
  dft <- dcp
  dft <- dft[!is.na(dft$pathway),]
  unik_pathways <- unique(dft$pathway)
  
  wb <- openxlsx::createWorkbook()
  for(pathway in unik_pathways){
    
    wdft <- dft[dft$pathway==pathway,c("pathway", "gene","fcscale","fc")]
    wdft$fc <- round(wdft$fc,2)
    colnames(wdft) <- c("Pathway","Gene","Scale","FC")
    
    sheetname <- pathway 
    if(pathway=="Receptor Tyrosine Kinase/Growth Factor Signaling"){sheetname <- "Receptor TK/GF Signaling"}
    if(pathway=="Mitogen Activated Protein (MAP) Kinase Signaling"){sheetname <- "MAP Kinase Signaling"}
    if(pathway=="Chromatin Remodeling/DNA Methylation"){sheetname <- "Chromatin/DNA Methylation"}
    if(pathway=="Janus Kinase (JAK)/ (STAT) Signaling"){sheetname <- "JAK/STAT Signaling"}
    if(pathway=="Cellular architecture and microenvironment"){sheetname <- "Cellular architecture"}
    if(pathway=="non-WNT/non-SHH medulloblastoma-related markers"){sheetname <- "Medulloblastoma markers"}
    
    sheetname <- gsub("/","-",sheetname)
    
    openxlsx::addWorksheet(wb, sheetName = sheetname)
    openxlsx::writeData(wb, sheet = sheetname, wdft)
  }
  
  openxlsx::saveWorkbook(wb, file=file.path(outdirpac,paste0(sampleName,"_",dataset,"_report.xlsx")), overwrite = TRUE, returnValue = TRUE)
  
  
}

################################################################################
# KONEC



