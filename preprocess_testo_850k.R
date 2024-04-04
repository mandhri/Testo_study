
suppressPackageStartupMessages({
  library("plyr")
  library("R.utils")
  library("missMethyl")
  library("limma")
  library("DMRcate")
  library("DMRcatedata")
  library("topconfects")
  library("minfi")
  library("IlluminaHumanMethylation450kmanifest")
  library("RColorBrewer")
  library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  library("GEOquery")
  library("eulerr")
  library("plyr")
  library("gplots")
  library("reshape2")
  library("forestplot")
  library("beeswarm")
  library("RCircos")
  library("qqman")
  library("ENmix")
})
# Loading libraries
library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(HelpersMG)
library(limma)
library(minfi)
library(readxl)
library(sva)
library(ChAMPdata)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2manifest)
library(DMRcate)
library(DMRcatedata)

#Making a directory for the files 
WORKING_DIR=getwd()
baseDir<-WORKING_DIR
dir.create("deakin")
ARRAY_DATA="deakin.tar"
DEST=paste(WORKING_DIR,"/",ARRAY_DATA,sep="")

list.files(path = "/mnt/vol1/sev_data/IDAT/MET2022-354-014/deakin/")
setwd("/mnt/vol1/sev_data/IDAT/MET2022-354-014/deakin/")
files <- list.files(getwd(),pattern = "20",recursive = TRUE)
files

targets<-read.csv("/mnt/vol1/sev_data/IDAT/MET2022-354-014/deakin/deakin_phenotype.csv")

rgSet <- read.metharray.exp(targets = targets)
rgSet
annotation(rgSet)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(rgSet)["annotation"] = "20a1.hg38"
rgSet

library(IlluminaHumanMethylationEPICv2manifest)
library("IlluminaHumanMethylationEPICv2manifest")
IlluminaHumanMethylationEPICv2manifest
library("IlluminaHumanMethylationEPICv2manifest")
getManifest(rgSet)
manifest <- getManifest(rgSet)
head(getProbeInfo(manifest))
mSet <- preprocessRaw(rgSet)
meth = getMeth(mSet)
unmeth = getUnmeth(mSet)
RSet <- ratioConvert(mSet, what = "both", keepCN = TRUE)
detP <- detectionP(rgSet)
GRset <- mapToGenome(RSet)


# include sex chromosomes
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw <- mSetSw[keep,]
betaValues <- getBeta(mSet)

#
myImport <- champ.import("/mnt/vol1/sev_data/IDAT/MET2022-354-014/deakin/", arraytype = "EPICv2")

Filter <- champ.filter(beta=myImport$beta, 
                         pd=myImport$pd, 
                         detP=myImport$detP,
                         Meth = myImport$Meth,
                         arraytype = "EPICv2") 


#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/home/mandhri/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = Filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/home/mandhri/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = Filter$pd

colnames(pheno)
champ.QC(beta = beta,
         pheno = pheno$Age,
         dendrogram = FALSE)
dev.off()

myNorm <- champ.norm(beta=beta,arraytype = "EPICv2")
which(is.na(myNorm))


champ.SVD(beta=myNorm,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/")

dev.off()

M <- logit2(myNorm)
myCombat=ComBat(dat=as.matrix(M),
                batch=pheno$Sample_Name,
                mod=NULL)

myCombat=ilogit2(myCombat)
champ.QC(beta = myCombat,
         pheno = targets$Age,
         dendrogram = FALSE,
         resultsDir="./CHAMP_QCimages/combat_QCimages/")

dev.off()

write.table(myCombat, "/mnt/vol1/sev_data/IDAT/MET2022-354-014/SARAH_DEAKIN.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")



