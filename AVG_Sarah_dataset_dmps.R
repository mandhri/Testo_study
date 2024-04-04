library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(HelpersMG)
library(limma)
library(minfi)
library(readr)

# Creating function called "create_summary" to put the resulted dmps in a summary
create_summary <- function(toptable = NULL,
                           dataset_label = NULL,
                           directory = getwd())
{
  CPG <- rownames(toptable)
  ALLELE1 <- rep(1,nrow(toptable))
  ALLELE2 <- rep(2,nrow(toptable))
  TESTSTAT <- toptable$t
  PVALUE <- toptable$P.Value
  EFFECTSIZE <- toptable$logFC
  SE <- toptable$SE
  results = data.frame(CPG,
                       ALLELE1,
                       ALLELE2,
                       TESTSTAT,
                       PVALUE,
                       EFFECTSIZE,
                       SE)
  write.table(results,
              file=paste0(directory,"/",dataset_label,".tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

B <- data.table::fread("/mnt/vol1/sev_data/IDAT/MET2022-354-014/SARAH_DEAKIN.txt")
k<-data.table::fread("/mnt/vol1/sev_data/IDAT/MET2022-354-014/deakin/")
pheno <- read_csv("/mnt/vol1/sev_data/IDAT/MET2022-354-014/deakin/deakin_phenotype.csv")

#Remove the row.names in the df and name it as row_id
#pheno <- pheno %>%
#tibble::rownames_to_column("row_id")

B <- as.data.frame(B)

rownames(B) <- B$V1

#Drop the column V1 in B.  
#Use dplyr:: since R mistake the select function to MASS package for dplyr
B <- B %>% dplyr::select(-V1)


B[1:3,1:3]

new_rownames <- sub("_.*$", "", rownames(B))
duplicates <- duplicated(new_rownames) | duplicated(new_rownames, fromLast = TRUE)
new_rownames[duplicates] <- paste(new_rownames[duplicates], seq_along(new_rownames[duplicates]), sep = "_")
rownames(B) <- new_rownames

illumina_ID = getManifestInfo(IlluminaHumanMethylationEPICv2manifest, "locusNames")
head(illumina_ID)
any(duplicated(illumina_ID))
probe_ID = gsub("_.*$", "", illumina_ID)
tb = table(probe_ID)
table(tb)
tb[tb == 10]
illumina_ID[probe_ID == "cg06373096"]

library(IlluminaHumanMethylationEPICmanifest)
probe1 = getManifestInfo(IlluminaHumanMethylationEPICmanifest, "locusNames")
probe2 = getManifestInfo(IlluminaHumanMethylationEPICv2manifest, "locusNames")

probe1 = unique(probe1)
probe2 = gsub("_.*$", "", probe2) 
probe2 = unique(probe2)


beta2 = do.call(rbind, tapply(1:nrow(B), gsub("_.*$", "", rownames(B)), function(ind) {
  colMeans(B[ind, , drop = FALSE])
}, simplify = FALSE))

head(beta2)

####
beta2 <- as.data.frame(beta2)


write.table(beta2, "/mnt/vol1/sev_data/IDAT/MET2022-354-014/SARAH_DEAKIN_beta_avg.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#Is there NA values ?
na_rows<-pheno[apply(pheno, 1, function(x) any(is.na(x))),]
colSums(is.na(pheno))

pheno$exercise_status<-str_extract(pheno$Sample_Name,"(?<=-)[A-Z]+")
pheno$exercise_status
pheno$exercise_status_numeric <- ifelse(pheno$exercise_status == "PRE", 1, 2)

#Calculating M values by getting logistic distribution to base 2
M <- logit2(beta2)

#Design matrix
glimpse(pheno)

pheno$participant_id <- gsub("-(PRE|POST)$", "", pheno$Sample_Name)

design=model.matrix(~exercise_status_numeric+ `Hormonal contraceptive use (0=no, 1=yes).x`,
                    pheno)

Msubset <- M[sample(nrow(M),replace=FALSE),]
corfit <- duplicateCorrelation(as.matrix(Msubset), design, block=pheno$participant_id)


corfit$consensus

fit1_M <- lmFit(M,
                design,
                block = pheno$participant_id,
                correlation = 0.2701366)

fit2_M <- eBayes(fit1_M)
Bsubset <- beta2[sample(nrow(beta2),replace=FALSE),]
corfitB <- duplicateCorrelation(Bsubset, #we fit limma models on M values and not beta values
                                design,
                                block=pheno$participant_id)

corfitB$consensus

fit1_B <- lmFit(beta2,
                design,
                block = pheno$participant_id,
                correlation = 0.2933381)
fit2_B <- eBayes(fit1_B)


#Observing top differentially methylated cpgs associated with age for both B values and M values
coef = "exercise_status_numeric"
results <- limma::topTable(fit2_M,
                           coef=coef,
                           number = Inf,
                           p.value = 1)
results_B <- limma::topTable(fit2_B,
                             coef=coef,
                             number=Inf,
                             p.value=1)



#Differential cpg site,
#the log-fold change associated with the comparison we are interested in (DNAM change per unit age)
#the average intensity of the probe set/gene across all chips,
#the (moderated) t-statistic for the hypothesis that the log-fold change is zero (or equivalently, that the fold change is one),
#the associated raw and adjusted p-values for the t-statistic,
##an estimated log-odds ratio for DE.

#Substitute the log FC values of M to the log FC vales of Beta values. In here, log FC values will be normally 
#distributed. 


results$logFC <- results_B[rownames(results),"logFC"]
is.na(results)
na_rows<-results[apply(results, 1, function(x) any(is.na(x))),]
na_rowsq<-results_B[apply(results_B , 1, function(x) any(is.na(x))),]

SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for age DMPs with CTC",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue",bins = 30)
dev.off()

getwd()
fdr = 0.05
results_age=limma::topTable(fit2_M,
                            coef = "exercise_status_numeric",
                            number = nrow(M),
                            adjust.method = "BH",
                            p.value = fdr)

dim(results_age)


directory = ("/home/mandhri/")
create_summary(toptable = results,
               dataset_label = "850K_Sarah",
               directory = directory)

