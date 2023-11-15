## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Additional phenotype file preparation - atrial fibrillation - consent, relatedness and medication 

###############################################
### STEP 1. READ IN AND PREPARE FILES  ########
###############################################

## set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/")

## Read in covariate files 
af=read.csv("Covariates/incident_af.csv")
names(af)[1] <- "ID"

## Remove individuals whose GP surgeries have not consented to linkage analyses 

## Remove incomplete cases (with exception of drug columns as we have lots of missingness here)
## ECG 
n=af[,-c(39:40)]
n=n[complete.cases(n),]
af1=af[which(af$ID %in% n$ID),]


########################################
### STEP 2. RELATEDNESS ANALYSES #######
########################################

###### For each dataframe, random extract only one member from each family and remove first/second degree relatives if still present

## Read in pedigree data 
## Extract FID information 
ped <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/pedigree.csv")
fids <- ped[,c("volid", "famid")]
names(fids) <- c("ID", "FID")

## Merge in FIDs
af1 <- merge(af1, fids, by  ="ID")


## Write out FID and IID for relatedness analysis with plink 
### ECG - with medication ###
tmp1 <- af1[,c("FID","ID")]
write.table(tmp1, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/af.txt", row.names = F, col.names = F, quote = F, sep = "\t")


## In command line: relatedness analysis via plink 
 plink19 --genome --bfile /Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS --keep /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/af.txt --out /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/af


##################################################
### STEP 3. IDENTIFY INDIVIDUALS TO REMOVE #######
##################################################

# read in IBD information and subset for convenience of analysis
rel1=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/af.genome"))
rel1=rel1[rel1$PI_HAT > 0.12,]
# set criterion 
otherCriterion=af1[,c("ID","Case")]
## run relatedness filtering and subset af dataframe accordingly 
x1=relatednessFilter(relatedness=rel1, otherCriterion=otherCriterion, relatednessTh=0.187, otherCriterionTh=0.5, otherCriterionThDirection="lt", otherCriterionIID="ID",otherCriterionMeasure="Case", relatednessIID1 = "IID1", relatednessIID2 = "IID2", relatednessFID1 = "FID1", relatednessFID2 = "FID2")
af2 <- af1[-which(af1$ID %in% x1$failIDs$IID),]



## Separate into a group with medication and a group without 
af_med1 = af2[-which(is.na(af2$Digoxin)),]
af_med1$HTN<- NULL 
af_no_med1 = af2[which(is.na(af2$Digoxin)),]
af_no_med1$Digoxin <- NULL 

#################################################################
### STEP 4. MAKE PHENOTYPE + COVARIATE FILES FOR ANALYSES #######
#################################################################

########## COVARIATES ############


### AF - with medication ###
## Create continuous covariate file
tmp <- af_med1
names(tmp)[1] <- "IID"

quant_cov <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       age = tmp$age, 
                       height = tmp$height,
                       weight = tmp$weight,
                       NK = tmp$NK,
                       Bcell = tmp$Bcell,
                       Mono = tmp$Mono,
                       CD4T = tmp$CD4T,
                       CD8T = tmp$CD8T,
                       SBP = tmp$avg_sys,
 					   PC1 = tmp$PC1,
                       PC2 = tmp$PC2,
                       PC3 = tmp$PC3,
                       PC4 = tmp$PC4,
                       PC5 = tmp$PC5,
                       PC6 = tmp$PC6,
                       PC7 = tmp$PC7,
                       PC8 = tmp$PC8,
                       PC9 = tmp$PC9,
                       PC10 = tmp$PC10,
                       PC11 = tmp$PC11,
                       PC12 = tmp$PC12,
                       PC13 = tmp$PC13,
                       PC14 = tmp$PC14,
                       PC15 = tmp$PC15,
                       PC16 = tmp$PC16,
                       PC17 = tmp$PC17,
                       PC18 = tmp$PC18,
                       PC19 = tmp$PC19,
                       PC20 = tmp$PC20)
write.table(quant_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/AF/quant.qcov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - with drugs included 
fact_cov <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       sex = tmp$sex,
                       batch = tmp$Batch,
                       smoking = tmp$ever_smoke,
                       diabetes = tmp$diabetes_Y,
                       digoxin = tmp$Digoxin)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/AF/fact.cov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - no sex 
fact_cov2 <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       batch = tmp$Batch,
                       smoking = tmp$ever_smoke,
                       diabetes = tmp$diabetes_Y,
                       digoxin = tmp$Digoxin)
write.table(fact_cov2, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/AF/fact_nosex.cov", row.names=F, sep=' ', quote = F)


### ECG - without medication ###
## Create continuous covariate file
tmp <- af_no_med1
names(tmp)[1] <- "IID"

quant_cov <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       age = tmp$age, 
                       height = tmp$height,
                       weight = tmp$weight,
                       NK = tmp$NK,
                       Bcell = tmp$Bcell,
                       Mono = tmp$Mono,
                       CD4T = tmp$CD4T,
                       CD8T = tmp$CD8T,
                       SBP = tmp$avg_sys,
 						PC1 = tmp$PC1,
                       PC2 = tmp$PC2,
                       PC3 = tmp$PC3,
                       PC4 = tmp$PC4,
                       PC5 = tmp$PC5,
                       PC6 = tmp$PC6,
                       PC7 = tmp$PC7,
                       PC8 = tmp$PC8,
                       PC9 = tmp$PC9,
                       PC10 = tmp$PC10,
                       PC11 = tmp$PC11,
                       PC12 = tmp$PC12,
                       PC13 = tmp$PC13,
                       PC14 = tmp$PC14,
                       PC15 = tmp$PC15,
                       PC16 = tmp$PC16,
                       PC17 = tmp$PC17,
                       PC18 = tmp$PC18,
                       PC19 = tmp$PC19,
                       PC20 = tmp$PC20)
write.table(quant_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/AF/quant.qcov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - without drugs included 
fact_cov <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       sex = tmp$sex,
                       batch = tmp$Batch,
                       diabetes = tmp$diabetes_Y,
                       HTN = tmp$HTN,
                       smoking = tmp$ever_smoke)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/AF/fact.cov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - no sex 
fact_cov2 <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       batch = tmp$Batch,
                       diabetes = tmp$diabetes_Y,
                       HTN = tmp$HTN,
                       smoking = tmp$ever_smoke) 
write.table(fact_cov2, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/AF/fact_nosex.cov", row.names=F, sep=' ', quote = F)

