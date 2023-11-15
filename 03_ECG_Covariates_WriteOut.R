## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Additional phenotype file preparation - medication, relatedness and sex 

###############################################
### STEP 1. READ IN AND PREPARE FILES  ########
###############################################

######## ECG TRAITS + RR ###########

## set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/")

## Read in covariate files 
ecg=read.csv("Covariates/ecg_traits.csv")
rr=read.csv("Covariates/rr_interval.csv")

## Remove offending Inf values in RR interval column 
ecg=ecg[-which(is.infinite(ecg$RR_interval)),]
ecg=ecg[-which(ecg$JT==0),]

## Rename phenotypes to what is requested in final files 
names(ecg)[2:5] <- c("QRS", "QT", "PR", "JT")

## Create set column 
ecg$Set <- gsub("_.*", "", ecg$Batch)
rr$Set <- gsub("_.*", "", rr$Batch)

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
### ECG ###
ecg1 <- merge(ecg, fids, by  ="ID")
### RR  ###
rr1 <- merge(rr, fids, by  ="ID")


## Write out FID and IID for relatedness analysis with plink 
### ECG - with medication ###
tmp1 <- ecg1[,c("FID","ID")]
write.table(tmp1, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/ecg.txt", row.names = F, col.names = F, quote = F, sep = "\t")
### ECG - without medication ###
tmp2 <- rr1[,c("FID","ID")]
write.table(tmp2, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/rr.txt", row.names = F, col.names = F, quote = F, sep = "\t")


## In command line: relatedness analysis via plink 
 plink19 --genome --bfile /Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS --keep /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/ecg.txt --out /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/ecg
 plink19 --genome --bfile /Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS --keep /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/rr.txt --out /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/rr


##################################################
### STEP 3. IDENTIFY INDIVIDUALS TO REMOVE #######
##################################################

## Remove incomplete cases (with exception of drug columns as we have lots of missingness here)
## ECG 
n=ecg[,-c(40:43)]
n=n[complete.cases(n),]
ecg1=ecg1[which(ecg1$ID %in% n$ID),]
## RR 
m=rr[,-c(36:39)]
m=m[complete.cases(m),]
rr1=rr1[which(rr1$ID %in% m$ID),] 


### ECG - with medication ###
rel1=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/ecg.genome"))
rel1=rel1[rel1$PI_HAT > 0.12,]
x1=relatednessFilter(relatedness=rel1, relatednessTh=0.187, relatednessIID1 = "IID1", relatednessIID2 = "IID2", relatednessFID1 = "FID1", relatednessFID2 = "FID2")
ecg2 <- ecg1[-which(ecg1$ID %in% x1$failIDs$IID),]
g=rel1[which(rel1$IID1 %in% ecg2$ID & rel1$IID2 %in% ecg2$ID) ,]
x3=relatednessFilter(relatedness=g, relatednessTh=0.187, relatednessIID1 = "IID1", relatednessIID2 = "IID2", relatednessFID1 = "FID1", relatednessFID2 = "FID2")
ecg3 <- ecg2[-which(ecg2$ID %in% x3$failIDs$IID),]


### RR - without medication ###
rel2=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Relatedness/rr.genome"))
rel2=rel2[rel2$PI_HAT > 0.12,]
x2=relatednessFilter(relatedness=rel2, relatednessTh=0.187, relatednessIID1 = "IID1", relatednessIID2 = "IID2", relatednessFID1 = "FID1", relatednessFID2 = "FID2")
rr2 <- rr1[-which(rr1$ID %in% x2$failIDs$IID),]
v=rel2[which(rel2$IID1 %in% rr2$ID & rel2$IID2 %in% rr2$ID) ,]
x4=relatednessFilter(relatedness=v, relatednessTh=0.187, relatednessIID1 = "IID1", relatednessIID2 = "IID2", relatednessFID1 = "FID1", relatednessFID2 = "FID2")
rr3 <- rr2[-which(rr2$ID %in% x4$failIDs$IID),]



## Separate into a group with medication and a group without 
## ECG
ecg_med1 = ecg3[-which(is.na(ecg3$Beta_Agonist)),]
ecg_no_med1 = ecg3[which(is.na(ecg3$Beta_Agonist)),]
## RR 
rr_med1 = rr3[-which(is.na(rr3$Beta_Agonist)),]
rr_no_med1 = rr3[which(is.na(rr3$Beta_Agonist)),]



#################################################################
### STEP 4. MAKE PHENOTYPE + COVARIATE FILES FOR ANALYSES #######
#################################################################

########## COVARIATES ############


### ECG - with medication ###
## Create continuous covariate file
tmp <- ecg_med1
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
                       RR = tmp$RR_interval,
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
write.table(quant_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/ECG/quant.qcov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - with drugs included 
fact_cov <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       sex = tmp$sex,
                       batch = tmp$Set,
                       smoking = tmp$ever_smoke,
                       beta_blocker = tmp$Beta_Blocker,
                       calcium_channel_blocker = tmp$Calcium_Channel_Blocker,
                       beta_agonist = tmp$Beta_Agonist,
                       digoxin = tmp$Digoxin)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/ECG/fact.cov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - no sex 
fact_cov2 <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       batch = tmp$Set,
                       smoking = tmp$ever_smoke,
                       beta_blocker = tmp$Beta_Blocker,
                       calcium_channel_blocker = tmp$Calcium_Channel_Blocker,
                       beta_agonist = tmp$Beta_Agonist,
                       digoxin = tmp$Digoxin)
write.table(fact_cov2, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/ECG/fact_nosex.cov", row.names=F, sep=' ', quote = F)


### ECG - without medication ###
## Create continuous covariate file
tmp <- ecg_no_med1
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
                       RR = tmp$RR_interval,
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
write.table(quant_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/ECG/quant.qcov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - without drugs included 
fact_cov <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       sex = tmp$sex,
                       batch = tmp$Set,
                       smoking = tmp$ever_smoke)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/ECG/fact.cov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - no sex 
fact_cov2 <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       batch = tmp$Set,
                       smoking = tmp$ever_smoke) 
write.table(fact_cov2, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/ECG/fact_nosex.cov", row.names=F, sep=' ', quote = F)




### RR - with medication ###
## Create continuous covariate file
tmp <- rr_med1
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
write.table(quant_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/RR/quant.qcov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - with drugs included 
fact_cov <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       sex = tmp$sex,
                       batch = tmp$Set,
                       smoking = tmp$ever_smoke)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/RR/fact.cov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - no sex 
fact_cov2 <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       batch = tmp$Set,
                       smoking = tmp$ever_smoke)
write.table(fact_cov2, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/Medication/RR/fact_nosex.cov", row.names=F, sep=' ', quote = F)


### RR - without medication ###
## Create continuous covariate file
tmp <- rr_no_med1
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
write.table(quant_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/RR/quant.qcov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - without drugs included 
fact_cov <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       sex = tmp$sex,
                       batch = tmp$Set,
                       smoking = tmp$ever_smoke)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/RR/fact.cov", row.names=F, sep=' ', quote = F)

## Create factor covariate file - no sex 
fact_cov2 <- data.frame(FID =tmp$FID,
                       IID = tmp$IID,
                       batch = tmp$Set,
                       smoking = tmp$ever_smoke)
write.table(fact_cov2, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Covariates/No_Medication/RR/fact_nosex.cov", row.names=F, sep=' ', quote = F)



########## PHENOTYPES ############

## Loop through ecg traits (excluding RR interval) to make phenotype files - Medication 
## extract variables 
var <- c("QRS", "QT", "PR", "JT")
## loop step 
for(i in var){ 
## Create phenotype data frame 
phen <- data.frame(FID =ecg_med1$FID,
                       IID = ecg_med1$ID,
                       phen = ecg_med1[,i])
write.table(phen, file=paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/Medication/ECG/", i, ".phen"), row.names=F, sep=' ', quote = F)
}


## Loop through ecg traits (excluding RR interval) to make phenotype files - No Medication 
## extract variables 
var <- c("QRS", "QT", "PR", "JT")
## loop step 
for(i in var){ 
## Create phenotype data frame 
phen <- data.frame(FID =ecg_no_med1$FID,
                       IID = ecg_no_med1$ID,
                       phen = ecg_no_med1[,i])
write.table(phen, file=paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/No_Medication/ECG/", i, ".phen"), row.names=F, sep=' ', quote = F)
}


## RR phenotype file - Medication 
## Create phenotype data frame 
phen <- data.frame(FID =rr_med1$FID,
                       IID = rr_med1$ID,
                       phen = rr_med1$RR)
write.table(phen, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/Medication/RR/RR.phen", row.names=F, sep=' ', quote = F)



## RR phenotype file - No Medication 
## Create phenotype data frame 
phen <- data.frame(FID =rr_no_med1$FID,
                       IID = rr_no_med1$ID,
                       phen = rr_no_med1$RR)
write.table(phen, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/No_Medication/RR/RR.phen", row.names=F, sep=' ', quote = F)



#### SEX-SPECIFIC ANALYSES ######

## Create .txt files for keep flags to keep males and females in sex-specific analyses 

## ECG 
# Merge with FIDs and tabulate sex information
ecg_keep <- merge(ecg,fids,by="ID")
table(ecg_keep$sex)
# Males
ecg_males <- ecg_keep[which(ecg_keep$sex %in% "M"),c("FID","ID")]
write.table(ecg_males, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Keep/ECG/males.txt", row.names = F, col.names = F, quote = F, sep = "\t")
# Females 
ecg_females <- ecg_keep[which(ecg_keep$sex %in% "F"),c("FID","ID")]
write.table(ecg_females, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Keep/ECG/females.txt", row.names = F, col.names = F, quote = F, sep = "\t")

## RR
# Merge with FIDs and tabulate sex information
rr_keep <- merge(rr,fids,by="ID")
table(rr_keep$sex)
# Males
rr_males <- rr_keep[which(rr_keep$sex %in% "M"),c("FID","ID")]
write.table(rr_males, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Keep/RR/males.txt", row.names = F, col.names = F, quote = F, sep = "\t")
# Females 
rr_females <- rr_keep[which(rr_keep$sex %in% "F"),c("FID","ID")]
write.table(rr_females, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Keep/RR/females.txt", row.names = F, col.names = F, quote = F, sep = "\t")
