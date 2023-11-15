## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Cox EWAS regression for Atrial Fibrillation 

## Load requisite libraries 
library(data.table)
library(lumi)
library(survival)

## set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/")

## choose chromosomes
chrx=21
chry=22

## Read in covariate files 
af=read.csv("Covariates/incident_af.csv")
names(af)[1] <- "ID"
af$Set=gsub("_.*","",af$Batch)

## Subset to the medication only data 
sub=read.table("osca_files/Covariates/No_Medication/AF/quant.qcov",header=T)
a1=af[which(af$ID %in% sub$IID),]


## Read in methylation data - first chromosome 
tmp=readRDS(paste0("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/GS20k_chr",chrx,"_mvals.rds"))

## Subset to those in atrial fibrillation file 
tmp1=tmp[,which(colnames(tmp) %in% a1$SampleID)]

## Change to beta values 
tmp1=m2beta(tmp1) 

## Read in methylation data - second chromosome 
tmp2=readRDS(paste0("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/GS20k_chr",chry,"_mvals.rds"))

## Subset to those in atrial fibrillation file 
tmp3=tmp2[,which(colnames(tmp2) %in% a1$SampleID)]

## Change to beta values 
tmp3=m2beta(tmp3) 

## Merge both matrices 
x=rbind(tmp1,tmp3)


## Match IDs between covariate and methylation dataframes 
ids = a1$SampleID
x1=x[,match(ids,colnames(x))]
## transpose methylation dataframe
x1=t(x1)

## Set up output for EWAS 
out=as.data.frame(matrix(nrow = ncol(x1),ncol=5))
names(out) <- c("Probe", "Beta", "SE", "P", "N")


## Set up loop for cox models 
for(i in 1:ncol(x1)){ 
## run cox model 
mod1=summary(coxph(Surv(a1$tte, a1$Case) ~ x1[,i] + a1$age + factor(a1$HTN) + factor(a1$sex) + factor(a1$diabetes_Y) + a1$avg_sys + a1$Mono + a1$NK + a1$Bcell + a1$CD4T + a1$CD8T + a1$PC1 + a1$PC2 + a1$PC3 + a1$PC4 + a1$PC5 + a1$PC6 + a1$PC7 + a1$PC8 + a1$PC9 + a1$PC10 + factor(a1$ever_smoke) + a1$height + a1$weight + factor(a1$Set)))
# store output
out[i,1] <- as.character(colnames(x1)[i])
out[i,2] <- mod1$coefficients[1,1]
out[i,3] <- mod1$coefficients[1,3]
out[i,4] <- mod1$coefficients[1,5]
out[i,5] <- mod1$n
print(i)
} 
fwrite(out, paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Cox_Outputs/No_Medication/ALL/chr_", chrx, "_", chry, "_probes.txt"),row.names = F)




## Subset to males 
a2 = a1[a1$sex %in% "M",]
x2=x1[which(row.names(x1) %in% a2$SampleID),]

out1=as.data.frame(matrix(nrow = ncol(x1),ncol=5))
names(out1) <- c("Probe", "Beta", "SE", "P", "N")

## Set up loop for cox models 
for(i in 1:ncol(x2)){ 
## run cox model 
mod1=summary(coxph(Surv(a2$tte, a2$Case) ~ x2[,i] + a2$age + factor(a2$HTN) + factor(a2$diabetes_Y) + a2$avg_sys + a2$Mono + a2$NK + a2$Bcell + a2$CD4T + a2$CD8T + a2$PC1 + a2$PC2 + a2$PC3 + a2$PC4 + a2$PC5 + a2$PC6 + a2$PC7 + a2$PC8 + a2$PC9 + a2$PC10 + factor(a2$ever_smoke) + a2$height + a2$weight + factor(a2$Set)))
# store output
out1[i,1] <- as.character(colnames(x2)[i])
out1[i,2] <- mod1$coefficients[1,1]
out1[i,3] <- mod1$coefficients[1,3]
out1[i,4] <- mod1$coefficients[1,5]
out1[i,5] <- mod1$n
print(i)
} 
fwrite(out1, paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Cox_Outputs/No_Medication/MALES/chr_", chrx, "_", chry, "_probes.txt"),row.names = F)



## Subset to females 
a2 = a1[a1$sex %in% "F",]
x2=x1[which(row.names(x1) %in% a2$SampleID),]

out2=as.data.frame(matrix(nrow = ncol(x1),ncol=5))
names(out2) <- c("Probe", "Beta", "SE", "P", "N")


## Set up loop for cox models 
for(i in 1:ncol(x2)){ 
## run cox model 
mod1=summary(coxph(Surv(a2$tte, a2$Case) ~ x2[,i] + a2$age + factor(a2$HTN) + factor(a2$diabetes_Y) + a2$avg_sys + a2$Mono + a2$NK + a2$Bcell + a2$CD4T + a2$CD8T + a2$PC1 + a2$PC2 + a2$PC3 + a2$PC4 + a2$PC5 + a2$PC6 + a2$PC7 + a2$PC8 + a2$PC9 + a2$PC10 + factor(a2$ever_smoke) + a2$height + a2$weight + factor(a2$Set)))
# store output
out2[i,1] <- as.character(colnames(x2)[i])
out2[i,2] <- mod1$coefficients[1,1]
out2[i,3] <- mod1$coefficients[1,3]
out2[i,4] <- mod1$coefficients[1,5]
out2[i,5] <- mod1$n
print(i)
} 
fwrite(out2, paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Cox_Outputs/No_Medication/FEMALES/chr_", chrx, "_", chry, "_probes.txt"),row.names = F)
