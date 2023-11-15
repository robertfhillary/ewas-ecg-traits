## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Covariates for Incident Atrial Fibrillation  

##################################
### STEP 1. READ IN COVARIATES ###
##################################

## we need: age, sex, set, batch, height, weight, smoking, cell proportions, genetic principal components, drugs  

## age, sex, set, batch  
demo=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
names(demo)[2] <- "SampleID"

## height and weight 
bod=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/body.csv")
bod <- bod[,c("id", "height", "weight")]


## smoking (ever smoke)
smk=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/updated_smoking_jan_2019/ever_smoke.csv")
# recode to never, former, current 
# current 
smk[which(smk$ever_smoke == 1),"ever_smoke"] <- "Current"
# former 
smk[which(smk$ever_smoke == 2),"ever_smoke"] <- "Former"
smk[which(smk$ever_smoke == 3),"ever_smoke"] <- "Former"
# never
smk[which(smk$ever_smoke == 4),"ever_smoke"] <- "Never"
## na 
smk[which(smk$ever_smoke == 5),"ever_smoke"] <- NA


## cell proportions 
cell=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS_20k_DNAmAge.csv")
cell=cell[,c("SampleID", "NK", "Mono", "CD4T", "CD8T", "Bcell")]


## genetic PCs 
pc=read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/GS20K_ALL_MAF5_PCA.eigenvec")
names(pc)[1:2] <- c("FID", "IID")
names(pc)[3:ncol(pc)] <- paste0("PC", 1:20)
pc$FID <- NULL 


## SBP 
sbp=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/BPHR.csv") 
sbp=sbp[,c("id","avg_sys")]


## diabetes 
dis=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/disease.csv")
dis <- dis[,c("ID", "diabetes_Y")]


##################################
### STEP 2. MERGE COVARIATES #####
##################################

## merge covariates up to this point (before drug entries which are much more complicated and less complete)
cov=merge(demo, bod, by.x="Sample_Name",by.y="id",all.x=T)
cov=merge(cov, smk, by = "Sample_Name", all.x= T)
cov=merge(cov, cell, by = "SampleID", all.x = T)
cov=merge(cov, pc, by.x="Sample_Name", by.y="IID", all.x = T)
cov=merge(cov, sbp, by.x="Sample_Name", by.y="id", all.x= T)
cov=merge(cov, dis, by.x="Sample_Name",by.y="ID", all.x= T)



##################################
### STEP 3. EXTRACT DRUG INFO  ###
##################################

med = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/medicationv2.csv")
med2 = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/medicationv5.csv")


## clear up med file 
med$pills_na <- NULL
med$ointment_na <- NULL
med$inhaler_na <- NULL 
med$injection_na <- NULL 


## Extract IDs with digoxin intake at baseline 
output <- list()
for(i in 2:ncol(med)){
output[[i]] <- med[grep("DIGOXIN",med[,i]),"ID"]
}

tmp=med[which(med$ID %in% unique(unlist(output))),]
drug.ids<-tmp$ID


## Extract IDs with anti-hypertensive drug intake at baseline 
bp <- med2[,c("ID","bp")]
drug.ids2 <- bp[bp$bp%in%1,"ID"]

##################################
### STEP 4. ADD IN DRUG INFO  ####
##################################

## Create columns for each drug class in covariate file 
cov$Digoxin = NA 
cov$HTN = NA 

## Digoxin 
cov[which(cov$Sample_Name %in% med$ID),"Digoxin"] <- 0
cov[which(cov$Sample_Name %in% drug.ids),"Digoxin"] <- 1

## HTN 
cov[which(cov$Sample_Name %in% bp$ID),"HTN"] <- 0
cov[which(cov$Sample_Name %in% drug.ids2),"HTN"] <- 1


##########################################
# STEP 5. SAVE OUT FINAL COV FILE FOR AF #
##########################################

af=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Phenotypes/incident_af.csv")
af_cov=merge(af, cov, by="Sample_Name",all.y=T)
af_cov=af_cov[-which(is.na(af_cov$tte)),] ## 18201
write.csv(af_cov, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Covariates/incident_af.csv", row.names = F) ## 18201


