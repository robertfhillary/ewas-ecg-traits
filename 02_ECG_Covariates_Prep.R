## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Covariates for 5 ECG traits 

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


##################################
### STEP 2. MERGE COVARIATES #####
##################################

## merge covariates up to this point (before drug entries which are much more complicated and less complete)
cov=merge(demo, bod, by.x="Sample_Name",by.y="id",all.x=T)
cov=merge(cov, smk, by = "Sample_Name", all.x= T)
cov=merge(cov, cell, by = "SampleID", all.x = T)
cov=merge(cov, pc, by.x="Sample_Name", by.y="IID", all.x = T)


##################################
### STEP 3. EXTRACT DRUG INFO  ###
##################################

med = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/medicationv2.csv")
med2 = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/medicationv5.csv")
cov_drugs = read.csv("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/ECG_drugs_covariates.csv")


## clear up med file 
med$pills_na <- NULL
med$ointment_na <- NULL
med$inhaler_na <- NULL 
med$injection_na <- NULL 


## Beta Blocker 
excl_tmp <- cov_drugs[which(cov_drugs$ATC_Group %in% "C07"),]

## set up loop to query 
output <- list()

## query each C07 agent 
for(j in excl_tmp$Agent_allcaps){
drug <- j 
print(drug)
## query each column with drug information (pills)
for(i in 2:ncol(med)){
output[[j]][[i]] <- med[grep(drug,med[,i]),"ID"]
}
}


## Obtain IDs for covariates  
var=c("PROPRANOLOL","TIMOLOL","SOTALOL","NADOLOL","METOPROLOL","ATENOLOL","BISOPROLOL","CELIPROLOL","NEBIVOLOL","LABETALOL","CARVEDILOL")
drug.ids <- list()
for(j in var){ 
tmp=med[which(med$ID %in% unique(unlist(output[[j]]))),]
drug.ids[[j]]<-tmp$ID
} 
drug.ids<-unique(unlist(drug.ids))



## Calcium Channel Blocker 
## subset to C08D drugs 
excl_tmp <- cov_drugs[which(cov_drugs$ATC_Group %in% "C08D"),]

## set up loop to query 
output <- list()

## query each C08D agent 
for(j in excl_tmp$Agent_allcaps){
drug <- j 
print(drug)
## query each column with drug information (pills)
for(i in 2:ncol(med)){
output[[j]][[i]] <- med[grep(drug,med[,i]),"ID"]
}
}


## Obtain IDs for covariates 
var=c("DILTIAZEM","VERAPAMIL")
drug.ids2 <- list()
for(j in var){ 
tmp=med[which(med$ID %in% unique(unlist(output[[j]]))),]
drug.ids2[[j]]<-tmp$ID
} 
drug.ids2<-unique(unlist(drug.ids2))

## Combine with existing drug list
total=unique(c(drug.ids, drug.ids2))


## Beta Agonist 
excl_tmp <- cov_drugs[which(cov_drugs$ATC_Group %in% "R03A" | cov_drugs$ATC_Group %in% "R03C"),]

## set up loop to query 
output <- list()

## query each R03A/R03C agent 
for(j in excl_tmp$Agent_allcaps){
drug <- j 
print(drug)
## query each column with drug information (pills)
for(i in 2:ncol(med)){
output[[j]][[i]] <- med[grep(drug,med[,i]),"ID"]
}
}

## Obtain IDs to exclude 
var=c("EPHEDRINE","FORMOTEROL","SALMETEROL","TERBUTALINE","SALBUTAMOL","EPINEPHRINE")
drug.ids3 <- list()
for(j in var){ 
tmp=med[which(med$ID %in% unique(unlist(output[[j]]))),]
drug.ids3[[j]]<-tmp$ID
} 
drug.ids3<-unique(unlist(drug.ids3))

## Combine with existing drug list
total=unique(c(total, drug.ids3))


## Digoxin 
output <- list()
for(i in 2:ncol(med)){
output[[i]] <- med[grep("DIGOXIN",med[,i]),"ID"]
}

tmp=med[which(med$ID %in% unique(unlist(output))),]
drug.ids4<-tmp$ID

## Combine with existing drug list
total=unique(c(total, drug.ids4)) ##724 


##################################
### STEP 4. ADD IN DRUG INFO  ####
##################################

## Create columns for each drug class in covariate file 
cov$Beta_Blocker = NA 
cov$Calcium_Channel_Blocker = NA 
cov$Beta_Agonist = NA 
cov$Digoxin = NA 

## Beta Blocker 
cov[which(cov$Sample_Name %in% med$ID),"Beta_Blocker"] <- 0
cov[which(cov$Sample_Name %in% drug.ids),"Beta_Blocker"] <- 1

## Calcium Channel Blocker 
cov[which(cov$Sample_Name %in% med$ID),"Calcium_Channel_Blocker"] <- 0
cov[which(cov$Sample_Name %in% drug.ids2),"Calcium_Channel_Blocker"] <- 1

## Beta Agonist 
cov[which(cov$Sample_Name %in% med$ID),"Beta_Agonist"] <- 0
cov[which(cov$Sample_Name %in% drug.ids3),"Beta_Agonist"] <- 1

## Digoxin 
cov[which(cov$Sample_Name %in% med$ID),"Digoxin"] <- 0
cov[which(cov$Sample_Name %in% drug.ids4),"Digoxin"] <- 1


##########################################
# STEP 5. SAVE OUT FINAL COV FILE FOR RR #
##########################################

rr=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Phenotypes/rr_interval.csv")
rr <- rr[,c("ID","RR")]
rr_cov=merge(rr, cov, by.x="ID",by.y="Sample_Name",all.y=T)
rr_cov=rr_cov[-which(is.na(rr_cov$RR)),] ## 17259
write.csv(rr_cov, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Covariates/rr_interval.csv", row.names = F) ## 17259


##########################################
# STEP 5. ADD IN RR FOR OTHER ECG TRAITS #
##########################################

ecg_phen = read.csv("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Phenotypes/ecg_traits_no_RR.csv") 
ecg_phen <- ecg_phen[,c("ID","QRS_duration", "QT_interval", "PR_interval", "JT")]

## add in RR interval
ecg = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/ECG.csv")
ecg <- ecg[,c("ID","Heart_Rate")]
ecg$RR_interval <- 60000/ecg$Heart_Rate
ecg_phen = merge(ecg_phen, ecg[,c("ID", "RR_interval")], by="ID", all.x= T)

## combine with other covariates and save out 
ecg_cov=merge(ecg_phen, cov, by.x="ID",by.y="Sample_Name",all.y=T)
ecg_cov=ecg_cov[-which(is.na(ecg_cov$QRS_duration)),] 
write.csv(ecg_cov, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Covariates/ecg_traits.csv", row.names = F)  ## 17486
