## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Quality control of 5 ECG traits 

##################################
### STEP 1. READ IN ECG TRAITS ###
##################################

## read in traits 
ecg = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/ECG.csv")

## extract four ECG traits of interest 
ecg1 = ecg[,c("ID", "age", "sex", "QRS_duration", "QT_interval", "PR_interval", "Heart_Rate", "Minn_Group8_1", "Minn_Group8_2") ]

## calculate fifth ECG trait of interest - JT Interval (QT-QRS)
ecg1$JT <- as.numeric(ecg1$QT_interval - ecg1$QRS_duration)



####################################################
### STEP 2. QC ECG TRAITS - REMOVE OUTLIERS ########
####################################################

## CHARGE QC steps 

## 1. QRS >= 120 
length(ecg1[which(ecg1$QRS_duration >= 120),"QRS_duration"]) ## 389
ecg1[which(ecg1$QRS_duration >= 120),"QRS_duration"] <- NA 


## 2. Heart_Rate < 40 or > 120 
length(ecg1[which(ecg1$Heart_Rate < 40),"Heart_Rate"]) ## 57
length(ecg1[which(ecg1$Heart_Rate > 120),"Heart_Rate"]) ## 3 ## 60 total

ecg1[which(ecg1$Heart_Rate < 40),"Heart_Rate"] <- NA 
ecg1[which(ecg1$Heart_Rate > 120),"Heart_Rate"] <- NA 

## 3. PR_interval >320  
length(ecg1[which(ecg1$PR_interval > 320),"PR_interval"]) ## 9 
ecg1[which(ecg1$PR_interval > 320),"PR_interval"] <- NA 


## Additional QC steps - some negative values detected for ECG traits 

## QRS 
length(ecg1[which(ecg1$QRS_duration < 0),"QRS_duration"]) ## 3 ## 392 total 
ecg1[which(ecg1$QRS_duration < 0),"QRS_duration"] <- NA 

## QT 
length(ecg1[which(ecg1$QT_interval < 0),"QT_interval"]) ## 3 ## 3 total 
ecg1[which(ecg1$QT_interval < 0),"QT_interval"] <- NA 

## PR 
length(ecg1[which(ecg1$PR_interval < 0),"PR_interval"]) ## 54 ## 63 total 
ecg1[which(ecg1$PR_interval < 0),"PR_interval"] <- NA 

## JT 
length(ecg1[which(ecg1$JT < 0),"JT"]) #- none we skip 

## Heart_Rate 
length(ecg1[which(ecg1$Heart_Rate < 0),"Heart_Rate"]) ## 0 ## 60 total 




####################################################
### STEP 3. EXCLUSION CRITERIA #####################
####################################################

excl = read.csv("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/2022-05-31_rob_ecg_exclusions.csv")

### 3.1 Prevalent Heart Failure, Myocardial Infarction, Heart Block, Pacemaker, WFW syndrome  

## heart failure 
hf=excl[,c("id","hf","ym")]
hf=hf[complete.cases(hf$hf),]
hf_ids=hf[which(hf$hf <= hf$ym),]
nrow(hf_ids) #101

nrow(ecg1) ##21283
ecg1=ecg1[-which(ecg1$ID %in% hf_ids$id),] ##21198 

## Myocardial Infarction 
mi=excl[,c("id","mi","ym")]
mi=mi[complete.cases(mi$mi),]
mi_ids=mi[which(mi$mi <= mi$ym),]
nrow(mi_ids) #261

ecg1=ecg1[-which(ecg1$ID %in% mi_ids$id),] ##21014


## 2nd degree heart block 
hb2d=excl[,c("id","hb2d","ym")]
hb2d=hb2d[complete.cases(hb2d$hb2d),]
hb2d_ids=hb2d[which(hb2d$hb2d <= hb2d$ym),]
nrow(hb2d_ids) #4

ecg1=ecg1[-which(ecg1$ID %in% hb2d_ids$id),] ##21011


## 3rd degree heart block 
hb3d=excl[,c("id","hb3d","ym")]
hb3d=hb3d[complete.cases(hb3d$hb3d),]
hb3d_ids=hb3d[which(hb3d$hb3d <= hb3d$ym),]
nrow(hb3d_ids) #6

ecg1=ecg1[-which(ecg1$ID %in% hb3d_ids$id),] ##21009


## Pacemaker 
pc=excl[,c("id","pacemaker","ym")]
pc=pc[complete.cases(pc$pacemaker),]
pc_ids=pc[which(pc$pacemaker <= pc$ym),]
nrow(pc_ids) #41

ecg1=ecg1[-which(ecg1$ID %in% pc_ids$id),] ##20990

## WFW 
wpw=excl[,c("id","wpw","ym")]
wpw=wpw[complete.cases(wpw$wpw),]
wpw_ids=wpw[which(wpw$wpw <= wpw$ym),]
nrow(wpw_ids) #6

ecg1=ecg1[-which(ecg1$ID %in% wpw_ids$id),] ##20984


### 3.2 Prevalent Atrial Fibrillation based on ECG data - Minnesota code 8-3-1
ecg1=ecg1[-which(ecg1$Minn_Group8_1 %in% "831"),] ##20867


## 3.3 Medication data 
med = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/medicationv2.csv")
med2 = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/medicationv5.csv")
excl_drugs = read.csv("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/ECG_drugs_to_exclude.csv")

## clear up med file 
med$pills_na <- NULL
med$ointment_na <- NULL
med$inhaler_na <- NULL 
med$injection_na <- NULL 

## C01B
## set up loop to query 
output <- list()

## query each C01B agent 
for(j in excl_drugs$Agent_allcaps){
drug <- j 
print(drug)
## query each column with drug information (pills)
for(i in 2:29){
output[[j]][[i]] <- med[grep(drug,med[,i]),"ID"]
}
}

## check output (visual inspection)
tmp1=lapply(output, function(x) do.call("cbind", x))

## FLECAINIDE, AMIODARONE - only ones present 

## FLECAINIDE: 143289, 29520, 98814
## AMIODARONE: 151605, 82977, 177008, 134500, 171862


## Diagnostics 
## AMIODARONE 
check1=med[med$ID %in% c(151605, 82977, 177008,134500,171862),]
## FLECAINIDE 
check2=med[med$ID %in% c(143289, 29520, 98814),]

## Exclude all - pass QC checks 
ecg1[-which(ecg1$ID %in% c(143289, 29520, 98814.151605, 82977, 177008,134500,171862)),] ##20865 



####################################################
## STEP 4. EXTRACT ECG TRAITS (EXCEPT HEART RATE) ##
####################################################

## All ECG traits except heart rate should have same analysis n 

ecg_final <- ecg1[,c("ID", "age", "sex", "QRS_duration", "QT_interval", "PR_interval","JT")] ##20865
ecg_final.2 <- ecg_final[complete.cases(ecg_final),] ##20489 
write.csv(ecg_final.2, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Phenotypes/ecg_traits_no_RR.csv",row.names=F)



####################################################
## STEP 5. ADDITIONAL QC FOR HEART RATE (RR INT) ###
####################################################

ecg_final <- ecg1[,c("ID", "age", "sex", "Heart_Rate")] ##20865

## Additional exclusion criteria based on other drugs - just for Heart Rate/RR Interval 
med = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/medicationv2.csv")
med2 = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/medicationv5.csv")
excl_drugs = read.csv("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/ECG_drugs_covariates.csv")

## clear up med file 
med$pills_na <- NULL
med$ointment_na <- NULL
med$inhaler_na <- NULL 
med$injection_na <- NULL 


## 5.1 C07
## subset to C07 drugs 
excl_tmp <- excl_drugs[which(excl_drugs$ATC_Group %in% "C07"),]

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


## check output (visual inspection)
lapply(output, function(x) do.call("cbind", x))

## drugs to check 
var=c("PROPRANOLOL","TIMOLOL","SOTALOL","NADOLOL","METOPROLOL","ATENOLOL","BISOPROLOL","CELIPROLOL","NEBIVOLOL","LABETALOL","CARVEDILOL")
output.tmp <- list()

## loop through each drug 
for(j in var){ 
tmp=med[which(med$ID %in% unique(unlist(output[[j]]))),]
for(i in 2:ncol(tmp)){
## get name of column 
nm = names(tmp)[i]
##subset df 
tmp1=tmp[,c(1,i)]
## store unique drug names where drug of interest is mentioned 
output.tmp[[j]][[i]] <- unique(tmp1[grep(j, tmp1[,nm]),nm])
}
} 

## Obtain IDs to exclude 
var=c("PROPRANOLOL","TIMOLOL","SOTALOL","NADOLOL","METOPROLOL","ATENOLOL","BISOPROLOL","CELIPROLOL","NEBIVOLOL","LABETALOL","CARVEDILOL")
drug.ids <- list()
for(j in var){ 
tmp=med[which(med$ID %in% unique(unlist(output[[j]]))),]
drug.ids[[j]]<-tmp$ID
} 
drug.ids<-unique(unlist(drug.ids))



## 5.2 C08D
## subset to C08D drugs 
excl_tmp <- excl_drugs[which(excl_drugs$ATC_Group %in% "C08D"),]

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


## check output (visual inspection)
lapply(output, function(x) do.call("cbind", x))

## drugs to check 
var=c("DILTIAZEM","VERAPAMIL")
output.tmp <- list()

## loop through each drug 
for(j in var){ 
tmp=med[which(med$ID %in% unique(unlist(output[[j]]))),]
for(i in 2:ncol(tmp)){
## get name of column 
nm = names(tmp)[i]
##subset df 
tmp1=tmp[,c(1,i)]
## store unique drug names where drug of interest is mentioned 
output.tmp[[j]][[i]] <- unique(tmp1[grep(j, tmp1[,nm]),nm])
}
} 

## Obtain IDs to exclude 
var=c("DILTIAZEM","VERAPAMIL")
drug.ids2 <- list()
for(j in var){ 
tmp=med[which(med$ID %in% unique(unlist(output[[j]]))),]
drug.ids2[[j]]<-tmp$ID
} 
drug.ids2<-unique(unlist(drug.ids2))

## Combine with existing drug list
total=unique(c(drug.ids, drug.ids2))



## 5.3 R03A and R03C
## subset to R03A/R03C drugs 
excl_tmp <- excl_drugs[which(excl_drugs$ATC_Group %in% "R03A" | excl_drugs$ATC_Group %in% "R03C"),]

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


## check output (visual inspection)
lapply(output, function(x) do.call("cbind", x))

## drugs to check 
var=c("EPHEDRINE","FORMOTEROL","SALMETEROL","TERBUTALINE","SALBUTAMOL","EPINEPHRINE")
output.tmp <- list()

## loop through each drug 
for(j in var){ 
tmp=med[which(med$ID %in% unique(unlist(output[[j]]))),]
for(i in 2:ncol(tmp)){
## get name of column 
nm = names(tmp)[i]
##subset df 
tmp1=tmp[,c(1,i)]
## store unique drug names where drug of interest is mentioned 
output.tmp[[j]][[i]] <- unique(tmp1[grep(j, tmp1[,nm]),nm])
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


## 5.4 DIGOXIN 
## query digoxin 
output <- list()
for(i in 2:ncol(med)){
output[[i]] <- med[grep("DIGOXIN",med[,i]),"ID"]
}

## check output (visual inspection)
output.tmp<-list()
tmp=med[which(med$ID %in% unique(unlist(output))),]
for(i in 2:ncol(tmp)){
## get name of column 
nm = names(tmp)[i]
##subset df 
tmp1=tmp[,c(1,i)]
## store unique drug names where drug of interest is mentioned 
output.tmp[[i]] <- unique(tmp1[grep("DIGOXIN", tmp1[,nm]),nm])
}

## Obtain IDs to exclude 
tmp=med[which(med$ID %in% unique(unlist(output))),]
drug.ids4<-tmp$ID

## Combine with existing drug list
total=unique(c(total, drug.ids4)) ##724



####################################################
## STEP 6. EXTRACT HEART RATE + CONVERT TO RR ######
####################################################

## Exclude individuals with self-reported drug use from heart rate analysis n 
nrow(ecg_final) ##20865
ecg_final.hr = ecg_final[-which(ecg_final$ID %in% total),]
ecg_final.hr <- ecg_final.hr[complete.cases(ecg_final.hr),]
nrow(ecg_final.hr) ##20255 

## Convert to RR Interval and write out file
ecg_final.hr$RR <- 60000/ecg_final.hr$Heart_Rate
write.csv(ecg_final.hr, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Phenotypes/rr_interval.csv",row.names=F)

