## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Preparation of Final Files to send 


##########################################
### STEP 1. PREPARE FINAL FILES  #########
##########################################

## load requisite libraries
library(data.table)

## set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Outputs/")

## get list of 'good' probes
lis=readRDS("../../common_cpgs.rds")

#### ECG + RR ####

### MEDICATION DATA #####

## extract files 
files=list.files("Medication/", ".linear")

## Loop through files to format for file sharing 
for(i in files){ 
## read in file
tmp = as.data.frame(fread(paste0("Medication/", i)))
## extract required columns
tmp1 <- tmp[,c("Probe","b", "se", "p", "NMISS")]
## rename columns to required format 
names(tmp1) <- c("Probe","Beta", "SE", "P", "N")
## reorder by top CpGs
tmp2 <- tmp1[order(tmp1$P),]
## extract file name
nm = gsub(".linear", "", i)
# remove bad cpgs 
tmp2 <- tmp2[which(tmp2$Probe %in% lis$Probe),]
## save out file 
fwrite(tmp2, paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Final_Files/Medication/", nm, ".csv"), row.names = F) 
## print i to denote completion
print(i)
}




### NO MEDICATION DATA #####

## extract files 
files=list.files("No_Medication/", ".linear")

## Loop through files to format for file sharing 
for(i in files){ 
## read in file
tmp = as.data.frame(fread(paste0("No_Medication/", i)))
## extract required columns
tmp1 <- tmp[,c("Probe","b", "se", "p", "NMISS")]
## rename columns to required format 
names(tmp1) <- c("Probe","Beta", "SE", "P", "N")
## reorder by top CpGs
tmp2 <- tmp1[order(tmp1$P),]
## extract file name
nm = gsub(".linear", "", i)
# remove bad cpgs 
tmp2 <- tmp2[which(tmp2$Probe %in% lis$Probe),]
## save out file 
fwrite(tmp2, paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Final_Files/No_Medication/", nm, ".csv"), row.names = F) 
## print i to denote completion
print(i)
}


## AF 

setwd("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Cox_Outputs/Medication")

### MEDICATION DATA #####

## extract files 
files1=list.files(".", ".")

for(i in 1:length(files1)){
## Loop through files to format for file sharing 
files=list.files(paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Cox_Outputs/Medication/", files1[i]), ".")
## change working directory to read in relevant files
## combine files and reorder by p value 
tables <- lapply(paste0(files1[i], "/", files), read.csv, header = TRUE)
combined.df <- do.call(rbind , tables)
combined.df <- combined.df[order(combined.df$P),]
##save out file
fwrite(combined.df, paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Final_Files/Medication/AF_", files1[i], "_GS_EA_RFH_090621.csv"), row.names = F) 
} 


setwd("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Cox_Outputs/No_Medication")

### NO MEDICATION DATA #####

## extract files 
files1=list.files(".", ".")

for(i in 1:length(files1)){
## Loop through files to format for file sharing 
files=list.files(paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Cox_Outputs/No_Medication/", files1[i]), ".")
## change working directory to read in relevant files
## combine files and reorder by p value 
tables <- lapply(paste0(files1[i], "/", files), read.csv, header = TRUE)
combined.df <- do.call(rbind , tables)
combined.df <- combined.df[order(combined.df$P),]
##save out file
fwrite(combined.df, paste0("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Final_Files/No_Medication/AF_", files1[i], "_GS_EA_RFH_090621.csv"), row.names = F) 
} 




##########################################
### STEP 2. GZIP FINAL FILES  ############
##########################################

gzip /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Final_Files/Medication/*
gzip /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Final_Files/No_Medication/*





