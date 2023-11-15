#### CHARGE Cardiovascular consortium 
## Generation Scotland 
## Preparation of OSCA methylation files 


##########################################
### STEP 1. MAKE .TXT METHYLATION FILE ###
##########################################

## Load in requisite libraries 
library(lumi)
library(data.table)

## Read in Methylation File 
gs20k <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/normalised-betas.rds")

## Subset to those passing QC 
demo=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
gs20=gs20k[,which(colnames(gs20k) %in% demo$Sample_Sentrix_ID)]

## Match order of IDs 
ids=demo$Sample_Sentrix_ID
gs20=gs20[,match(ids,colnames(gs20))]
table(colnames(gs20)==demo$Sample_Sentrix_ID)

## Rename colnames of Methylation file according to Sample_Name 
colnames(gs20) <- demo$Sample_Name 

## Extract FID information 
ped <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/pedigree.csv")
ped <- ped[which(ped$volid %in% demo$Sample_Name),]
ped <- ped[,c("famid","volid")]
names(ped) <- c("FID", "IID")
ids=demo$Sample_Name
ped <- ped[match(ids,ped$IID),]
table(colnames(gs20)==ped$IID)

## Create dataframe for .txt file needed for OSCA 
osca_dat <- data.frame(FID=ped$FID, IID=ped$IID)
tgs20 <- t(gs20)
rm(gs20)
gc()
osca_dat1 <- cbind(osca_dat, tgs20)
rm(tsg20)
gc()

## Export as .txt for OSCA make-bod
fwrite(osca_dat1, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/bvals-norm20k-18413-831733.txt", quote = FALSE, row.names = FALSE, sep = " ")


#############################################
### STEP 2. MAKE BINARY METHYLATION FILE ####
#############################################

### In command line
osca_Linux --efile /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/bvals-norm20k-18413-831733.txt --methylation-beta --make-bod --out /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/bvals-norm20k-18413-831733


#############################################
### STEP 3. REMAKE .OPI ANNOTATION FILE #####
#############################################

### Alter .opi file
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
osca_dat <- as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/bvals-norm20k-18413-831733.opi", header = FALSE))
rownames(osca_dat) <- osca_dat$V2
opi <- anno[rownames(osca_dat),c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- gsub("X", "23", opi$chr)
opi$chr <- gsub("Y", "24", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name  <- as.factor(opi$UCSC_RefGene_Name )
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name ==""), "UCSC_RefGene_Name"] <- NA
write.table(opi, file="/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/bvals-norm20k-18413-831733.opi", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')


