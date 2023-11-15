## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Running Linear Regression OSCA models 

## set screen 
#screen -S osca 

########### ECG TRAITS ############################

######### Medication ################ 

### ALL ####

cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/Medication/ECG/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--pheno $i \
--qcovar ../../../Covariates/Medication/ECG/quant.qcov \
--covar ../../../Covariates/Medication/ECG/fact.cov \
--fast-linear \
--out ../../../Outputs/Medication/${A}_ALL_GS_EA_RFH_090621 \
--methylation-beta
done 



### MALES ####
cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/Medication/ECG/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--keep ../../../Keep/ECG/males.txt \
--pheno $i \
--qcovar ../../../Covariates/Medication/ECG/quant.qcov \
--covar ../../../Covariates/Medication/ECG/fact_nosex.cov \
--fast-linear \
--out ../../../Outputs/Medication/${A}_MALE_GS_EA_RFH_090621 \
--methylation-beta
done 


### FEMALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/Medication/ECG/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--keep ../../../Keep/ECG/females.txt \
--pheno $i \
--qcovar ../../../Covariates/Medication/ECG/quant.qcov \
--covar ../../../Covariates/Medication/ECG/fact_nosex.cov \
--fast-linear \
--out ../../../Outputs/Medication/${A}_FEMALE_GS_EA_RFH_090621 \
--methylation-beta
done 



######### No Medication ################ 

### ALL ####

cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/No_Medication/ECG/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--pheno $i \
--qcovar ../../../Covariates/No_Medication/ECG/quant.qcov \
--covar ../../../Covariates/No_Medication/ECG/fact.cov \
--fast-linear \
--out ../../../Outputs/No_Medication/${A}_ALL_GS_EA_RFH_090621 \
--methylation-beta
done 


### MALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/No_Medication/ECG/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--keep ../../../Keep/ECG/males.txt \
--pheno $i \
--qcovar ../../../Covariates/No_Medication/ECG/quant.qcov \
--covar ../../../Covariates/No_Medication/ECG/fact_nosex.cov \
--fast-linear \
--out ../../../Outputs/No_Medication/${A}_MALE_GS_EA_RFH_090621 \
--methylation-beta
done 

### FEMALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/No_Medication/ECG/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--keep ../../../Keep/ECG/females.txt \
--pheno $i \
--qcovar ../../../Covariates/No_Medication/ECG/quant.qcov \
--covar ../../../Covariates/No_Medication/ECG/fact_nosex.cov \
--fast-linear \
--out ../../../Outputs/No_Medication/${A}_FEMALE_GS_EA_RFH_090621 \
--methylation-beta
done 




########### RR INTERVAL ############################

######### Medication ################ 

### ALL ####


cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/Medication/RR/
  
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--pheno RR.phen \
--qcovar ../../../Covariates/Medication/RR/quant.qcov \
--covar ../../../Covariates/Medication/RR/fact.cov \
--fast-linear \
--out ../../../Outputs/Medication/RR_ALL_GS_EA_RFH_090621 \
--methylation-beta


### MALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/Medication/RR/
  
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--pheno RR.phen \
--keep ../../../Keep/RR/males.txt \
--qcovar ../../../Covariates/Medication/RR/quant.qcov \
--covar ../../../Covariates/Medication/RR/fact_nosex.cov \
--fast-linear \
--out ../../../Outputs/Medication/RR_MALE_GS_EA_RFH_090621 \
--methylation-beta


### FEMALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/Medication/RR/
  
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--pheno RR.phen \
--keep ../../../Keep/RR/females.txt \
--qcovar ../../../Covariates/Medication/RR/quant.qcov \
--covar ../../../Covariates/Medication/RR/fact_nosex.cov \
--fast-linear \
--out ../../../Outputs/Medication/RR_FEMALE_GS_EA_RFH_090621 \
--methylation-beta



######### No Medication ################ 


cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/No_Medication/RR/
  
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--pheno RR.phen \
--qcovar ../../../Covariates/No_Medication/RR/quant.qcov \
--covar ../../../Covariates/No_Medication/RR/fact.cov \
--fast-linear \
--out ../../../Outputs/No_Medication/RR_ALL_GS_EA_RFH_090621 \
--methylation-beta

### MALES ####


cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/No_Medication/RR/
  
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--pheno RR.phen \
--keep ../../../Keep/RR/males.txt \
--qcovar ../../../Covariates/No_Medication/RR/quant.qcov \
--covar ../../../Covariates/No_Medication/RR/fact_nosex.cov \
--fast-linear \
--out ../../../Outputs/No_Medication/RR_MALE_GS_EA_RFH_090621 \
--methylation-beta



### FEMALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/osca_files/Phenotypes/No_Medication/RR/
  
osca_Linux \
--linear \
--befile ../../../bvals-norm20k-18413-831733 \
--pheno RR.phen \
--keep ../../../Keep/RR/females.txt \
--qcovar ../../../Covariates/No_Medication/RR/quant.qcov \
--covar ../../../Covariates/No_Medication/RR/fact_nosex.cov \
--fast-linear \
--out ../../../Outputs/No_Medication/RR_FEMALE_GS_EA_RFH_090621 \
--methylation-beta
