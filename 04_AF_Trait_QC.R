## CHARGE Cardiovascular consortium 
## Generation Scotland 
## Quality control of Incident Atrial Fibrillation 

##################################
### STEP 1. PREP INCIDENCE DATA ##
##################################

## Read in incidence data 
inc = read.csv("/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/2022-05-31_rob_ecg_exclusions.csv")
af <- inc[,c("id", "af", "ym")]
length(which(!is.na(af$af))) ##958 

## Separate prevalent from incident disease 

# Prevalent 
prev_af=af[which(af$af <= af$ym),] ##218 

## Incident 
inc_af=af[which(af$af > af$ym),] ##740 
inc_af=inc_af[order(inc_af$af),]



#######################################
## STEP 2. IDENTIFY CASES + CONTROLS ##
#######################################

## Read in everyone with DNAm data 
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k//GS20k_Targets.rds")

## Exclude those with prevalent Atrial Fibrillation 
samps1=samps[-which(samps$Sample_Name %in% prev_af$id),]

## Identify cases 
samps1$Case = NA 
samps1[which(samps1$Sample_Name %in% inc_af$id),"Case"] <- 1 

## Identify controls 
samps1[-which(samps1$Sample_Name %in% inc_af$id),"Case"] <- 0



#####################################
## STEP 3. CALCULATE TIME-TO-EVENT ##
#####################################

## Read in censoring info 
cens=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/censoring_oct2020.csv")
cens1 <- cens[,c("Sample_Name", "dead", "aged", "Age", "yob", "mob")]

## separate into alive and dead participants 
alive=cens1[cens1$dead==0,]
dead=cens1[cens1$dead==1,]

## Scenario 1) Case and Alive - diff between age at baseline and and age at event 
cases=merge(alive,samps1,by="Sample_Name")
cases=cases[cases$Case==1,]
cases=merge(cases,inc_af, by.x="Sample_Name",by.y="id", all.x=T)

## extract month and year of event 
cases$moe <- substring(cases$af, 5,6)
cases$yoe <- substring(cases$af, 1,4)
## substract from month and year of birth to get age at event 
cases$month_event1 = (as.numeric(cases$moe) - as.numeric(cases$mob))/12
cases$age_event1 = as.numeric(cases$yoe) - as.numeric(cases$yob)
cases$age_event = cases$age_event1 + cases$month_event1
## tidy up for harmonising with other groups 
cases$tte <- cases$age_event-cases$Age 
cases1 <- cases[,c("Sample_Name","Case","tte")]


## Scenario 2) Control and Alive - diff between age at baseline and age at censoring
controls=merge(alive,samps1,by="Sample_Name")
controls=controls[controls$Case==0,]
controls$tte <- controls$aged-controls$Age
controls1 <- controls[,c("Sample_Name","Case","tte")]


## Scenario 3) Participant has died - diff between age at baseline and age at death 
dead1=merge(dead,samps1,by="Sample_Name")
dead1$tte <- dead1$aged-dead1$Age
dead1 <- dead1[,c("Sample_Name","Case","tte")]


#####################################
## STEP 4. SAVE OUT PHENO FILE ######
#####################################

## Combine all three groups with different time-to-event criteria 
total=rbind(cases1,controls1)
total=rbind(total,dead1)
## save out file 
write.csv(total, "/Cluster_Filespace/Marioni_Group/Rob/CHARGE_ECG/Phenotypes/incident_af.csv",row.names = F)