# AnthropoAge, a novel approach to integrate body composition into the estimation of biological age 
# Data Analysis: Carlos A. Fermín-Martínez, Alejandro Márquez-Salinas & Omar Yaxmehen Bello-Chavolla
# Latest version of Analysis 16-Sept-2021
# Any question regarding analysis contact Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

####---- Database management ----####
pacman::p_load(haven, tidyverse, ggpubr, lmtest, nortest, gtools, data.table, caret, glmnet, survival, flextable, blandr, BlandAltmanLeh,
               rms, bestNormalize, flexsurv, pROC, timeROC, fmsb, factoextra, rgl, gridExtra,  nhanesA, wesanderson,forestmodel, ggedit,
               FactoMineR, fpc, NbClust, ggimage, glmnet, ggsci, survminer, cluster, ggplotify, UpSetR, nortest, viridis, officer, magrittr)

#Extra functions
conc95 <- function(x){
  y <- (summary(x)$concordance[2])*1.96
  c <- summary(x)$concordance[1]
  c_low <- c-y; c_up <- c+y
  `names<-`(c(c,c_low,c_up),c("Concordance","Lower 95-CI","Upper 95-CI"))}

#Para que se vean los acentos, guardar como WINDOWS-1252 (por favor)
setwd("C:/Users/investigacion/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/AnthropoAge")
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/AnthropoAge")
setwd("~/UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ALEJANDRO MARQUEZ SALINAS - Antropometría")
setwd("C:/Users/facmed/UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ALEJANDRO MARQUEZ SALINAS - Antropometría")

#Load NHANES III
NHANES<-read.csv("nhanes4.csv",na = c("", "N/A", "NA","na", "#N/A", "88888", "8888", "888888", "888","5555","999","9999", "99999", "999999"))
NHANES3<-read.csv("nhanes3.csv", na = c("", "N/A", "NA","na", "#N/A", "88888", "8888", "888888", "888","5555","999","9999", "99999", "999999"))
mortalidad3<-read.csv("nhanes3_mort.csv")
mortalidad<-read.csv("nhanes_mortalidad.csv")
source("functions_predict.R")
mortalidad$SEQN<-mortalidad$seqn
mortalidad3$SEQN<-mortalidad3$seqn

NHANES0<- NHANES3 %>% filter(HSAGEIR.x>=20)%>% 
  mutate(BMI=BMPBMI, Thigh_circumference=BMPTHICI, Arm_circumference=BMPARMC, Weight= BMPWT, Waist=BMPWAIST, Height=BMPHT,
         Triceps_skinfold = BMPTRI, Subscapular_skinfold = BMPSUB, Leg_length=BMPLEG, Arm_length=BMPARML, ICE=BMPWAIST/BMPHT, 
         Diabetes=HAD1, Asthma=HAC1E, Overweight=HAM15A, Arthritis=HAC1A, Heart_failure=HAC1D, 
         Heart_attack=HAF10, Stroke=HAC1D, Emphysema=HAC1D, Bronchitis=HAC1D, Malignancy=HAC1D)

NHANES0$Sex <- factor(NHANES0$HSSEX.x, levels= c(1,2),labels= c("Men", "Women"))
NHANES0$Ethnicity<-NHANES0$DMARETHN.x
NHANES0$Ethnicity[NHANES0$Ethnicity==4]<-5
NHANES0$Ethnicity[NHANES0$DMAETHNR==2]<-4
NHANES0$Age<-NHANES0$HSAGEIR.x

NHANES0$PhenoAge = 141.5 + ((log(-0.00553*log(1-(1-exp((-1.51714*exp(-19.907-0.0336*NHANES0$AMP+0.0095*NHANES0$CEPSI +0.1953*NHANES0$G1PSI+
                                                                       0.0954*log(NHANES0$CRP)-0.0120*NHANES0$LMPPCNT+0.0268*NHANES0$MVPSI+
                                                                       0.3306*NHANES0$RWP+0.00188*NHANES0$AP+0.0554*NHANES0$WCP+
                                                                       0.0804*NHANES0$HSAGEIR.x))/(0.0076927)))))))/(0.09165)
NHANES0 <- merge(NHANES0,mortalidad3,by="SEQN")

NHANES1 <- NHANES %>% 
  mutate(DXDTOFAT_N = (DXDTOFAT/1000)/((BMXWT/100)^2), DXXTRFAT_N = (DXXTRFAT/1000)/((BMXWT/100)^2), DXXHEFAT_N = (DXXHEFAT/1000)/((BMXWT/100)^2),
         DXXLAFAT_N = (DXXLAFAT/1000)/((BMXWT/100)^2), DXXRAFAT_N = (DXXRAFAT/1000)/((BMXWT/100)^2), DXXLLFAT_N = (DXXLLFAT/1000)/((BMXWT/100)^2),
         DXXRLFAT_N = (DXXRLFAT/1000)/((BMXWT/100)^2), DXXTRFAT_DXDTOFAT = (DXXTRFAT)/(DXDTOFAT), Weight = BMXWT, Waist = BMXWAIST, Height = BMXHT,
         BMI = BMXBMI, Calf_circumference = BMXCALF, Arm_circumference = BMXARMC, Thigh_circumference = BMXTHICR, Triceps_skinfold = BMXTRI, Subscapular_skinfold = BMXSUB,
         Leg_length=BMXLEG, ICE = BMXWAIST/BMXHT, Arm_length= BMXARML, METSIR = (log(Glucose*2+TAG)*BMXBMI)/log(HDL), EXTSUP_FAT = DXXLAFAT_N + DXXRAFAT_N, EXTINF_FAT = DXXLLFAT_N + DXXRLFAT_N,
         DXDHELE_N = (DXDHELE/1000), DXDLALE_M = (DXDLALE/1000), DXDRALE_N = (DXDRALE/1000), DXDRLLE_N = (DXDRLLE/1000), DXDLALE_N = (DXDLALE/1000),
         DXDLLLE_N = (DXDLLLE/1000), DXDTOLE_N=(DXDTOLE/1000),DXDTRLE_N = (DXDTRLE/1000)) %>%    
  mutate(METS_VF = 4.466 + 0.011*(log(METSIR)^3)+ 3.239*(log(ICE)^3)-0.319*(2-Sex) + 0.594*(log(Age)), EXTSUP_FAT_N = (EXTSUP_FAT/1000)/((Height/100)^2),
         EXTINF_FAT_N = (EXTINF_FAT/1000)/((Height/100)^2), EXTSUP_FAT_DXDTOFAT = (EXTSUP_FAT/1000)/(DXDTOFAT/1000), EXTINF_FAT_DXDTOFAT = (EXTINF_FAT/1000)/(DXDTOFAT/1000))


NHANES1$Sex <- factor(NHANES1$Sex, levels= c(1,2),labels= c("Men", "Women"))
NHANES1$Sex <- factor(NHANES1$Sex, levels= c("Men","Women"),labels= c("Men", "Women"))
NHANES1 <- merge(NHANES1,mortalidad,by="SEQN")
NHANES1 <- NHANES1%>%filter(Age>=20)

NHANES1$PhenoAge = 141.5 + ((log(-0.00553*log(1-(1-exp((-1.51714*exp(-19.907-0.0336*NHANES1$Albumin+0.0095*NHANES1$Cr +0.1953*NHANES1$Glucose+
                                                                       0.0954*log(NHANES1$CRP)-0.0120*NHANES1$LymP+0.0268*NHANES1$MCV+
                                                                       0.3306*NHANES1$RDW+0.00188*NHANES1$ALP+0.0554*NHANES1$WBC+
                                                                       0.0804*NHANES1$Age))/(0.0076927)))))))/(0.09165)

#Recode ethnicity
RecodeEth <- NHANES1$Ethnicity; NHANES1$Ethnicity[RecodeEth==5] <- 5
NHANES1$Ethnicity[RecodeEth==1] <- 3; NHANES1$Ethnicity[RecodeEth==2] <- 4
NHANES1$Ethnicity[RecodeEth==3] <- 1; NHANES1$Ethnicity[RecodeEth==4] <- 2

d1<-dummies::dummy(NHANES1$ucod_leading)
colnames(d1)<-c("Alive","Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                "Nephone_diseases","Other_causes")
NHANES1<-cbind(NHANES1, d1)

d1<-dummies::dummy(NHANES0$ucod_leading)
colnames(d1)<-c("Alive","Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                "Nephone_diseases","Other_causes")
NHANES0<-cbind(NHANES0, d1)

nhanes0<-NHANES0 %>% dplyr::select(SEQN, Age, Sex, Ethnicity, mortstat,permth_int,PhenoAge,BMI, Thigh_circumference, Arm_circumference,Triceps_skinfold,Subscapular_skinfold,Leg_length, Arm_length, Waist, Height, Weight,
                                   ICE,"Alive","Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                                   "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                                   "Nephone_diseases","Other_causes",Diabetes,Asthma,Arthritis,Heart_failure, 
                                   Heart_attack,Stroke,Emphysema,Bronchitis,Malignancy)
nhanes0$id<-rep(1, nrow(nhanes0))
NHANES1$Heart_failure<-NHANES1$Heart_Failure
nhanes1<-NHANES1 %>% dplyr::select(SEQN, Age, Sex, Ethnicity, mortstat,permth_int,PhenoAge,BMI, Thigh_circumference, Arm_circumference,Triceps_skinfold,Subscapular_skinfold,Leg_length, Arm_length, Waist, Height, Weight,
                                   ICE,"Alive","Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                                   "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                                   "Nephone_diseases","Other_causes",Diabetes,Asthma,Arthritis,Heart_failure, 
                                   Heart_attack,Stroke,Emphysema,Bronchitis,Malignancy)

nhanes1$id<-rep(2, nrow(nhanes1))
nhanes<-rbind(nhanes0, nhanes1)
nhanes$Diabetes<-ifelse(nhanes$Diabetes==2, 0, 1)
nhanes$Asthma<-ifelse(nhanes$Asthma==2, 0, 1)
nhanes$Arthritis<-ifelse(nhanes$Arthritis==2, 0, 1)
nhanes$Heart_failure<-ifelse(nhanes$Heart_failure==2, 0, 1)
nhanes$Heart_attack<-ifelse(nhanes$Heart_attack==2, 0, 1)
nhanes$Emphysema<-ifelse(nhanes$Emphysema==2, 0, 1)
nhanes$Bronchitis<-ifelse(nhanes$Bronchitis==2, 0, 1)
nhanes$Malignancy<-ifelse(nhanes$Malignancy==2, 0, 1)
nhanes$Stroke<-ifelse(nhanes$Stroke==2, 0, 1)
nhanes$num_comorb<-nhanes$Diabetes+nhanes$Asthma+nhanes$Arthritis+nhanes$Heart_failure+nhanes$Heart_attack+nhanes$Emphysema+nhanes$Bronchitis+nhanes$Malignancy+nhanes$Stroke

nhanes$tr_weight<-log(nhanes$Weight)
nhanes$tr_waist<-log(nhanes$Waist)
nhanes$tr_armc<-log(nhanes$Arm_circumference)
nhanes$tr_ice<-log(nhanes$ICE)
nhanes$tr_imc<-log(nhanes$BMI)
nhanes$tr_tric<-sqrt(nhanes$Triceps_skinfold)
nhanes$tr_subs<-sqrt(nhanes$Subscapular_skinfold)
nhanes$tr_arm<-log(nhanes$Arm_length)
nhanes$tr_thigh<-log(nhanes$Thigh_circumference)
nhanes$Ethnicity<-factor(nhanes$Ethnicity, labels = c("White", "Black", "Mexican-American", "Hispanic", "Other"))

nhanes_0<-nhanes %>% filter(!duplicated(SEQN))
nhanes0<-nhanes_0[,-c(7)] %>% drop_na()
nhanes0<-(merge(nhanes0,nhanes_0[,c(1,7)],by = "SEQN"))
nhanes0$PhenoAge[is.infinite(nhanes0$PhenoAge)]<-NA

apply(apply(nhanes0,2,is.na),2,sum)
apply(apply(nhanes0,2,is.infinite),2,sum)
nhanes_men<-nhanes0 %>% filter(Sex=="Men")
nhanes_women<-nhanes0 %>% filter(Sex=="Women")

####---- Descriptive statistics ----####

##Number of cases per dataset
table(nhanes0$id)
nrow(nhanes0)

## Sex
table(nhanes0$Sex)
table(nhanes0$Sex) %>% prop.table() %>% round(3)
table(nhanes0$Sex, nhanes0$id) %>% prop.table(2)
table(nhanes0$id, nhanes0$Sex) %>% prop.test()

##Age
hist(nhanes0$Age); ad.test(nhanes0$Age)
quantile(nhanes0$Age)
tapply(nhanes0$PhenoAge, nhanes0$id, quantile, na.rm=T)
wilcox.test(nhanes0$Age~nhanes0$id)

## PhenoAge
hist(nhanes0$PhenoAge); ad.test(nhanes0$PhenoAge)
quantile(nhanes0$PhenoAge, na.rm = T)
tapply(nhanes0$PhenoAge, nhanes0$id, quantile, na.rm=T)
wilcox.test(nhanes0$PhenoAge~nhanes0$id)

#Ethnicity
table(nhanes0$Ethnicity)
table(nhanes0$Ethnicity) %>% prop.table() %>% round(3)
table(nhanes0$Ethnicity, nhanes0$id) %>% prop.table(2)
table(nhanes0$id, nhanes0$Ethnicity) %>% chisq.test()

#BMI
tapply(nhanes_men$BMI, nhanes_men$id, quantile, na.rm=T)
wilcox.test(nhanes_men$BMI~nhanes_men$id)
tapply(nhanes_women$BMI, nhanes_women$id, quantile, na.rm=T)
wilcox.test(nhanes_women$BMI~nhanes_women$id)

#WHtR
tapply(nhanes_men$ICE, nhanes_men$id, quantile, na.rm=T)
wilcox.test(nhanes_men$ICE~nhanes_men$id)
tapply(nhanes_women$ICE, nhanes_women$id, quantile, na.rm=T)
wilcox.test(nhanes_women$ICE~nhanes_women$id)

#Thigh circumference
tapply(nhanes_men$Thigh_circumference, nhanes_men$id, quantile, na.rm=T)
wilcox.test(nhanes_men$Thigh_circumference~nhanes_men$id)
tapply(nhanes_women$Thigh_circumference, nhanes_women$id, quantile, na.rm=T)
wilcox.test(nhanes_women$Thigh_circumference~nhanes_women$id)

#Arm circumference
tapply(nhanes_men$Arm_circumference, nhanes_men$id, quantile, na.rm=T)
wilcox.test(nhanes_men$Arm_circumference~nhanes_men$id)
tapply(nhanes_women$Arm_circumference, nhanes_women$id, quantile, na.rm=T)
wilcox.test(nhanes_women$Arm_circumference~nhanes_women$id)

#Arm length
tapply(nhanes_men$Arm_length, nhanes_men$id, quantile, na.rm=T)
wilcox.test(nhanes_men$Arm_length~nhanes_men$id)
tapply(nhanes_women$Arm_length, nhanes_women$id, quantile, na.rm=T)
wilcox.test(nhanes_women$Arm_length~nhanes_women$id)

#Triceps skinfold
tapply(nhanes_men$Triceps_skinfold, nhanes_men$id, quantile, na.rm=T)
wilcox.test(nhanes_men$Triceps_skinfold~nhanes_men$id)
tapply(nhanes_women$Triceps_skinfold, nhanes_women$id, quantile, na.rm=T)
wilcox.test(nhanes_women$Triceps_skinfold~nhanes_women$id)

#Subscapular skinfold
tapply(nhanes_men$Subscapular_skinfold, nhanes_men$id, quantile, na.rm=T)
wilcox.test(nhanes_men$Subscapular_skinfold~nhanes_men$id)
tapply(nhanes_women$Subscapular_skinfold, nhanes_women$id, quantile, na.rm=T)
wilcox.test(nhanes_women$Subscapular_skinfold~nhanes_women$id)

## Deaths
table(nhanes0$mortstat)
table(nhanes0$mortstat) %>% prop.table()
table(nhanes0$mortstat, nhanes0$id)
table(nhanes0$mortstat, nhanes0$id) %>% prop.table(2)
table(nhanes0$id, nhanes0$mortstat) %>% chisq.test()

##TABLE 1
# Sex
Male0 <- table(nhanes0$Sex)[2]
pMale0 <- round(((table(nhanes0$Sex)%>%prop.table())[2])*100,1)
Male1 <- table((nhanes0%>% filter(id==1))$Sex)[2]
pMale1 <- round(((table((nhanes0%>% filter(id==1))$Sex)%>%prop.table())[2])*100,1)
Male2 <- table((nhanes0%>% filter(id==2))$Sex)[2]
pMale2 <- round(((table((nhanes0%>% filter(id==2))$Sex)%>%prop.table())[2])*100,1)
p1 <- format.pval(prop.test(x=c(Male0,Male1), n=c(nrow((nhanes0%>% filter(id==1))),nrow((nhanes0%>% filter(id==2)))))$p.value, eps = .001, digits = 3) 

# Age
Age0 <- round(summary(nhanes0$Age)[3],2)
Age0.1 <- round(summary(nhanes0$Age)[2],2)
Age0.3 <- round(summary(nhanes0$Age)[5],2)
Age1 <- round(summary((nhanes0%>% filter(id==1))$Age)[3],2)
Age1.1 <- round(summary((nhanes0%>% filter(id==1))$Age)[2],2)
Age1.3 <- round(summary((nhanes0%>% filter(id==1))$Age)[5],2)
Age2 <- round(summary((nhanes0%>% filter(id==2))$Age)[3],2)
Age2.1 <- round(summary((nhanes0%>% filter(id==2))$Age)[2],2)
Age2.3 <- round(summary((nhanes0%>% filter(id==2))$Age)[5],2)
p2 <- format.pval(wilcox.test(nhanes0$Age~nhanes0$id)$p.value, eps = .001, digits = 3) 

tapply(nhanes0$Age,nhanes0$id,summary)

# PhenoAge
PhenoAge0 <- round(summary(nhanes0$PhenoAge)[3],2)
PhenoAge0.1 <- round(summary(nhanes0$PhenoAge)[2],2)
PhenoAge0.3 <- round(summary(nhanes0$PhenoAge)[5],2)
PhenoAge1 <- round(summary((nhanes0%>% filter(id==1))$PhenoAge)[3],2)
PhenoAge1.1 <- round(summary((nhanes0%>% filter(id==1))$PhenoAge)[2],2)
PhenoAge1.3 <- round(summary((nhanes0%>% filter(id==1))$PhenoAge)[5],2)
PhenoAge2 <- round(summary((nhanes0%>% filter(id==2))$PhenoAge)[3],2)
PhenoAge2.1 <- round(summary((nhanes0%>% filter(id==2))$PhenoAge)[2],2)
PhenoAge2.3 <- round(summary((nhanes0%>% filter(id==2))$PhenoAge)[5],2)
p3<- format.pval(wilcox.test((nhanes0%>% filter(id==1))$PhenoAge,(nhanes0%>% filter(id==2))$PhenoAge)$p.value, eps = .001, digits = 3)

tapply(nhanes0$PhenoAge,nhanes0$id,summary)

# >= 1 comorb
nhanes0$comorb<-ifelse(nhanes0$num_comorb>=1, 1, 0)
Comorb0 <- table(nhanes0$comorb)[2]
pComorb0 <- round(((table(nhanes0$comorb)%>%prop.table())[2])*100,1)
Comorb1 <- table((nhanes0%>% filter(id==1))$comorb)[2]
pComorb1 <- round(((table((nhanes0%>% filter(id==1))$comorb)%>%prop.table())[2])*100,1)
Comorb2 <- table((nhanes0%>% filter(id==2))$comorb)[2]
pComorb2 <- round(((table((nhanes0%>% filter(id==2))$comorb)%>%prop.table())[2])*100,1)
p4 <- format.pval(prop.test(x=c(Comorb0,Comorb1), n=c(nrow((nhanes0%>% filter(id==1))),nrow((nhanes0%>% filter(id==2)))))$p.value, eps = .001, digits = 3) 

# Mortality
Mort0 <- table(nhanes0$mortstat)[2]
pMort0 <- round(((table(nhanes0$mortstat)%>%prop.table())[2])*100,1)
Mort1 <- table((nhanes0%>% filter(id==1))$mortstat)[2]
pMort1 <- round(((table((nhanes0%>% filter(id==1))$mortstat)%>%prop.table())[2])*100,1)
Mort2 <- table((nhanes0%>% filter(id==2))$mortstat)[2]
pMort2 <- round(((table((nhanes0%>% filter(id==2))$mortstat)%>%prop.table())[2])*100,1)
p5 <- format.pval(prop.test(x=c(Mort1,Mort2), n=c(nrow((nhanes0%>% filter(id==1))),nrow((nhanes0%>% filter(id==2)))))$p.value, eps = .001, digits = 3) 

#Follow-up time
fut0 <- round(summary(nhanes0$permth_int)[3],2)
fut0.1 <- round(summary(nhanes0$permth_int)[2],2)
fut0.3 <- round(summary(nhanes0$permth_int)[5],2)
fut1 <- round(summary((nhanes0%>% filter(id==1))$permth_int)[3],2)
fut1.1 <- round(summary((nhanes0%>% filter(id==1))$permth_int)[2],2)
fut1.3 <- round(summary((nhanes0%>% filter(id==1))$permth_int)[5],2)
fut2 <- round(summary((nhanes0%>% filter(id==2))$permth_int)[3],2)
fut2.1 <- round(summary((nhanes0%>% filter(id==2))$permth_int)[2],2)
fut2.3 <- round(summary((nhanes0%>% filter(id==2))$permth_int)[5],2)
p6<- format.pval(wilcox.test(nhanes0$permth_int, nhanes0$id)$p.value, eps = .001, digits = 3)

#Ethnicity
table(nhanes0$Ethnicity)
table(nhanes0$Ethnicity) %>% prop.table() %>% round(3)
table(nhanes0$Ethnicity, nhanes0$id) %>% prop.table(2)
table(nhanes0$id, nhanes0$Ethnicity) %>% chisq.test()

#BMI Women
bmi_w0 <- round(summary(nhanes_women$BMI)[3],2)
bmi_w0.1 <- round(summary(nhanes_women$BMI)[2],2)
bmi_w0.3 <- round(summary(nhanes_women$BMI)[5],2)
bmi_w1 <- round(summary((nhanes_women%>% filter(id==1))$BMI)[3],2)
bmi_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$BMI)[2],2)
bmi_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$BMI)[5],2)
bmi_w2 <- round(summary((nhanes_women%>% filter(id==2))$BMI)[3],2)
bmi_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$BMI)[2],2)
bmi_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$BMI)[5],2)
p7<- format.pval(wilcox.test(nhanes_women$BMI, nhanes_women$id)$p.value, eps = .001, digits = 3)

#BMI Men
bmi_m0 <- round(summary(nhanes_men$BMI)[3],2)
bmi_m0.1 <- round(summary(nhanes_men$BMI)[2],2)
bmi_m0.3 <- round(summary(nhanes_men$BMI)[5],2)
bmi_m1 <- round(summary((nhanes_men%>% filter(id==1))$BMI)[3],2)
bmi_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$BMI)[2],2)
bmi_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$BMI)[5],2)
bmi_m2 <- round(summary((nhanes_men%>% filter(id==2))$BMI)[3],2)
bmi_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$BMI)[2],2)
bmi_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$BMI)[5],2)
p8<- format.pval(wilcox.test(nhanes_men$BMI, nhanes_men$id)$p.value, eps = .001, digits = 3)

#ICE Women
ice_w0 <- round(summary(nhanes_women$ICE)[3],2)
ice_w0.1 <- round(summary(nhanes_women$ICE)[2],2)
ice_w0.3 <- round(summary(nhanes_women$ICE)[5],2)
ice_w1 <- round(summary((nhanes_women%>% filter(id==1))$ICE)[3],2)
ice_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$ICE)[2],2)
ice_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$ICE)[5],2)
ice_w2 <- round(summary((nhanes_women%>% filter(id==2))$ICE)[3],2)
ice_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$ICE)[2],2)
ice_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$ICE)[5],2)
p9<- format.pval(wilcox.test(nhanes_women$ICE, nhanes_women$id)$p.value, eps = .001, digits = 3)

#ICE Men
ice_m0 <- round(summary(nhanes_men$ICE)[3],2)
ice_m0.1 <- round(summary(nhanes_men$ICE)[2],2)
ice_m0.3 <- round(summary(nhanes_men$ICE)[5],2)
ice_m1 <- round(summary((nhanes_men%>% filter(id==1))$ICE)[3],2)
ice_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$ICE)[2],2)
ice_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$ICE)[5],2)
ice_m2 <- round(summary((nhanes_men%>% filter(id==2))$ICE)[3],2)
ice_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$ICE)[2],2)
ice_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$ICE)[5],2)
p10<- format.pval(wilcox.test(nhanes_men$ICE, nhanes_men$id)$p.value, eps = .001, digits = 3)

#Thigh circumference Women
thigh_w0 <- round(summary(nhanes_women$Thigh_circumference)[3],2)
thigh_w0.1 <- round(summary(nhanes_women$Thigh_circumference)[2],2)
thigh_w0.3 <- round(summary(nhanes_women$Thigh_circumference)[5],2)
thigh_w1 <- round(summary((nhanes_women%>% filter(id==1))$Thigh_circumference)[3],2)
thigh_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Thigh_circumference)[2],2)
thigh_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Thigh_circumference)[5],2)
thigh_w2 <- round(summary((nhanes_women%>% filter(id==2))$Thigh_circumference)[3],2)
thigh_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Thigh_circumference)[2],2)
thigh_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Thigh_circumference)[5],2)
p11<- format.pval(wilcox.test(nhanes_women$Thigh_circumference, nhanes_women$id)$p.value, eps = .001, digits = 3)

#Thigh circumference Men
thigh_m0 <- round(summary(nhanes_men$Thigh_circumference)[3],2)
thigh_m0.1 <- round(summary(nhanes_men$Thigh_circumference)[2],2)
thigh_m0.3 <- round(summary(nhanes_men$Thigh_circumference)[5],2)
thigh_m1 <- round(summary((nhanes_men%>% filter(id==1))$Thigh_circumference)[3],2)
thigh_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Thigh_circumference)[2],2)
thigh_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Thigh_circumference)[5],2)
thigh_m2 <- round(summary((nhanes_men%>% filter(id==2))$Thigh_circumference)[3],2)
thigh_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Thigh_circumference)[2],2)
thigh_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Thigh_circumference)[5],2)
p12<- format.pval(wilcox.test(nhanes_men$Thigh_circumference, nhanes_men$id)$p.value, eps = .001, digits = 3)

#Arm circumference Women
armc_w0 <- round(summary(nhanes_women$Arm_circumference)[3],2)
armc_w0.1 <- round(summary(nhanes_women$Arm_circumference)[2],2)
armc_w0.3 <- round(summary(nhanes_women$Arm_circumference)[5],2)
armc_w1 <- round(summary((nhanes_women%>% filter(id==1))$Arm_circumference)[3],2)
armc_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Arm_circumference)[2],2)
armc_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Arm_circumference)[5],2)
armc_w2 <- round(summary((nhanes_women%>% filter(id==2))$Arm_circumference)[3],2)
armc_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Arm_circumference)[2],2)
armc_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Arm_circumference)[5],2)
p13<- format.pval(wilcox.test(nhanes_women$Arm_circumference, nhanes_women$id)$p.value, eps = .001, digits = 3)

#Arm circumference Men
armc_m0 <- round(summary(nhanes_men$Arm_circumference)[3],2)
armc_m0.1 <- round(summary(nhanes_men$Arm_circumference)[2],2)
armc_m0.3 <- round(summary(nhanes_men$Arm_circumference)[5],2)
armc_m1 <- round(summary((nhanes_men%>% filter(id==1))$Arm_circumference)[3],2)
armc_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Arm_circumference)[2],2)
armc_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Arm_circumference)[5],2)
armc_m2 <- round(summary((nhanes_men%>% filter(id==2))$Arm_circumference)[3],2)
armc_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Arm_circumference)[2],2)
armc_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Arm_circumference)[5],2)
p14<- format.pval(wilcox.test(nhanes_men$Arm_circumference, nhanes_men$id)$p.value, eps = .001, digits = 3)

#Arm Length Women
arml_w0 <- round(summary(nhanes_women$Arm_length)[3],2)
arml_w0.1 <- round(summary(nhanes_women$Arm_length)[2],2)
arml_w0.3 <- round(summary(nhanes_women$Arm_length)[5],2)
arml_w1 <- round(summary((nhanes_women%>% filter(id==1))$Arm_length)[3],2)
arml_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Arm_length)[2],2)
arml_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Arm_length)[5],2)
arml_w2 <- round(summary((nhanes_women%>% filter(id==2))$Arm_length)[3],2)
arml_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Arm_length)[2],2)
arml_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Arm_length)[5],2)
p15<- format.pval(wilcox.test(nhanes_women$Arm_length, nhanes_women$id)$p.value, eps = .001, digits = 3)

#Arm Length Men
arml_m0 <- round(summary(nhanes_men$Arm_length)[3],2)
arml_m0.1 <- round(summary(nhanes_men$Arm_length)[2],2)
arml_m0.3 <- round(summary(nhanes_men$Arm_length)[5],2)
arml_m1 <- round(summary((nhanes_men%>% filter(id==1))$Arm_length)[3],2)
arml_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Arm_length)[2],2)
arml_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Arm_length)[5],2)
arml_m2 <- round(summary((nhanes_men%>% filter(id==2))$Arm_length)[3],2)
arml_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Arm_length)[2],2)
arml_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Arm_length)[5],2)
p16<- format.pval(wilcox.test(nhanes_men$Arm_length, nhanes_men$id)$p.value, eps = .001, digits = 3)


#Triceps skinfold Women
tric_w0 <- round(summary(nhanes_women$Triceps_skinfold)[3],2)
tric_w0.1 <- round(summary(nhanes_women$Triceps_skinfold)[2],2)
tric_w0.3 <- round(summary(nhanes_women$Triceps_skinfold)[5],2)
tric_w1 <- round(summary((nhanes_women%>% filter(id==1))$Triceps_skinfold)[3],2)
tric_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Triceps_skinfold)[2],2)
tric_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Triceps_skinfold)[5],2)
tric_w2 <- round(summary((nhanes_women%>% filter(id==2))$Triceps_skinfold)[3],2)
tric_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Triceps_skinfold)[2],2)
tric_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Triceps_skinfold)[5],2)
p17<- format.pval(wilcox.test(nhanes_women$Triceps_skinfold, nhanes_women$id)$p.value, eps = .001, digits = 3)

#Triceps skinfold Men
tric_m0 <- round(summary(nhanes_men$Triceps_skinfold)[3],2)
tric_m0.1 <- round(summary(nhanes_men$Triceps_skinfold)[2],2)
tric_m0.3 <- round(summary(nhanes_men$Triceps_skinfold)[5],2)
tric_m1 <- round(summary((nhanes_men%>% filter(id==1))$Triceps_skinfold)[3],2)
tric_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Triceps_skinfold)[2],2)
tric_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Triceps_skinfold)[5],2)
tric_m2 <- round(summary((nhanes_men%>% filter(id==2))$Triceps_skinfold)[3],2)
tric_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Triceps_skinfold)[2],2)
tric_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Triceps_skinfold)[5],2)
p18<- format.pval(wilcox.test(nhanes_men$Triceps_skinfold, nhanes_men$id)$p.value, eps = .001, digits = 3)


#Subscapular skinfold Women
subs_w0 <- round(summary(nhanes_women$Subscapular_skinfold)[3],2)
subs_w0.1 <- round(summary(nhanes_women$Subscapular_skinfold)[2],2)
subs_w0.3 <- round(summary(nhanes_women$Subscapular_skinfold)[5],2)
subs_w1 <- round(summary((nhanes_women%>% filter(id==1))$Subscapular_skinfold)[3],2)
subs_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Subscapular_skinfold)[2],2)
subs_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Subscapular_skinfold)[5],2)
subs_w2 <- round(summary((nhanes_women%>% filter(id==2))$Subscapular_skinfold)[3],2)
subs_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Subscapular_skinfold)[2],2)
subs_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Subscapular_skinfold)[5],2)
p19<- format.pval(wilcox.test(nhanes_women$Subscapular_skinfold, nhanes_women$id)$p.value, eps = .001, digits = 3)

#Triceps skinfold Men
subs_m0 <- round(summary(nhanes_men$Subscapular_skinfold)[3],2)
subs_m0.1 <- round(summary(nhanes_men$Subscapular_skinfold)[2],2)
subs_m0.3 <- round(summary(nhanes_men$Subscapular_skinfold)[5],2)
subs_m1 <- round(summary((nhanes_men%>% filter(id==1))$Subscapular_skinfold)[3],2)
subs_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Subscapular_skinfold)[2],2)
subs_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Subscapular_skinfold)[5],2)
subs_m2 <- round(summary((nhanes_men%>% filter(id==2))$Subscapular_skinfold)[3],2)
subs_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Subscapular_skinfold)[2],2)
subs_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Subscapular_skinfold)[5],2)
p20<- format.pval(wilcox.test(nhanes_men$Subscapular_skinfold, nhanes_men$id)$p.value, eps = .001, digits = 3)


tab0<-data.frame("Characteristics"=c("Female (%)","Age (years)", "PhenoAge (years)", "???1 comorbidity (%)", "Mortality (%)", "Follow-up (months)",
                                     "BMI Females (kg/m2)", "BMI Males (kg/m2)", "WHtR Females", "WHtR Males", "Thigh circumference Females (cm)","Thigh circumference Males (cm)",
                                     "Arm circumference Females (cm)", "Arm circumference Males (cm)","Arm length Females (cm)", "Arm length Males (cm)",
                                     "Triceps skinfold Females (mm)", "Triceps skinfold Females (mm)", "Subscapular skinfold Females (cm)", "Subscapcular skinfold Males (mm)"),
                 "Overall (n=18,930)"=c(paste0(Male0," (",pMale0,")"),paste0(Age0," (",Age0.1,"-",Age0.3,")"),paste0(PhenoAge0," (",PhenoAge0.1,"-",PhenoAge0.3,")"),
                                      paste0(Comorb0," (",pComorb0,")"),paste0(Mort0," (",pMort0,")"),paste0(fut0," (",fut0.1,"-",fut0.3,")"),
                                      paste0(bmi_w0," (",bmi_w0.1,"-",bmi_w0.3,")"),paste0(bmi_m0," (",bmi_m0.1,"-",bmi_m0.3,")"),
                                      paste0(ice_w0," (",ice_w0.1,"-",ice_w0.3,")"),paste0(ice_m0," (",ice_m0.1,"-",ice_m0.3,")"),paste0(thigh_w0," (",thigh_w0.1,"-",thigh_w0.3,")"),paste0(thigh_m0," (",thigh_m0.1,"-",thigh_m0.3,")"),
                                      paste0(armc_w0," (",armc_w0.1,"-",armc_w0.3,")"),paste0(armc_m0," (",armc_m0.1,"-",armc_m0.3,")"),paste0(arml_w0," (",arml_w0.1,"-",arml_w0.3,")"),paste0(arml_m0," (",arml_m0.1,"-",arml_m0.3,")"),
                                      paste0(subs_w0," (",subs_w0.1,"-",subs_w0.3,")"),paste0(subs_m0," (",subs_m0.1,"-",subs_m0.3,")"),paste0(subs_w0," (",subs_w0.1,"-",subs_w0.3,")"),paste0(subs_m0," (",subs_m0.1,"-",subs_m0.3,")")),
                 "NHANES-III (n=11,865)"=c(paste0(Male1," (",pMale1,")"),paste0(Age1," (",Age1.1,"-",Age1.3,")"),paste0(PhenoAge1," (",PhenoAge1.1,"-",PhenoAge1.3,")"),
                                            paste0(Comorb1," (",pComorb1,")"),paste0(Mort1," (",pMort1,")"),paste0(fut1," (",fut1.1,"-",fut1.3,")"),
                                            paste0(bmi_w1," (",bmi_w1.1,"-",bmi_w1.3,")"),paste0(bmi_m1," (",bmi_m1.1,"-",bmi_m1.3,")"),
                                            paste0(ice_w1," (",ice_w1.1,"-",ice_w1.3,")"),paste0(ice_m1," (",ice_m1.1,"-",ice_m1.3,")"),paste0(thigh_w1," (",thigh_w1.1,"-",thigh_w1.3,")"),paste0(thigh_m1," (",thigh_m1.1,"-",thigh_m1.3,")"),
                                            paste0(armc_w1," (",armc_w1.1,"-",armc_w1.3,")"),paste0(armc_m1," (",armc_m1.1,"-",armc_m1.3,")"),paste0(arml_w1," (",arml_w1.1,"-",arml_w1.3,")"),paste0(arml_m1," (",arml_m1.1,"-",arml_m1.3,")"),
                                            paste0(subs_w1," (",subs_w1.1,"-",subs_w1.3,")"),paste0(subs_m1," (",subs_m1.1,"-",subs_m1.3,")"),paste0(subs_w1," (",subs_w1.1,"-",subs_w1.3,")"),paste0(subs_m1," (",subs_m1.1,"-",subs_m1.3,")")),
                 "NHANES-IV (n=7,165)"=c(paste0(Male2," (",pMale2,")"),paste0(Age2," (",Age2.1,"-",Age2.3,")"),paste0(PhenoAge2," (",PhenoAge2.1,"-",PhenoAge2.3,")"),
                                         paste0(Comorb2," (",pComorb2,")"),paste0(Mort2," (",pMort2,")"),paste0(fut2," (",fut2.1,"-",fut2.3,")"),
                                         paste0(bmi_w2," (",bmi_w2.1,"-",bmi_w2.3,")"),paste0(bmi_m2," (",bmi_m2.1,"-",bmi_m2.3,")"),
                                         paste0(ice_w2," (",ice_w2.1,"-",ice_w2.3,")"),paste0(ice_m2," (",ice_m2.1,"-",ice_m2.3,")"),paste0(thigh_w2," (",thigh_w2.1,"-",thigh_w2.3,")"),paste0(thigh_m2," (",thigh_m2.1,"-",thigh_m2.3,")"),
                                         paste0(armc_w2," (",armc_w2.1,"-",armc_w2.3,")"),paste0(armc_m2," (",armc_m2.1,"-",armc_m2.3,")"),paste0(arml_w2," (",arml_w2.1,"-",arml_w2.3,")"),paste0(arml_m2," (",arml_m2.1,"-",arml_m2.3,")"),
                                         paste0(subs_w2," (",subs_w2.1,"-",subs_w2.3,")"),paste0(subs_m2," (",subs_m2.1,"-",subs_m2.3,")"),paste0(subs_w2," (",subs_w2.1,"-",subs_w2.3,")"),paste0(subs_m2," (",subs_m2.1,"-",subs_m2.3,")")),
                 "P-value"=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20))

tab0<-`names<-`(tab0,c("Characteristics","Overall (n=18,930)","NHANES-III (n=11,865)","NHANES-IV (n=7,165)","P-value"))
tab0<-align(flextable(tab0,cwidth = c(2,1.5,1.5,2,1)),align = "center",part = "all")
save_as_docx(tab0,path="tabla1.docx")


####---- Models for anthropometric measures ---####
### Models in men ###
gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+poly(BMI,3)+num_comorb+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk1<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+poly(ICE,2)+num_comorb+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk2<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Subscapular_skinfold,2)+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk3<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Triceps_skinfold,3)+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk4<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Thigh_circumference,4)+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk5<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+Arm_length+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk6<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Arm_circumference,2)+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk7<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+Leg_length+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk8<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+Height+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk9<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Weight,2)+shape(Ethnicity),dist="gompertz",data=nhanes_men)
p1F<-1-predict(gomp1, newdata = nhanes_men ,type="survival", ci=F, times = c(120))
nhanes_men$risk10<-p1F$.pred



### Models in women ###
gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(BMI,2)+shape(Ethnicity),dist="gompertz",data=nhanes_women) #Podria ser 3
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk1<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+ICE+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk2<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Subscapular_skinfold,2)+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk3<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Triceps_skinfold,2)+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk4<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Thigh_circumference,3)+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk5<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+Arm_length+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk6<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Arm_circumference,2)+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk7<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+Leg_length+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk8<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+Height+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk9<-p1F$.pred

gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(Weight,2)+shape(Ethnicity),dist="gompertz",data=nhanes_women)
p1F<-1-predict(gomp1, newdata = nhanes_women ,type="survival", ci=F, times = c(120))
nhanes_women$risk10<-p1F$.pred




### Run Figure 1 ###
nhanes_fin<-rbind(nhanes_women, nhanes_men)

f1A<-ggplot(nhanes_fin, aes(x=BMI, y=risk1, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+xlim(17,45)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Body-mass index (kg/m2)")+scale_fill_manual(values=c("#994455","#6699CC"))

f1B<-ggplot(nhanes_fin, aes(x=ICE, y=risk2, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Waist-to-height ratio")+xlim(0.39,0.8)+scale_fill_manual(values=c("#994455","#6699CC"))

f1C<-ggplot(nhanes_fin, aes(x=Subscapular_skinfold, y=risk3, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+xlim(5,45)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Subscapular skinfold (cm)")+scale_fill_manual(values=c("#994455","#6699CC"))

f1D<-ggplot(nhanes_fin, aes(x=Triceps_skinfold, y=risk4, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+xlim(4,43)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Triceps skinfold (cm)")+scale_fill_manual(values=c("#994455","#6699CC"))

f1E<-ggplot(nhanes_fin, aes(x=Thigh_circumference, y=risk5, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+xlim(36,71)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Thigh circumference (cm)")+scale_fill_manual(values=c("#994455","#6699CC"))

f1F<-ggplot(nhanes_fin, aes(x=Arm_circumference, y=risk7, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Arm circumference (cm)")+xlim(21,45)+scale_fill_manual(values=c("#994455","#6699CC"))

f1G<-ggplot(nhanes_fin, aes(x=Arm_length, y=risk6, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Arm length (cm)")+xlim(25,45)+scale_fill_manual(values=c("#994455","#6699CC"))

f1H<-ggplot(nhanes_fin, aes(x=Weight, y=risk10, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Weight (Kg)")+xlim(35,130)+scale_fill_manual(values=c("#994455","#6699CC"))

f1I<-ggplot(nhanes_fin, aes(x=Leg_length, y=risk8, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Leg length (cm)")+xlim(28,51)+scale_fill_manual(values=c("#994455","#6699CC"))

f1J<-ggplot(nhanes_fin, aes(x=Height, y=risk9, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+
  theme_pubclean()+ylab("Risk of 10-year mortality")+xlab("Height (cm)")+xlim(140,193)+scale_fill_manual(values=c("#994455","#6699CC"))

fig1<-ggarrange(f1A, f1B, f1H, f1G, f1C, f1D, f1E, f1F, ncol = 4, nrow=2, labels = letters[1:8], common.legend = TRUE)

ggsave(fig1,filename = "Figure1.jpg", 
       width = 35, 
       height = 20,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)



####---- Development of AnthropoAge ----####

### Training and validation ###
train_nhanes <- nhanes0[nhanes0$id==1,]
test_nhanes <- nhanes0[nhanes0$id==2,]

###Gompertz models###
options(scipen=10)
gomp1aM<-flexsurvreg(Surv(permth_int,mortstat)~Age+poly(tr_thigh,2)+tr_armc+BMI+
                       poly(tr_ice,2)+tr_arm+shape(Ethnicity),dist="gompertz",data=train_nhanes%>%filter(Sex=="Men"))
gomp1aM;BIC(gomp1aM)

gomp1aF<-flexsurvreg(Surv(permth_int,mortstat)~Age+tr_weight+poly(tr_thigh,2)+
                       tr_tric+tr_ice+poly(tr_subs,2)+
                       shape(Ethnicity),dist="gompertz",data=train_nhanes%>%filter(Sex=="Women"))
gomp1aF;BIC(gomp1aF)
forest
gomp1bM<-flexsurvreg(Surv(permth_int,mortstat)~Age,dist="gompertz",data=train_nhanes%>%filter(Sex=="Men"))
gomp1bM;BIC(gomp1bM)
sM<-1/((exp(coef(gomp1bM)[1]*120)-1)/((coef(gomp1bM)[1])))
b0M<-coef(gomp1bM)[2]
b1M<-coef(gomp1bM)[3]

gomp1bF<-flexsurvreg(Surv(permth_int,mortstat)~Age,dist="gompertz",data=train_nhanes%>%filter(Sex=="Women"))
gomp1bF;BIC(gomp1bF)
sW<-1/((exp(coef(gomp1bF)[1]*120)-1)/((coef(gomp1bF)[1])))
b0W<-coef(gomp1bF)[2]
b1W<-coef(gomp1bF)[3]


### Table 2 ###
coef1M<-round(coef(gomp1aM),4)
confint1M<-round(confint(gomp1aM)[,1],4)
confint2M<-round(confint(gomp1aM)[,2],4)
bic1<-round(BIC(gomp1aM),1)

coef1F<-round(coef(gomp1aF),4)
confint1F<-round(confint(gomp1aF)[,1],4)
confint2F<-round(confint(gomp1aF)[,2],4)
bic2<-round(BIC(gomp1aF),1)


tab2 <- data.frame("Model"=c(c(rep(" ",5)),"AnthropoAge Males",paste0("BIC ",bic1),c(rep(" ",13)),"AnthropoAge Females",paste0("BIC ",bic2),c(rep(" ",7))), 
                   "Parameter"=c("Shape", "Rate", "Chronological Age", "Thigh circumference (OP,1)","Thigh circumference (OP,2)","Arm circumference", "BMI", "WHtR (OP,1)", "WHtR (OP,2)",
                                 "Arm length", "Shape - Black ethnicity", "Shape - Mexican American ehtnicity", "Shape - Hispanic ethnicity", "Shape - Other ethnicity",
                                 " ", "Shape", "Rate", "Chronological Age", "Weight", "Thigh circumference (OP,1)","Thigh circumference (OP,2)", "Tricipital skinfold", "WHtR",
                                 "Subscapular skinfold (OP,1)", "Subscapular skinfold (OP,2)","Shape - Black ethnicity", "Shape - Mexican American ehtnicity", "Shape - Hispanic ethnicity", "Shape - Other ethnicity"),
                   "B-coefficient"=c(coef1M," ",coef1F),
                   "Lower 95%CI"=c(confint1M, " ",confint1F),
                   "Upper 95%CI"=c(confint2M, " ",confint1F))
rownames(tab2)<-NULL

tab2<-`names<-`(tab2,c("Model","Parameter","B-coefficient","Lower 95%CI","Upper 95%CI"))
tab2 <-align(flextable(tab2,cwidth = c(1.5,2.5,1,1,1)),align = "center",part = "all")

save_as_docx(tab2,path="tabla2.docx")


#### Training dataset####
p1F<-predict(gomp1aF, newdata =train_nhanes %>%filter(Sex=="Women"),type="survival", ci=F, times = c(120))
p1M<-predict(gomp1aM, newdata =train_nhanes %>%filter(Sex=="Men"),type="survival", ci=F, times = c(120))

train_nhanes$pred[train_nhanes$Sex=="Women"]<-as.numeric(1-p1F$.pred)
train_nhanes$pred[train_nhanes$Sex=="Men"]<-as.numeric(1-p1M$.pred)
train_nhanes$AnthropoAge[train_nhanes$Sex=="Women"]<-(log(-sW*log(1-train_nhanes$pred[train_nhanes$Sex=="Women"]))-b0W)/b1W
train_nhanes$AnthropoAge[train_nhanes$Sex=="Men"]<-(log(-sM*log(1-train_nhanes$pred[train_nhanes$Sex=="Men"]))-b0M)/b1M

### Accelerated metrics ###
train_nhanes$AnthropoAgeAccel <- NULL
m1F<-lm(AnthropoAge~Age, data=train_nhanes %>% filter(Sex=="Women"))
train_nhanes$AnthropoAgeAccel[train_nhanes$Sex=="Women"]<-m1F$residuals
m1M<-lm(AnthropoAge~Age, data=train_nhanes %>% filter(Sex=="Men"))
train_nhanes$AnthropoAgeAccel[train_nhanes$Sex=="Men"]<-m1M$residuals

train_nhanes$PhenoAgeAccel <- NULL
m1Ph<-lm(PhenoAge~Age, data=train_nhanes %>% filter(!is.infinite(PhenoAge)))
train_nhanes$PhenoAgeAccel[!is.na(train_nhanes$PhenoAge)]<-m1Ph$residuals

save(gomp1aM, file="Shiny App/anthropoage/Models/1_CAnthropo_M.rda")
save(gomp1aF, file="Shiny App/anthropoage/Models/2_CAnthropo_F.rda")
save(gomp1bM, file="Shiny App/anthropoage/Models/3_CAge_M.rda")
save(gomp1bF, file="Shiny App/anthropoage/Models/4_CAge_F.rda")
save(m1M, file="Shiny App/anthropoage/Models/5_CAccel_M.rda")
save(m1F, file="Shiny App/anthropoage/Models/6_CAccel_F.rda")
save(m1Ph, file="Shiny App/anthropoage/Models/13_PhenoAgeAccel.rda")

tapply(train_nhanes$AnthropoAge, train_nhanes$Sex, quantile)
wilcox.test(train_nhanes$AnthropoAge~as.numeric(train_nhanes$Sex))
tapply(train_nhanes$AnthropoAgeAccel, train_nhanes$Sex, quantile)
wilcox.test(train_nhanes$AnthropoAgeAccel~as.numeric(train_nhanes$Sex)) 
tapply(train_nhanes$PhenoAgeAccel, train_nhanes$Sex, quantile, na.rm=T)
wilcox.test(train_nhanes$PhenoAgeAccel~as.numeric(train_nhanes$Sex))


#### Testing dataset ####
p1F<-predict(gomp1aF, newdata = test_nhanes %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M<-predict(gomp1aM, newdata = test_nhanes %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

test_nhanes$pred[test_nhanes$Sex=="Women"]<-as.numeric(1-p1F$.pred)
test_nhanes$pred[test_nhanes$Sex=="Men"]<-as.numeric(1-p1M$.pred)
test_nhanes$AnthropoAge[test_nhanes$Sex=="Women"]<-(log(-sW*log(1-test_nhanes$pred[test_nhanes$Sex=="Women"]))-b0W)/b1W
test_nhanes$AnthropoAge[test_nhanes$Sex=="Men"]<-(log(-sM*log(1-test_nhanes$pred[test_nhanes$Sex=="Men"]))-b0M)/b1M

##Accelerated metrics ###
test_nhanes$AnthropoAgeAccel <- NULL
m1.1F<-lm(AnthropoAge~Age, data=test_nhanes %>% filter(Sex=="Women"))
test_nhanes$AnthropoAgeAccel[test_nhanes$Sex=="Women"]<-m1.1F$residuals
m1.1M<-lm(AnthropoAge~Age, data=test_nhanes %>% filter(Sex=="Men"))
test_nhanes$AnthropoAgeAccel[test_nhanes$Sex=="Men"]<-m1.1M$residuals

test_nhanes$PhenoAgeAccel <- NULL
m1.1Ph<-lm(PhenoAge~Age, data=test_nhanes %>% filter(!is.infinite(PhenoAge)))
test_nhanes$PhenoAgeAccel[!is.na(test_nhanes$PhenoAge)]<-m1.1Ph$residuals

tapply(test_nhanes$AnthropoAge, test_nhanes$Sex, quantile)
wilcox.test(test_nhanes$AnthropoAge~as.numeric(test_nhanes$Sex))
tapply(test_nhanes$AnthropoAgeAccel, test_nhanes$Sex, quantile)
wilcox.test(test_nhanes$AnthropoAgeAccel~as.numeric(test_nhanes$Sex))
tapply(test_nhanes$PhenoAgeAccel, test_nhanes$Sex, quantile, na.rm=T)
wilcox.test(test_nhanes$PhenoAgeAccel~as.numeric(test_nhanes$Sex))


##Figure 2
f1<-ggplot(train_nhanes, aes(x=AnthropoAge, y=Age, col=Sex))+geom_jitter(size=0.75,alpha=0.6)+geom_smooth(method="lm")+
  theme_classic()+theme(legend.position = "bottom")+scale_color_manual(values=c("#994455","#8DCBE4"))

f2<-ggplot(train_nhanes, aes(x=PhenoAge, y=Age, col=Sex))+geom_jitter(size=0.75,alpha=0.6)+geom_smooth(method="lm")+
  theme_classic()+theme(legend.position = "bottom")+scale_color_manual(values=c("#994455","#8DCBE4"))

f3<-ggplot(test_nhanes, aes(x=AnthropoAge, y=Age, col=Sex))+geom_jitter(size=0.75,alpha=0.6)+geom_smooth(method="lm")+
  theme_classic()+theme(legend.position = "bottom")+scale_color_manual(values=c("#994455","#8DCBE4"))

f4<-ggplot(test_nhanes, aes(x=PhenoAge, y=Age, col=Sex))+geom_jitter(size=0.75,alpha=0.6)+geom_smooth(method="lm")+
  theme_classic()+theme(legend.position = "bottom")+scale_color_manual(values=c("#994455","#8DCBE4"))

f5 <- ggplot(train_nhanes,aes(x=AnthropoAgeAccel, fill=Sex)) + scale_fill_manual(values=c("#994455","#8DCBE4"))+
  geom_density(alpha=0.75,size=0.3) + theme_classic()+ylab("Density")+xlab("AnthropoAgeAccel (years)")+theme(legend.position = "bottom")

f6 <- ggplot(train_nhanes,aes(x=PhenoAgeAccel, fill=Sex)) + scale_fill_manual(values=c("#994455","#8DCBE4"))+
  geom_density(alpha=0.75,size=0.3) + theme_classic()+ylab("Density")+xlab("PhenoAgeAccel (years)")+theme(legend.position = "bottom")

f7 <- ggplot(test_nhanes,aes(x=AnthropoAgeAccel, fill=Sex)) + scale_fill_manual(values=c("#994455","#8DCBE4"))+
  geom_density(alpha=0.75,size=0.3) + theme_classic()+ylab("Density")+xlab("AnthropoAgeAccel (years)")+theme(legend.position = "bottom")

f8 <- ggplot(test_nhanes,aes(x=PhenoAgeAccel, fill=Sex)) + scale_fill_manual(values=c("#994455","#8DCBE4"))+
  geom_density(alpha=0.75,size=0.3) + theme_classic()+ylab("Density")+xlab("PhenoAgeAccel (years)")+theme(legend.position = "bottom")


fig1<-ggarrange(annotate_figure(ggarrange(f1,f2,f5,f6, labels=letters[c(1,2,5,6)], common.legend = TRUE), top=text_grob("Training cohort", face = "bold", size=12)),
                annotate_figure(ggarrange(f3,f4,f7,f8, labels=letters[c(3,4,7,8)], common.legend = TRUE), top=text_grob("Validation cohort", face = "bold", size=12)), nrow=1, ncol=2)

ggsave(file = "Figure2.jpg", 
       fig1,
       bg = "transparent",
       width = 27, 
       height = 15,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)





####---- Simplified AnthropoAge (S-AnthropoAge) ----####
#Data management
nhanes_0<-nhanes %>% filter(!duplicated(SEQN))
nhanes0<-nhanes_0[,-c(7)] %>% drop_na()
nhanes0<-(merge(nhanes0,nhanes_0[,c(1,7)],by = "SEQN"))
nhanes0$PhenoAge[is.infinite(nhanes0$PhenoAge)]<-NA
nhanes0$Ethnicity<-factor(nhanes0$Ethnicity)

train_nhanes1 <- nhanes0[nhanes0$id==1,]
test_nhanes1 <- nhanes0[nhanes0$id==2,]

#Gompertz w/anthropometric variables
gomp1aM1<-flexsurvreg(Surv(permth_int,mortstat)~Age+poly(BMI,3)+
                        poly(ICE,2)+shape(Ethnicity),dist="gompertz",data=train_nhanes1%>%filter(Sex=="Men"))
gomp1aM1; BIC(gomp1aM1)

gomp1aF1<-flexsurvreg(Surv(permth_int,mortstat)~Age+poly(BMI,3)+
                        poly(ICE,2)+shape(Ethnicity),dist="gompertz",data=train_nhanes1%>%filter(Sex=="Women"))
gomp1aF1; BIC(gomp1aF1)

#Gompertz w/age
gomp1bM1<-flexsurvreg(Surv(permth_int,mortstat)~Age,dist="gompertz",data=train_nhanes1%>%filter(Sex=="Men"))
gomp1bM1;BIC(gomp1bM1)
sM1<-1/((exp(coef(gomp1bM1)[1]*120)-1)/((coef(gomp1bM1)[1])))
b0M1<-coef(gomp1bM1)[2]
b1M1<-coef(gomp1bM1)[3]

gomp1bF1<-flexsurvreg(Surv(permth_int,mortstat)~Age,dist="gompertz",data=train_nhanes1%>%filter(Sex=="Women"))
gomp1bF1;BIC(gomp1bF1)
sW1<-1/((exp(coef(gomp1bF1)[1]*120)-1)/((coef(gomp1bF1)[1])))
b0W1<-coef(gomp1bF1)[2]
b1W1<-coef(gomp1bF1)[3]


## Predictions and S-AnthropoAge calculation (train)
p1F1<-predict(gomp1aF1, newdata = train_nhanes1 %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M1<-predict(gomp1aM1, newdata = train_nhanes1 %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

train_nhanes1$pred[train_nhanes1$Sex=="Women"]<-as.numeric(1-p1F1$.pred)
train_nhanes1$pred[train_nhanes1$Sex=="Men"]<-as.numeric(1-p1M1$.pred)
train_nhanes1$AnthropoAge2[train_nhanes1$Sex=="Women"]<-(log(-sW1*log(1-train_nhanes1$pred[train_nhanes1$Sex=="Women"]))-b0W1)/b1W1
train_nhanes1$AnthropoAge2[train_nhanes1$Sex=="Men"]<-(log(-sM1*log(1-train_nhanes1$pred[train_nhanes1$Sex=="Men"]))-b0M1)/b1M1

## Predictions and S-AnthropoAge calculation (test)
p1F1<-predict(gomp1aF1, newdata = test_nhanes1 %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M1<-predict(gomp1aM1, newdata = test_nhanes1 %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

test_nhanes1$pred[test_nhanes1$Sex=="Women"]<-as.numeric(1-p1F1$.pred)
test_nhanes1$pred[test_nhanes1$Sex=="Men"]<-as.numeric(1-p1M1$.pred)
test_nhanes1$AnthropoAge2[test_nhanes1$Sex=="Women"]<-(log(-sW1*log(1-test_nhanes1$pred[test_nhanes1$Sex=="Women"]))-b0W1)/b1W1
test_nhanes1$AnthropoAge2[test_nhanes1$Sex=="Men"]<-(log(-sM1*log(1-test_nhanes1$pred[test_nhanes1$Sex=="Men"]))-b0M1)/b1M1

#S-AnthropoAgeAccel (train)
train_nhanes1$AnthropoAgeAccel2 <- NULL
m2F<-lm(AnthropoAge2~Age, data=train_nhanes1 %>% filter(Sex=="Women"))
train_nhanes1$AnthropoAgeAccel2[train_nhanes1$Sex=="Women"]<-m2F$residuals
m2M<-lm(AnthropoAge2~Age, data=train_nhanes1 %>% filter(Sex=="Men"))
train_nhanes1$AnthropoAgeAccel2[train_nhanes1$Sex=="Men"]<-m2M$residuals

#S-AnthropoAgeAccel (test)
test_nhanes1$AnthropoAgeAccel2 <- NULL
m2.1F<-lm(AnthropoAge2~Age, data=test_nhanes1 %>% filter(Sex=="Women"))
test_nhanes1$AnthropoAgeAccel2[test_nhanes1$Sex=="Women"]<-m2.1F$residuals
m2.1M<-lm(AnthropoAge2~Age, data=test_nhanes1 %>% filter(Sex=="Men"))
test_nhanes1$AnthropoAgeAccel2[test_nhanes1$Sex=="Men"]<-m2.1M$residuals

save(gomp1aM1, file="Shiny App/anthropoage/Models/7_SAnthropo_M.rda")
save(gomp1aF1, file="Shiny App/anthropoage/Models/8_SAnthropo_F.rda")
save(gomp1bM1, file="Shiny App/anthropoage/Models/9_SAge_M.rda")
save(gomp1bF1, file="Shiny App/anthropoage/Models/10_SAge_F.rda")
save(m2M, file="Shiny App/anthropoage/Models/11_SAccel_M.rda")
save(m2F, file="Shiny App/anthropoage/Models/12_SAccel_F.rda")

par(mfrow=c(1,2))
hist(train_nhanes1$AnthropoAgeAccel2); hist(test_nhanes1$AnthropoAgeAccel2)


#S-AnthropoAge performance
train_nhanes1 <- train_nhanes1 %>% filter(!is.na(PhenoAge))
test_nhanes1 <- test_nhanes1 %>% filter(!is.na(PhenoAge))

#Correlation with Age
cor(train_nhanes1$Age, train_nhanes1$PhenoAge)
cor(test_nhanes1$Age, test_nhanes1$PhenoAge)

cor(train_nhanes1$Age, train_nhanes1$AnthropoAge2)
cor(test_nhanes1$Age, test_nhanes1$AnthropoAge2)


#Train
r1 <- pROC::roc(train_nhanes1$mortstat, train_nhanes1$PhenoAge, ci=T)
pROC::roc(train_nhanes1$mortstat, train_nhanes1$PhenoAge, ci=T)$auc
pROC::roc(train_nhanes1$mortstat, train_nhanes1$PhenoAge, ci=T)$ci

r2 <- pROC::roc(train_nhanes1$mortstat, train_nhanes1$AnthropoAge2, ci=T) 
pROC::roc(train_nhanes1$mortstat, train_nhanes1$AnthropoAge2, ci=T)$auc
pROC::roc(train_nhanes1$mortstat, train_nhanes1$AnthropoAge2, ci=T)$ci

r3 <- pROC::roc(train_nhanes1$mortstat, train_nhanes1$Age, ci=T)
pROC::roc(train_nhanes1$mortstat, train_nhanes1$Age, ci=T)$auc
pROC::roc(train_nhanes1$mortstat, train_nhanes1$Age, ci=T)$ci

#roc.test(r1, r2, method = "boot"); roc.test(r2, r3, method = "boot")


#Test
r1 <- pROC::roc(test_nhanes1$mortstat, test_nhanes1$PhenoAge, ci=T)
pROC::roc(test_nhanes1$mortstat, test_nhanes1$PhenoAge, ci=T)$auc
pROC::roc(test_nhanes1$mortstat, test_nhanes1$PhenoAge, ci=T)$ci

r2 <- pROC::roc(test_nhanes1$mortstat, test_nhanes1$AnthropoAge2, ci=T) 
pROC::roc(test_nhanes1$mortstat, test_nhanes1$AnthropoAge2, ci=T)$auc
pROC::roc(test_nhanes1$mortstat, test_nhanes1$AnthropoAge2, ci=T)$ci

r3 <- pROC::roc(test_nhanes1$mortstat, test_nhanes1$Age, ci=T)
pROC::roc(test_nhanes1$mortstat, test_nhanes1$Age, ci=T)$auc
pROC::roc(test_nhanes1$mortstat, test_nhanes1$Age, ci=T)$ci

#roc.test(r1, r2, method = "boot"); roc.test(r2, r3, method = "boot")


s3a<-ggplot(train_nhanes, aes(x=AnthropoAge, y=Age, col=Ethnicity))+geom_point()+geom_smooth(method="lm")+theme_classic()+ggthemes::scale_color_colorblind()+ylab("Chronological Age (years)")
s3b<-ggplot(test_nhanes, aes(x=AnthropoAge, y=Age, col=Ethnicity))+geom_point()+geom_smooth(method="lm")+theme_classic()+ggthemes::scale_color_colorblind()+ylab("Chronological Age (years)")
s3c<-ggplot(train_nhanes, aes(x=PhenoAge, y=Age, col=Ethnicity))+geom_point()+geom_smooth(method="lm")+theme_classic()+ggthemes::scale_color_colorblind()+ylab("Chronological Age (years)")
s3d<-ggplot(test_nhanes, aes(x=PhenoAge, y=Age, col=Ethnicity))+geom_point()+geom_smooth(method="lm")+theme_classic()+ggthemes::scale_color_colorblind()+ylab("Chronological Age (years)")
s3e<-ggplot(train_nhanes1, aes(x=AnthropoAge2, y=Age, col=Ethnicity))+geom_point()+geom_smooth(method="lm")+theme_classic()+ggthemes::scale_color_colorblind()+xlab("S-AnthropoAge")+ylab("Chronological Age (years)")
s3f<-ggplot(test_nhanes1, aes(x=AnthropoAge2, y=Age, col=Ethnicity))+geom_point()+geom_smooth(method="lm")+theme_classic()+ggthemes::scale_color_colorblind()+xlab("S-AnthropoAge")+ylab("Chronological Age (years)")

supp3<-ggarrange(s3a,s3b,s3c,s3d,s3e,s3f, labels = letters[1:6], common.legend = TRUE)
ggsave(file = "SuppFig3.jpg", supp3, bg = "transparent",
       width = 45, height = 20, units=c("cm"),
       dpi = 500, limitsize = FALSE)



#### Bland-Altman ####

nhanes_0<-nhanes %>% filter(!duplicated(SEQN))
nhanes0<-nhanes_0[,-c(7)] %>% drop_na()
nhanes0<-(merge(nhanes0,nhanes_0[,c(1,7)],by = "SEQN"))
nhanes0$PhenoAge[is.infinite(nhanes0$PhenoAge)]<-NA

### Calculation of BA metrics ###
##AnthropoAge y AnthropoAgeAccel
p1F<-predict(gomp1aF, newdata = nhanes0 %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M<-predict(gomp1aM, newdata = nhanes0 %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

nhanes0$pred[nhanes0$Sex=="Women"]<-as.numeric(1-p1F$.pred)
nhanes0$pred[nhanes0$Sex=="Men"]<-as.numeric(1-p1M$.pred)
nhanes0$AnthropoAge[nhanes0$Sex=="Women"]<-(log(-sW*log(1-nhanes0$pred[nhanes0$Sex=="Women"]))-b0W)/b1W
nhanes0$AnthropoAge[nhanes0$Sex=="Men"]<-(log(-sM*log(1-nhanes0$pred[nhanes0$Sex=="Men"]))-b0M)/b1M

nhanes0$AnthropoAgeAccel[nhanes0$Sex=="Women"&nhanes0$id==1] <- (nhanes0%>%filter(Sex=="Women"&id==1))$AnthropoAge-predict(m1F,newdata=nhanes0%>%filter(Sex=="Women"&id==1))
nhanes0$AnthropoAgeAccel[nhanes0$Sex=="Men"&nhanes0$id==1] <- (nhanes0%>%filter(Sex=="Men"&id==1))$AnthropoAge-predict(m1M,newdata=nhanes0%>%filter(Sex=="Men"&id==1))
nhanes0$AnthropoAgeAccel[nhanes0$Sex=="Women"&nhanes0$id==2] <- (nhanes0%>%filter(Sex=="Women"&id==2))$AnthropoAge-predict(m1.1F,newdata=nhanes0%>%filter(Sex=="Women"&id==2))
nhanes0$AnthropoAgeAccel[nhanes0$Sex=="Men"&nhanes0$id==2] <- (nhanes0%>%filter(Sex=="Men"&id==2))$AnthropoAge-predict(m1.1M,newdata=nhanes0%>%filter(Sex=="Men"&id==2))


##S-AnthropoAge and S-AnthropoAgeAccel
p1F1<-predict(gomp1aF1, newdata = nhanes0 %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M1<-predict(gomp1aM1, newdata = nhanes0 %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

nhanes0$pred[nhanes0$Sex=="Women"]<-as.numeric(1-p1F1$.pred)
nhanes0$pred[nhanes0$Sex=="Men"]<-as.numeric(1-p1M1$.pred)
nhanes0$AnthropoAge2[nhanes0$Sex=="Women"]<-(log(-sW1*log(1-nhanes0$pred[nhanes0$Sex=="Women"]))-b0W1)/b1W1
nhanes0$AnthropoAge2[nhanes0$Sex=="Men"]<-(log(-sM1*log(1-nhanes0$pred[nhanes0$Sex=="Men"]))-b0M1)/b1M1

nhanes0$AnthropoAgeAccel2[nhanes0$Sex=="Women"&nhanes0$id==1] <- (nhanes0%>%filter(Sex=="Women"&id==1))$AnthropoAge2-predict(m2F,newdata=nhanes0%>%filter(Sex=="Women"&id==1))
nhanes0$AnthropoAgeAccel2[nhanes0$Sex=="Men"&nhanes0$id==1] <- (nhanes0%>%filter(Sex=="Men"&id==1))$AnthropoAge2-predict(m2M,newdata=nhanes0%>%filter(Sex=="Men"&id==1))
nhanes0$AnthropoAgeAccel2[nhanes0$Sex=="Women"&nhanes0$id==2] <- (nhanes0%>%filter(Sex=="Women"&id==2))$AnthropoAge2-predict(m2.1F,newdata=nhanes0%>%filter(Sex=="Women"&id==2))
nhanes0$AnthropoAgeAccel2[nhanes0$Sex=="Men"&nhanes0$id==2] <- (nhanes0%>%filter(Sex=="Men"&id==2))$AnthropoAge2-predict(m2.1M,newdata=nhanes0%>%filter(Sex=="Men"&id==2))


##PhenoAgeAccel
nhanes0 <- nhanes0 %>% filter(!is.na(PhenoAge)); nrow(nhanes0)
sum(is.na(nhanes0$PhenoAge)); sum(is.infinite(nhanes0$PhenoAge))
m1<-lm(PhenoAge~Age, data=nhanes0); nhanes0$PhenoAgeAccel<-m1$residuals

m1<-coxph(Surv(permth_int, mortstat) ~ AnthropoAgeAccel + PhenoAgeAccel + Age, data=nhanes0)
summary(m1); cox.zph(m1)
nhanes0$PhenoAgeAccel


#Bland-Altman AnthropoAge vs S-AnthropoAge
stats1<-blandr::blandr.statistics(nhanes0$AnthropoAge, nhanes0$AnthropoAge2, sig.level = 0.95)
s<-round(stats1$bias,3); sU<-round(stats1$biasUpperCI,3); sL<-round(stats1$biasLowerCI,3)
ann1 <- paste0("Bias = ",s," (",sL," - ",sU,")")

antro <- nhanes0%>%dplyr::select(AnthropoAge, AnthropoAge2)
ICC <- irr::icc(antro, model = "twoway", type = "agreement", unit = "single")
i1 <- round(ICC$value,4); i2 <- format.pval(ICC$p.value,eps="0.001")
ann2 <- paste0("ICC = ", i1, ", 95%CI 0.9928-0.9930, "," p", i2)

s1a <- BlandAltmanLeh::bland.altman.plot(nhanes0$AnthropoAge, nhanes0$AnthropoAge2, graph.sys = "ggplot2") %>% 
  ggedit::remove_geom(geom = "point") +
  geom_point(aes(color=nhanes0$Sex), alpha=0.6) + theme_pubclean() +
  scale_color_manual(name="Sex", values = c("#994455","#6699CC")) +
  labs(title="Bland-Altman Comparison\nAnthropoAge vs. S-AnthropoAge", x="Mean", y="Difference") +
  theme(plot.title = element_text(face="bold", hjust=0.5), legend.position = "bottom")+
  annotate("text",x=75,y=-27,label=ann1)

s1b <- ggplot(nhanes0, aes(x=AnthropoAge, y=AnthropoAge2, col=Sex))+geom_jitter(alpha=0.6)+
  theme_pubclean() + scale_color_manual(values=c("#994455","#6699CC"))+ylab("S-AnthropoAge (years)")+
  xlab("AnthropoAge (years)") + theme(legend.position = "bottom") + annotate("text", x=40, y=95, label=ann2)

supp2<-ggarrange(s1a, s1b, labels = c("a", "b"))
ggsave(file = "SuppFig2.jpg", supp2, bg = "transparent", width = 40, 
       height = 15, units=c("cm"), dpi = 500, limitsize = FALSE)


#### ROC curves ####
train_nhanes <- nhanes0 %>% filter(id==1); test_nhanes <- nhanes0 %>% filter(id==2)

### ROC curves in the training cohort (NHANES-III) ###
train_nhanes<-train_nhanes %>% filter(!is.na(PhenoAge) & !is.infinite(PhenoAge))

r1<-pROC::roc(train_nhanes$mortstat, train_nhanes$PhenoAge, ci=T)
r2<-pROC::roc(train_nhanes$mortstat, train_nhanes$AnthropoAge, ci=T)
r3<-pROC::roc(train_nhanes$mortstat, train_nhanes$Age, ci=T)
r4<-pROC::roc(train_nhanes$mortstat, train_nhanes$ICE, ci=T)
r5<-pROC::roc(train_nhanes$mortstat, train_nhanes$Thigh_circumference, ci=T)
r6<-pROC::roc(train_nhanes$mortstat, train_nhanes$Triceps_skinfold, ci=T)
r7<-pROC::roc(train_nhanes$mortstat, train_nhanes$BMI, ci=T)
r8<-pROC::roc(train_nhanes$mortstat, train_nhanes$AnthropoAge2, ci=T)

#roc.test(r1, r2, method = "boot"); roc.test(r3, r2, method = "boot")
#roc.test(r4, r2, method = "boot"); roc.test(r5, r2, method = "boot")
#roc.test(r6, r2, method = "boot"); roc.test(r7, r2, method = "boot")
#roc.test(r8, r1, method = "boot"); roc.test(r8, r2, method = "boot") #S-AnthropoAge

roc.list<-list(r1, r2, r3, r4, r5, r6, r7)

table1<-matrix(c(round(r1$ci[1],3), round(r1$auc,3), round(r1$ci[3],3),
                 round(r2$ci[1],3), round(r2$auc,3), round(r2$ci[3],3),
                 round(r3$ci[1],3), round(r3$auc,3), round(r3$ci[3],3),
                 round(r4$ci[1],3), round(r4$auc,3), round(r4$ci[3],3),
                 round(r5$ci[1],3), round(r5$auc,3), round(r5$ci[3],3),
                 round(r6$ci[1],3), round(r6$auc,3), round(r6$ci[3],3),
                 round(r7$ci[1],3), round(r7$auc,3), round(r7$ci[3],3)),ncol=3,byrow=T)

table1<-table1[,c(2,1,3)]
colnames(table1)<-c("AUC", "Lower-CI","Upper-CI")
rownames(table1)<-c("PhenoAge", "AnthropoAge", "CA", "WHtR", "Thigh circ.", "Tric. skinfold", "BMI")
names(roc.list)<-c("PhenoAge", "AnthropoAge", "CA", "WHtR", "Thigh circ.", "Tric. skinfold", "BMI")

fig3a<-ggroc(roc.list, size=1.5, aes("color","linetype"))+
  theme(legend.position="none")+ theme_minimal() + ggtitle("Training cohort") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid")+
  xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)") +
  labs(col="Scores", linetype="Scores")+ theme_pubclean()+
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size =17, face = 'plain'), axis.title.y = element_text(size =17, face = 'plain'),
        axis.text.x = element_text(size = 14, face ='plain'), axis.text.y = element_text(size = 14, face ='plain'),
        legend.text = element_text(size = 17))+
  annotation_custom(tableGrob(table1, theme=ttheme_default(base_size = 14)), xmin=-0.5, xmax=0, ymin=0.2, ymax=0.4)+
  scale_color_manual(values=c("#994455","#6699CC","#EECC66","#364B9A","#F67B4E","#EE99AA","#C2E4EF"))+
  theme(legend.key.size = unit(12,"mm"))+scale_linetype_manual(values=c("solid","11","dashed","longdash","dotted","twodash","solid"))+
  theme(text=element_text(size=17),axis.text=element_text(size=15),legend.text=element_text(size=17), legend.position = "bottom",
        plot.title = element_text(size=18, face="bold", hjust=0.5))


### ROC curves in the validation cohort (NHANES-IV) ###
test_nhanes<-test_nhanes %>% filter(!is.na(PhenoAge) & !is.infinite(PhenoAge))

r1<-pROC::roc(test_nhanes$mortstat, test_nhanes$PhenoAge, ci=T)
r2<-pROC::roc(test_nhanes$mortstat, test_nhanes$AnthropoAge, ci=T)
r3<-pROC::roc(test_nhanes$mortstat, test_nhanes$Age, ci=T)
r4<-pROC::roc(test_nhanes$mortstat, test_nhanes$ICE, ci=T)
r5<-pROC::roc(test_nhanes$mortstat, test_nhanes$Thigh_circumference, ci=T)
r6<-pROC::roc(test_nhanes$mortstat, test_nhanes$Triceps_skinfold, ci=T)
r7<-pROC::roc(test_nhanes$mortstat, test_nhanes$BMI, ci=T)
r8<-pROC::roc(test_nhanes$mortstat, test_nhanes$AnthropoAge2, ci=T)

#roc.test(r2, r1, method = "boot"); roc.test(r3, r2, method = "boot")
#roc.test(r4, r2, method = "boot"); roc.test(r5, r2, method = "boot")
#roc.test(r6, r2, method = "boot"); roc.test(r7, r2, method = "boot")
#roc.test(r8, r1, method = "boot"); roc.test(r8, r2, method = "boot") #S-AnthropoAge

roc.list<-list(r1, r2, r3, r4, r5, r6, r7)

table1<-matrix(c(round(r1$ci[1],3), round(r1$auc,3), round(r1$ci[3],3),
                 round(r2$ci[1],3), round(r2$auc,3), round(r2$ci[3],3),
                 round(r3$ci[1],3), round(r3$auc,3), round(r3$ci[3],3),
                 round(r4$ci[1],3), round(r4$auc,3), round(r4$ci[3],3),
                 round(r5$ci[1],3), round(r5$auc,3), round(r5$ci[3],3),
                 round(r6$ci[1],3), round(r6$auc,3), round(r6$ci[3],3),
                 round(r7$ci[1],3), round(r7$auc,3), round(r7$ci[3],3)),ncol=3,byrow=T)

table1<-table1[,c(2,1,3)]
colnames(table1)<-c("AUC", "Lower-CI","Upper-CI")
rownames(table1)<-c("PhenoAge", "AnthropoAge", "CA", "WHtR", "Thigh circ.", "Tric. skinfold", "BMI")
names(roc.list)<-c("PhenoAge", "AnthropoAge", "CA", "WHtR", "Thigh circ.", "Tric. skinfold", "BMI")

fig3b<-ggroc(roc.list, size=1.5, aes("color","linetype"))+
  theme(legend.position="none")+ theme_minimal() + ggtitle("Validation cohort") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid")+
  xlab("False Positive Rate (1-Specificity)") + 
  ylab("True Positive Rate (Sensitivity)") +
  labs(col="Scores", linetype="Scores")+
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size =17, face = 'plain'),
        axis.title.y = element_text(size =17, face = 'plain'),
        axis.text.x = element_text(size = 14, face ='plain'),
        axis.text.y = element_text(size = 14, face ='plain'),
        legend.text = element_text(size = 17))+theme_pubclean()+
  annotation_custom(tableGrob(table1, theme=ttheme_default(base_size = 14)), xmin=-0.5, xmax=0, ymin=0.2, ymax=0.4)+
  scale_color_manual(values=c("#994455","#6699CC","#EECC66","#364B9A","#F67B4E","#EE99AA","#C2E4EF"))+
  theme(legend.key.size = unit(12,"mm"))+scale_linetype_manual(values=c("solid","11","dashed","longdash","dotted","twodash","solid"))+
  theme(text=element_text(size=17),axis.text=element_text(size=15),legend.text=element_text(size=17), legend.position = "bottom",
        plot.title = element_text(size=18, face="bold", hjust=0.5))


f3.1<-ggarrange(fig3a," ",fig3b, nrow=1, ncol = 3,labels = c("a","","b"), widths=c(10,1,10), font.label = list("size"=17.5),
                common.legend = T, legend = "bottom")



#### ROC by sex and comorbidities ####
### ROC curves by sex ###

RA_Antropo <- OptimalCutpoints::optimal.cutpoints(X = AnthropoAge~mortstat, tag.healthy = 0, categorical.cov = "Sex",
                                                 methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                                 conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()

RA_Antropo2 <- OptimalCutpoints::optimal.cutpoints(X = AnthropoAge2~mortstat, tag.healthy = 0, categorical.cov = "Sex",
                                                  methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                                  conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()

RA_Pheno <- OptimalCutpoints::optimal.cutpoints(X = PhenoAge~mortstat, tag.healthy = 0, categorical.cov = "Sex",
                                               methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                               conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()


Sex1A <- RA_Antropo$p.table[["Men"]][["AUC_CI"]]; Sex1B <- RA_Antropo$p.table[["Women"]][["AUC_CI"]]
Sex2A <- RA_Antropo2$p.table[["Men"]][["AUC_CI"]]; Sex2B <- RA_Antropo2$p.table[["Women"]][["AUC_CI"]]
Sex3A <- RA_Pheno$p.table[["Men"]][["AUC_CI"]]; Sex3B <- RA_Pheno$p.table[["Women"]][["AUC_CI"]]


### ROC curves by sex ###
test_nhanes <- test_nhanes %>% mutate(comorb_cat = num_comorb %>% cut(c(-Inf,0,1, Inf))) %>%
  mutate(comorb_cat=factor(comorb_cat, labels = c("0", "1", ">=2")))

RB_Antropo <- OptimalCutpoints::optimal.cutpoints(X = AnthropoAge~mortstat, tag.healthy = 0, categorical.cov = "comorb_cat",
                                                  methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                                  conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()

RB_Antropo2 <- OptimalCutpoints::optimal.cutpoints(X = AnthropoAge2~mortstat, tag.healthy = 0, categorical.cov = "comorb_cat",
                                                   methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                                   conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()

RB_Pheno <- OptimalCutpoints::optimal.cutpoints(X = PhenoAge~mortstat, tag.healthy = 0, categorical.cov = "comorb_cat",
                                                methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                                conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()


Com1A <- RB_Antropo$p.table[["0"]][["AUC_CI"]]
Com1B <- RB_Antropo$p.table[["1"]][["AUC_CI"]]
Com1C <- RB_Antropo$p.table[[">=2"]][["AUC_CI"]]

Com2A <- RB_Antropo2$p.table[["0"]][["AUC_CI"]]
Com2B <- RB_Antropo2$p.table[["1"]][["AUC_CI"]]
Com2C <- RB_Antropo2$p.table[[">=2"]][["AUC_CI"]]

Com3A <- RB_Pheno$p.table[["0"]][["AUC_CI"]]
Com3B <- RB_Pheno$p.table[["1"]][["AUC_CI"]]
Com3C <- RB_Pheno$p.table[[">=2"]][["AUC_CI"]]


#Table
StratVar <- c("Sex","Men","Women","Comorbidities","0","1",">=2")
AnthropoVar <- c("",Sex1A,Sex1B,"",Com1A,Com1B,Com1C)
Anthropo2Var <- c("",Sex2A,Sex2B,"",Com2A,Com2B,Com2C)
PhenoVar <- c("",Sex3A,Sex3B,"",Com3A,Com3B,Com3C)

tabX <- data.frame("Strata"=StratVar,"AnthropoAge"=AnthropoVar,"S-AnthropoAge"=Anthropo2Var,"PhenoAge"=PhenoVar)
tableX <- align(flextable(tabX,cwidth = c(2,1.5,1.5,2,1)),align = "center",part = "all") %>% autofit()
save_as_docx(tableX,path="tablaX.docx")


##Figure
ROC_Sex1A <- test_nhanes%>%filter(Sex=="Men") %>% pROC::roc(mortstat, AnthropoAge, ci=T)
ROC_Sex1B <- test_nhanes%>%filter(Sex=="Women") %>% pROC::roc(mortstat, AnthropoAge, ci=T)
ROC_Com1A <- test_nhanes%>%filter(comorb_cat=="0") %>% pROC::roc(mortstat, AnthropoAge, ci=T)
ROC_Com1B <- test_nhanes%>%filter(comorb_cat=="1") %>% pROC::roc(mortstat, AnthropoAge, ci=T)
ROC_Com1C <- test_nhanes%>%filter(comorb_cat==">=2") %>% pROC::roc(mortstat, AnthropoAge, ci=T)

ROC_Sex2A <- test_nhanes%>%filter(Sex=="Men") %>% pROC::roc(mortstat, AnthropoAge2, ci=T)
ROC_Sex2B <- test_nhanes%>%filter(Sex=="Women") %>% pROC::roc(mortstat, AnthropoAge2, ci=T)
ROC_Com2A <- test_nhanes%>%filter(comorb_cat=="0") %>% pROC::roc(mortstat, AnthropoAge2, ci=T)
ROC_Com2B <- test_nhanes%>%filter(comorb_cat=="1") %>% pROC::roc(mortstat, AnthropoAge2, ci=T)
ROC_Com2C <- test_nhanes%>%filter(comorb_cat==">=2") %>% pROC::roc(mortstat, AnthropoAge2, ci=T)

ROC_Sex3A <- test_nhanes%>%filter(Sex=="Men") %>% pROC::roc(mortstat, PhenoAge, ci=T)
ROC_Sex3B <- test_nhanes%>%filter(Sex=="Women") %>% pROC::roc(mortstat, PhenoAge, ci=T)
ROC_Com3A <- test_nhanes%>%filter(comorb_cat=="0") %>% pROC::roc(mortstat, PhenoAge, ci=T)
ROC_Com3B <- test_nhanes%>%filter(comorb_cat=="1") %>% pROC::roc(mortstat, PhenoAge, ci=T)
ROC_Com3C <- test_nhanes%>%filter(comorb_cat==">=2") %>% pROC::roc(mortstat, PhenoAge, ci=T)

#AnthropoAge vs PhenoAge
#roc.test(ROC_Sex1A, ROC_Sex3A, method = "boot") #Males
#roc.test(ROC_Sex1B, ROC_Sex3B, method = "boot") #Females
#roc.test(ROC_Com1A, ROC_Com3A, method = "boot") #0 comorb
#roc.test(ROC_Com1B, ROC_Com3B, method = "boot") #1 comorb
#roc.test(ROC_Com1C, ROC_Com3C, method = "boot") #>=2 comorb

#S-AnthropoAge vs PhenoAge
#roc.test(ROC_Sex2A, ROC_Sex3A, method = "boot") #Males
#roc.test(ROC_Sex2B, ROC_Sex3B, method = "boot") #Females
#roc.test(ROC_Com2A, ROC_Com3A, method = "boot") #0 comorb
#roc.test(ROC_Com2B, ROC_Com3B, method = "boot") #1 comorb
#roc.test(ROC_Com2C, ROC_Com3C, method = "boot") #>=2 comorb


fig3c <- data.frame(rep(c("1A","2SA","3P"),2), c(rep("1Male",3),rep("2Female",3)),
                    c(ROC_Sex3A$ci[2],ROC_Sex1A$ci[2],ROC_Sex2A$ci[2],ROC_Sex3B$ci[2],ROC_Sex1B$ci[2],ROC_Sex2B$ci[2]),
                    c(ROC_Sex3A$ci[1],ROC_Sex1A$ci[1],ROC_Sex2A$ci[1],ROC_Sex3B$ci[1],ROC_Sex1B$ci[1],ROC_Sex2B$ci[1]),
                    c(ROC_Sex3A$ci[3],ROC_Sex1A$ci[3],ROC_Sex2A$ci[3],ROC_Sex3B$ci[3],ROC_Sex1B$ci[3],ROC_Sex2B$ci[3])) %>% 
  `names<-`(c("Variable","Sex","AUROC","Low","Up")) %>% 
  ggplot(aes(x=Sex, y=AUROC, fill=Variable)) +
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2, size=1.25, position=position_dodge(.9))+
  geom_col(position = "dodge", color="black", size=1.25, alpha=0.6) + 
  coord_cartesian(ylim=c(0.75,0.92)) + theme_pubclean() + 
  scale_fill_manual(values=c("#994455","#6699CC","#EECC66"), labels=c("PhenoAge", "AnthropoAge", "S-AnthropoAge"))+
  scale_x_discrete(labels=c("Male","Female")) + ggtitle("AUROC by sex\n(validation cohort)")+theme(legend.key.size = unit(12,"mm"))+
  theme(text=element_text(size=17),axis.text=element_text(size=15),legend.text=element_text(size=17), legend.position = "bottom",
        plot.title = element_text(size=18, face="bold", hjust=0.5))


fig3d <- data.frame(rep(c("1A","2SA","3P"),3), c(rep("A0",3),rep("B1",3),rep("C>=2",3)),
                    c(ROC_Com3A$ci[2],ROC_Com1A$ci[2],ROC_Com2A$ci[2],ROC_Com3B$ci[2],ROC_Com1B$ci[2],ROC_Com2B$ci[2],
                      ROC_Com3C$ci[2],ROC_Com1C$ci[2],ROC_Com2C$ci[2]),
                    c(ROC_Com3A$ci[1],ROC_Com1A$ci[1],ROC_Com2A$ci[1],ROC_Com3B$ci[1],ROC_Com1B$ci[1],ROC_Com2B$ci[1],
                      ROC_Com3C$ci[1],ROC_Com1C$ci[1],ROC_Com2C$ci[1]),
                    c(ROC_Com3A$ci[3],ROC_Com1A$ci[3],ROC_Com2A$ci[3],ROC_Com3B$ci[3],ROC_Com1B$ci[3],ROC_Com2B$ci[3],
                      ROC_Com3C$ci[3],ROC_Com1C$ci[3],ROC_Com2C$ci[3])) %>% 
  `names<-`(c("Variable","Comorbitidies","AUROC","Low","Up")) %>% 
  ggplot(aes(x=Comorbitidies, y=AUROC, fill=Variable)) +
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2, size=1.25, position=position_dodge(.9))+
  geom_col(position = "dodge", color="black", size=1.25, alpha=0.6) + 
  coord_cartesian(ylim=c(0.75,0.92)) + theme_pubclean() + 
  scale_fill_manual(values=c("#994455","#6699CC","#EECC66"), labels=c("PhenoAge", "AnthropoAge", "S-AnthropoAge"))+
  scale_x_discrete(labels=c("0","1", ">=2")) + ggtitle("AUROC by number of comorbidities\n(validation cohort)")+theme(legend.key.size = unit(12,"mm"))+
  theme(text=element_text(size=17),axis.text=element_text(size=15),legend.text=element_text(size=17), legend.position = "bottom",
        plot.title = element_text(size=18, face="bold", hjust=0.5))


f3.2<-ggarrange(fig3c," ",fig3d, nrow=1, ncol = 3,labels = c("c","","d"), widths=c(10,1,10), font.label = list("size"=17.5),
                common.legend = T, legend = "bottom")

f3 <- ggarrange(f3.1,"",f3.2, nrow=3, ncol=1, heights = c(10,1,7))


ggsave(file = "Figure3.jpg", f3, bg = "transparent",
       width = 40, height = 39*0.8, units=c("cm"),
       dpi = 300, limitsize = FALSE)




#### Cause-specific mortality ####
### Cardiovascular mortality ###
nhanes0$cv_mort <- NULL
nhanes0$cv_mort[nhanes0$Heart_diseases==1 & nhanes0$mortstat==1]<-1
nhanes0$cv_mort[nhanes0$Heart_diseases==0 & nhanes0$mortstat==1]<-2
nhanes0$cv_mort[nhanes0$Heart_diseases==0 & nhanes0$mortstat==0]<-0
nhanes0$cv_mort<-factor(nhanes0$cv_mort, labels = c("Censored", "CV mortality", "Other causes"))
fgdata_cv <- finegray(Surv(permth_int, cv_mort) ~ ., data=nhanes0, na.action=na.pass)

m1_cv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cv)
m2_cv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cv)
m3_cv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cv)

C1.1<-conc95(m1_cv)
C1.2<-conc95(m2_cv)
C1.3<-conc95(m3_cv)

B1.1<-BIC(m1_cv)-BIC(m2_cv)
B1.2<-BIC(m1_cv)-BIC(m3_cv)
B1.3<-BIC(m2_cv)-BIC(m3_cv)


### Diabetes mortality ###
nhanes0$db_mort <- NULL
nhanes0$db_mort[nhanes0$Diabetes_mellitus==1 & nhanes0$mortstat==1]<-1
nhanes0$db_mort[nhanes0$Diabetes_mellitus==0 & nhanes0$mortstat==1]<-2
nhanes0$db_mort[nhanes0$Diabetes_mellitus==0 & nhanes0$mortstat==0]<-0
nhanes0$db_mort<-factor(nhanes0$db_mort, labels = c("Censored", "Diabetes mortality", "Other causes"))

fgdata_db <- finegray(Surv(permth_int, db_mort) ~ ., data=nhanes0, na.action=na.pass)

m1_db <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_db)
m2_db <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_db)
m3_db <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_db)

C2.1<-conc95(m1_db)
C2.2<-conc95(m2_db)
C2.3<-conc95(m3_db)

B2.1<-BIC(m1_db)-BIC(m2_db)
B2.2<-BIC(m1_db)-BIC(m3_db)
B2.3<-BIC(m2_db)-BIC(m3_db)


### Stroke mortality ###
nhanes0$cer_mort <- NULL
nhanes0$cer_mort[nhanes0$Cerebrovascular_diseases==1 & nhanes0$mortstat==1]<-1
nhanes0$cer_mort[nhanes0$Cerebrovascular_diseases==0 & nhanes0$mortstat==1]<-2
nhanes0$cer_mort[nhanes0$Cerebrovascular_diseases==0 & nhanes0$mortstat==0]<-0
nhanes0$cer_mort<-factor(nhanes0$cer_mort, labels = c("Censored", "Diabetes mortality", "Other causes"))

fgdata_cer <- finegray(Surv(permth_int, cer_mort) ~ ., data=nhanes0, na.action=na.pass)

m1_cer <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cer)
m2_cer <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cer)
m3_cer <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cer)

C3.1<-conc95(m1_cer)
C3.2<-conc95(m2_cer)
C3.3<-conc95(m3_cer)

B3.1<-BIC(m1_cer)-BIC(m2_cer)
B3.2<-BIC(m1_cer)-BIC(m3_cer)
B3.3<-BIC(m2_cer)-BIC(m3_cer)


#Cancer-related mortality
nhanes0$can_mort <- NULL
nhanes0$can_mort[nhanes0$Malignant_neoplasms==1 & nhanes0$mortstat==1]<-1
nhanes0$can_mort[nhanes0$Malignant_neoplasms==0 & nhanes0$mortstat==1]<-2
nhanes0$can_mort[nhanes0$Malignant_neoplasms==0 & nhanes0$mortstat==0]<-0
nhanes0$can_mort<-factor(nhanes0$can_mort, labels = c("Censored", "Diabetes mortality", "Other causes"))

fgdata_can <- finegray(Surv(permth_int, can_mort) ~ ., data=nhanes0, na.action=na.pass)

m1_can <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_can)
m2_can <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_can)
m3_can <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_can)

C4.1<-conc95(m1_can)
C4.2<-conc95(m2_can)
C4.3<-conc95(m3_can)

B4.1<-BIC(m1_can)-BIC(m2_can)
B4.2<-BIC(m1_can)-BIC(m3_can)
B4.3<-BIC(m2_can)-BIC(m3_can)


#Mortalidad por Influenza/Pneumonia
nhanes0$inf <- NULL
nhanes0$inf[nhanes0$Influenza_or_pneumonia==1 & nhanes0$mortstat==1]<-1
nhanes0$inf[nhanes0$Influenza_or_pneumonia==0 & nhanes0$mortstat==1]<-2
nhanes0$inf[nhanes0$Influenza_or_pneumonia==0 & nhanes0$mortstat==0]<-0
nhanes0$inf<-factor(nhanes0$inf, labels = c("Censored", "Influenza", "Other causes"))

fgdata_inf <- finegray(Surv(permth_int, inf) ~ ., data=nhanes0, na.action=na.pass)

m1_inf <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_inf)
m2_inf <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_inf)
m3_inf <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_inf)

C5.1<-conc95(m1_inf)
C5.2<-conc95(m2_inf)
C5.3<-conc95(m3_inf)

B5.1<-BIC(m1_inf)-BIC(m2_inf)
B5.2<-BIC(m1_inf)-BIC(m3_inf)
B5.3<-BIC(m2_inf)-BIC(m3_inf)


#Mortalidad por ERC
nhanes0$ckd <- NULL
nhanes0$ckd[nhanes0$Nephone_diseases==1 & nhanes0$mortstat==1]<-1
nhanes0$ckd[nhanes0$Nephone_diseases==0 & nhanes0$mortstat==1]<-2
nhanes0$ckd[nhanes0$Nephone_diseases==0 & nhanes0$mortstat==0]<-0
nhanes0$ckd<-factor(nhanes0$ckd, labels = c("Censored", "CKD", "Other causes"))

fgdata_ckd <- finegray(Surv(permth_int, ckd) ~ ., data=nhanes0, na.action=na.pass)

m1_ckd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_ckd)
m2_ckd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_ckd)
m3_ckd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_ckd)

C6.1<-conc95(m1_ckd)
C6.2<-conc95(m2_ckd)
C6.3<-conc95(m3_ckd)

B6.1<-BIC(m1_ckd)-BIC(m2_ckd)
B6.2<-BIC(m1_ckd)-BIC(m3_ckd)
B6.3<-BIC(m2_ckd)-BIC(m3_ckd)


#Mortalidad por Alzheimer
nhanes0$alz <- NULL
nhanes0$alz[nhanes0$Alzheimer_disease==1 & nhanes0$mortstat==1]<-1
nhanes0$alz[nhanes0$Alzheimer_disease==0 & nhanes0$mortstat==1]<-2
nhanes0$alz[nhanes0$Alzheimer_disease==0 & nhanes0$mortstat==0]<-0
nhanes0$alz<-factor(nhanes0$alz, labels = c("Censored", "Alz", "Other causes"))

fgdata_alz <- finegray(Surv(permth_int, alz) ~ ., data=nhanes0, na.action=na.pass)

m1_alz <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_alz)
m2_alz <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_alz)
m3_alz <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_alz)

C7.1<-conc95(m1_alz)
C7.2<-conc95(m2_alz)
C7.3<-conc95(m3_alz)

B7.1<-BIC(m1_alz)-BIC(m2_alz)
B7.2<-BIC(m1_alz)-BIC(m3_alz)
B7.3<-BIC(m2_alz)-BIC(m3_alz)


#Mortalidad por epoc
nhanes0$copd <- NULL
nhanes0$copd[nhanes0$Chronic_lower_respiratory_diseases==1 & nhanes0$mortstat==1]<-1
nhanes0$copd[nhanes0$Chronic_lower_respiratory_diseases==0 & nhanes0$mortstat==1]<-2
nhanes0$copd[nhanes0$Chronic_lower_respiratory_diseases==0 & nhanes0$mortstat==0]<-0
nhanes0$copd<-factor(nhanes0$copd, labels = c("Censored", "COPD", "Other causes"))

fgdata_copd <- finegray(Surv(permth_int, copd) ~ ., data=nhanes0, na.action=na.pass)

m1_copd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_copd)
m2_copd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_copd)
m3_copd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge+Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_copd)

C8.1<-conc95(m1_copd)
C8.2<-conc95(m2_copd)
C8.3<-conc95(m3_copd)

B8.1<-BIC(m1_copd)-BIC(m2_copd)
B8.2<-BIC(m1_copd)-BIC(m3_copd)
B8.3<-BIC(m2_copd)-BIC(m3_copd)


### TABLE 3 ###

antro<-c(round(C1.1[1],3), round(C2.1[1],3), round(C3.1[1],3), round(C4.1[1],3), 
         round(C5.1[1],3), round(C6.1[1],3), round(C7.1[1],3), round(C8.1[1],3))
antro_s<-c(round(C1.2[1],3), round(C2.2[1],3), round(C3.2[1],3), round(C4.2[1],3), 
           round(C5.2[1],3), round(C6.2[1],3), round(C7.2[1],3), round(C8.2[1],3))
pheno<-c(round(C1.3[1],3), round(C2.3[1],3), round(C3.3[1],3), round(C4.3[1],3), 
         round(C5.3[1],3), round(C6.3[1],3), round(C7.3[1],3), round(C8.3[1],3))
antro_ci<-c(paste0(round(C1.1[2],3),"-",round(C1.1[3],3)),paste0(round(C2.1[2],3),"-",round(C2.1[3],3)),
            paste0(round(C3.1[2],3),"-",round(C3.1[3],3)),paste0(round(C4.1[2],3),"-",round(C4.1[3],3)),
            paste0(round(C5.1[2],3),"-",round(C5.1[3],3)),paste0(round(C6.1[2],3),"-",round(C6.1[3],3)),
            paste0(round(C7.1[2],3),"-",round(C7.1[3],3)),paste0(round(C8.1[2],3),"-",round(C8.1[3],3)))
antros_ci<-c(paste0(round(C1.2[2],3),"-",round(C1.2[3],3)),paste0(round(C2.2[2],3),"-",round(C2.2[3],3)),
            paste0(round(C3.2[2],3),"-",round(C3.2[3],3)),paste0(round(C4.2[2],3),"-",round(C4.2[3],3)),
            paste0(round(C5.2[2],3),"-",round(C5.2[3],3)),paste0(round(C6.2[2],3),"-",round(C6.2[3],3)),
            paste0(round(C7.2[2],3),"-",round(C7.2[3],3)),paste0(round(C8.2[2],3),"-",round(C8.2[3],3)))
pheno_ci<-c(paste0(round(C1.3[2],3),"-",round(C1.3[3],3)),paste0(round(C2.3[2],3),"-",round(C2.3[3],3)),
            paste0(round(C3.3[2],3),"-",round(C3.3[3],3)),paste0(round(C4.3[2],3),"-",round(C4.3[3],3)),
            paste0(round(C5.3[2],3),"-",round(C5.3[3],3)),paste0(round(C6.3[2],3),"-",round(C6.3[3],3)),
            paste0(round(C7.3[2],3),"-",round(C7.3[3],3)),paste0(round(C8.3[2],3),"-",round(C8.3[3],3)))
antro_b<-c(round(B1.1, 3),round(B2.1, 3),round(B3.1, 3),round(B4.1, 3),
           round(B5.1, 3),round(B6.1, 3),round(B7.1, 3),round(B8.1, 3))
antro_s_b<-c(round(B1.2, 3),round(B2.2, 3),round(B3.2, 3),round(B4.2, 3),
           round(B5.2, 3),round(B6.2, 3),round(B7.2, 3),round(B8.2, 3))
pheno_b<-c(round(B1.3, 3),round(B2.3, 3),round(B3.3, 3),round(B4.3, 3),
           round(B5.3, 3),round(B6.3, 3),round(B7.3, 3),round(B8.3, 3))
causes<-c("Cardiovascular", "Diabetes Mellitus", "Stroke", "Cancer", "Influenza/Pneumonia", "Nephritis/Nephrosis", "Alzheimer", "COPD")

tab3<-data.frame("Cause-specific mortality"=c(causes),
                 "AnthropoAge \n c-statistic (95%CI)"=c(paste0(antro," (",antro_ci,")")),
                 "S-AnthropoAge \n c-statistic (95%CI)"=c(paste0(antro_s," (",antros_ci,")")),
                 "PhenoAge \n c-statistic (95%CI)"=c(paste0(pheno," (",pheno_ci,")")),
                 "AnthropoAge vs. \n S-AnthropoAge"=c(antro_b),
                 "AnthropoAge vs. \n PhenoAge"=c(antro_s_b),
                 "S-AnthropoAge vs. \n PhenoAge"=c(pheno_b)
                 )

tab3<-`names<-`(tab3,c("Cause-specific mortality",
                       "AnthropoAge \n c-statistic (95%CI)",
                       "S-AnthropoAge \n c-statistic (95%CI)",
                       "PhenoAge \n c-statistic (95%CI)",
                       "AnthropoAge vs. \n S-AnthropoAge",
                       "AnthropoAge vs. \n PhenoAge",
                       "S-AnthropoAge vs. \n PhenoAge"))
tab3<-align(flextable(tab3),align = "center",part = "all") %>% 
  autofit()

doc <- read_docx() %>%
  body_add_flextable(value = tab3, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print( target = "tabla3.docx" )


#### Multidimensional aging ####
nhanes0$accel1<-ifelse(nhanes0$AnthropoAgeAccel>0, 1, 0)
nhanes0$accel2<-ifelse(nhanes0$AnthropoAgeAccel2>0, 1, 0)
nhanes0$accel3<-ifelse(nhanes0$PhenoAgeAccel>0, 1, 0)

### Kaplan-Meier
colors_border<-c("#6699CC","#994455")
colors_in<- rgb(col2rgb(colors_border)[1,],col2rgb(colors_border)[2,],
                col2rgb(colors_border)[3,],max=255,alpha=(50)*(255/100))

mod2_kma<-survfit(Surv(permth_int, mortstat) ~ accel1, data = nhanes0)
f3a<-ggsurvplot(mod2_kma, data = nhanes0, size = 1,palette = colors_border,conf.int = T, risk.table = T,pval = TRUE, xlab="Time (Months)",
                ylab="Survival probability", title="AnthropoAgeAccel", pval.coord = c(0, 0.8), legend.labs = c("Phys-aging","Accel-aging"),
                ylim= c(0.7,1.0), xlim=c(0, 120), break.y.by= c(0.1), break.x.by= c(30), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), text = element_text(hjust=0.5, family = "sans", size=15),
                         legend.text = element_text(hjust=0.5, size=12)) ) )
f3a$table <- f3a$table + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5, color="black", margin = margin(r=20), face="bold", size="12", angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size="15")) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
fig3a<-ggarrange(f3a$plot, f3a$table, heights = c(2, 0.7), ncol = 1, nrow = 2)


mod2_kmb<-survfit(Surv(permth_int, mortstat) ~ accel2, data = nhanes0)
f3b<-ggsurvplot(mod2_kmb, data = nhanes0, size = 1,palette = colors_border,conf.int = T, risk.table = T,pval = TRUE, xlab="Time (Months)",
                ylab="Survival probability", title="S-AnthropoAgeAccel", pval.coord = c(0, 0.8), legend.labs = c("Phys-aging","Accel-aging"),
                ylim= c(0.7,1.0), xlim=c(0, 120), break.y.by= c(0.1), break.x.by= c(30), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), text = element_text(hjust=0.5, family = "sans", size=15),
                         legend.text = element_text(hjust=0.5, size=12)) ))
f3b$table <- f3b$table + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5, color="black", margin = margin(r=20), face="bold", size="12", angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size="15")) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
fig3b<-ggarrange(f3b$plot, f3b$table, heights = c(2, 0.7), ncol = 1, nrow = 2)


mod2_kmc<-survfit(Surv(permth_int, mortstat) ~ accel3, data = nhanes0)
f3c<-ggsurvplot(mod2_kmc, data = nhanes0, size = 1,palette = colors_border,conf.int = T, risk.table = T,pval = TRUE, xlab="Time (Months)",
                ylab="Survival probability", title="PhenoAgeAccel", pval.coord = c(0, 0.8), legend.labs = c("Phys-aging","Accel-aging"),
                ylim= c(0.7,1.0), xlim=c(0, 120), break.y.by= c(0.1), break.x.by= c(30), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), text = element_text(hjust=0.5, family = "sans", size=15),
                         legend.text = element_text(hjust=0.5, size=12)) ))
f3c$table <- f3c$table + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5, color="black", margin = margin(r=20), face="bold", size="12", angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size="15")) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
fig3c<-ggarrange(f3c$plot, f3c$table, heights = c(2, 0.7), ncol = 1, nrow = 2)


#### Mortality per accelerated combination###
m_fin1<-flexsurvreg(Surv(permth_int, mortstat)~AnthropoAgeAccel+PhenoAgeAccel+Sex+Age+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")

nhanes0$accel_comb<-factor( 2*(nhanes0$PhenoAgeAccel>0) + (nhanes0$AnthropoAgeAccel>0)); table(nhanes0$accel_comb)
mod2_kmd<-survfit(Surv(permth_int, mortstat) ~ factor(accel_comb), data = nhanes0)
#summary(mod2_kmd, time=120)

colors_border2<-c("#BDDFFF", "#EBEF6A", "#D49400", "#6C2E38")
f3d<-ggsurvplot(mod2_kmd, data = nhanes0, size = 1,palette = colors_border2,conf.int = T, risk.table = T,pval = TRUE, xlab="Time (Months)",
                ylab="Survival probability", title="Multidomain aging", pval.coord = c(0, 0.8),
                legend.labs = c("Phys-aging", "Accel-AnthropoAge", "Accel-PhenoAge", "Multi-accel"),
                ylim= c(0.7,1.0), xlim=c(0, 120), break.y.by= c(0.1), break.x.by= c(30), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), text = element_text(hjust=0.5, family = "sans", size=15),
                         legend.text = element_text(hjust=0.5, size=12))) )
f3d$table <- f3d$table + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5, color="black", margin = margin(r=20), face="bold", size=12, angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size=15)) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
f3d$plot <- f3d$plot + scale_color_manual(values=colors_border2, labels=c("Phys-aging", "Accel\nAnthropoAge", "Accel\nPhenoAge", "Multi-accel")) +
  scale_fill_manual(values=colors_border2, labels=c("Phys-aging", "Accel\nAnthropoAge", "Accel\nPhenoAge", "Multi-accel"))
fig3d<-ggarrange(f3d$plot, f3d$table, heights = c(2, 0.7), ncol = 1, nrow = 2)

fig5<-ggarrange(fig3a, fig3b, "", "", fig3c, fig3d, labels=c("a","b","","","c","d"), heights=c(10,1,10), ncol=2, nrow=3)
ggsave(file = "Figure5.jpg", fig5, bg = "transparent",
       width = 45*0.8, height = 30*0.8, units=c("cm"),
       dpi = 500, limitsize = FALSE)

m_fin<-flexsurvreg(Surv(permth_int, mortstat)~accel_comb+Sex+Age+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")
round(cbind("coef"=exp(coef(m_fin))[3:5], exp(confint(m_fin))[3:5,]),2)

#Accelerated AnthropoAge-S
nhanes0$accel_comb2<-factor(2*(nhanes0$PhenoAgeAccel>0)+(nhanes0$AnthropoAgeAccel2>0)); table(nhanes0$accel_comb2)
m_fin2<-flexsurvreg(Surv(permth_int, mortstat)~accel_comb2+Sex+Age+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")
round(cbind("coef"=exp(coef(m_fin2))[3:5], exp(confint(m_fin2))[3:5,]),2)



#### Comorbidity Asessment ####

#AnthropoAge
nhanes0<- nhanes0 %>% mutate(comorb_cat = num_comorb %>% 
                               cut(c(-Inf,0,1, Inf))) %>%
  mutate(comorb_cat=factor(comorb_cat, labels = c("0", "1", ">=2"))) %>%
  mutate(age_cat = Age %>% cut(c(-Inf,49,69, Inf))) %>%
  mutate(age_cat=factor(age_cat, labels=c("<50y", "50-70y", ">=70y")))

options(scipen = 0)
s5a<-nhanes0 %>% ggplot(aes(x=comorb_cat, y=AnthropoAgeAccel, fill=comorb_cat))+geom_boxplot()+
  facet_wrap(~age_cat)+theme_pubclean()+labs(fill="Number of comorbidities", x="")+ theme(legend.position="top")+
  ggpubr::stat_compare_means(size=3.5, label.x.npc = 0.35) + theme(text = element_text(hjust=0.5)) +
  scale_fill_manual(values=c("#EBEF6A", "#D49400", "#6C2E38"))

s5b <- nhanes0 %>% group_by(age_cat,accel_comb, comorb_cat) %>%
  summarise(Counts = n()) %>% mutate(freq = Counts / sum(Counts)) %>% as.data.frame %>% 
  mutate("accel_comb"=factor(accel_comb, labels = c("Phys-Aging", "Anthropo-Accel", "Pheno-Accel", "Both-Accel"))) %>% 
  ggplot(aes(y=freq, x=comorb_cat,fill=comorb_cat)) + geom_bar(stat="identity", position='dodge', color="black")+
  scale_fill_manual(values=c("#EBEF6A", "#D49400", "#6C2E38"))+
  scale_y_continuous(labels = scales::percent_format(), limits = c(0,1))+facet_wrap(~Cohort)+
  theme_pubclean()+ylab("Frequency (%)")+xlab("")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_text(aes(label=round(freq*100,1)), vjust=-0.5, color="black", size=3.5)+
  facet_wrap(~age_cat+accel_comb)+labs(fill="Number of comorbidities")

supp5<-ggarrange(s5a, s5b, labels=letters[1:2])

ggsave(file = "SuppFig5.jpg", supp5, bg = "transparent",
       width = 37,height = 20, units=c("cm"),
       dpi = 500, limitsize = FALSE)


#S-AnthropoAge
s6a<-nhanes0 %>% ggplot(aes(x=comorb_cat, y=AnthropoAgeAccel2, fill=comorb_cat))+geom_boxplot()+
  facet_wrap(~age_cat)+theme_pubclean()+labs(fill="Number of comorbidities", x="")+ theme(legend.position="top")+
  ggpubr::stat_compare_means(size=3.5, label.x.npc = 0.35) + theme(text = element_text(hjust=0.5)) +
  scale_fill_manual(values=c("#EBEF6A", "#D49400", "#6C2E38")) + ylab("S-AnthropoAgeAccel")

s6b <- nhanes0 %>% group_by(age_cat,accel_comb2, comorb_cat) %>%
  summarise(Counts = n()) %>% mutate(freq = Counts / sum(Counts)) %>% as.data.frame %>% 
  mutate("accel_comb2"=factor(accel_comb2, labels = c("Phys-Aging", "Anthropo-Accel", "Pheno-Accel", "Both-Accel"))) %>% 
  ggplot(aes(y=freq, x=comorb_cat,fill=comorb_cat)) + geom_bar(stat="identity", position='dodge', color="black")+
  scale_fill_manual(values=c("#EBEF6A", "#D49400", "#6C2E38"))+
  scale_y_continuous(labels = scales::percent_format(), limits = c(0,1))+facet_wrap(~Cohort)+
  theme_pubclean()+ylab("Frequency (%)")+xlab("")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_text(aes(label=round(freq*100,1)), vjust=-0.5, color="black", size=3.5)+
  facet_wrap(~age_cat+accel_comb2)+labs(fill="Number of comorbidities")

supp6<-ggarrange(s6a, s6b, labels=letters[1:2])
ggsave(file = "SuppFig6.jpg", supp6, bg = "transparent",
       width = 37,height = 20, units=c("cm"),
       dpi = 500, limitsize = FALSE)




####---- AnthropoAgeAccel>0 Spiderplots ----####
 
### Data management: Anthropometry and DXA ###
dxa<-NHANES1 %>% filter(!is.infinite(PhenoAge)&!is.na(PhenoAge)) %>% 
  dplyr::select(SEQN,Age,mortstat,permth_int,Sex,Ethnicity,SEQN, PhenoAge,
                BMI,ICE,Weight,Height,Arm_length,Triceps_skinfold,Subscapular_skinfold,Thigh_circumference,Arm_circumference,
                DXXHEFAT,DXDHELE,DXXLAFAT,DXDLALE,DXXLLFAT,DXDLLLE,DXXRAFAT,DXDRALE,DXXRLFAT,DXDRLLE,DXXLSBMD,
                DXXTRFAT,DXDTRLE,DXDTOFAT,DXDTOLE,DXDTOBMD) %>% filter(!duplicated(SEQN)) %>% drop_na() %>% 
  mutate("Ethnicity"=factor(Ethnicity, labels = c("White", "Black", "Mexican-American", "Hispanic", "Other")),
         "tr_weight"=log(Weight),
         "tr_armc"=log(Arm_circumference),
         "tr_ice"=log(ICE),
         "tr_imc"=log(BMI),
         "tr_tric"=sqrt(Triceps_skinfold),
         "tr_subs"=sqrt(Subscapular_skinfold),
         "tr_arm"=log(Arm_length),
         "tr_thigh"=log(Thigh_circumference),
         "arm_fmi"=(DXXLAFAT/1000+DXXRAFAT/1000)/((Height/100)^2), #fmi=fat mass index
         "arm_lmi"=(DXDLALE/1000+DXDRALE/1000)/((Height/100)^2), #lmi=lean mass index
         "leg_fmi"=(DXXLLFAT/1000+DXXRLFAT/1000)/((Height/100)^2),
         "leg_lmi"=(DXDLLLE/1000+DXDRLLE/1000)/((Height/100)^2),
         "afmi"=(DXXLAFAT/1000+DXXRAFAT/1000+DXXLLFAT/1000+DXXRLFAT/1000)/((Height/100)^2),
         "almi"=(DXDLALE/1000+DXDRALE/1000+DXDLLLE/1000+DXDRLLE/1000)/((Height/100)^2),
         "head_fmi"=(DXXHEFAT/1000)/((Height/100)^2),
         "head_lmi"=(DXDHELE/1000)/((Height/100)^2),
         "trunk_fmi"=(DXXTRFAT/1000)/((Height/100)^2),
         "trunk_lmi"=(DXDTRLE/1000)/((Height/100)^2),
         "total_fmi"=(DXDTOFAT/1000)/((Height/100)^2),
         "total_lmi"=(DXDTOLE/1000)/((Height/100)^2),
         "total_fap"=(DXDTOFAT/1000)/Weight,
         "total_lep"=(DXDTOLE/1000)/Weight,
         "app_ratio"=(DXXLAFAT+DXXRAFAT+DXXLLFAT+DXXRLFAT)/(DXDLALE+DXDRALE+DXDLLLE+DXDRLLE),
         "trunk_ratio"=DXXTRFAT/DXDTRLE,
         "total_ratio"=DXDTOFAT/DXDTOLE,
         "lumbar_bmd"=DXXLSBMD,
         "total_bmd"=DXDTOBMD); nrow(dxa)

#AnthropoAge
p1F<-predict(gomp1aF, newdata = dxa %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M<-predict(gomp1aM, newdata = dxa %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

dxa$pred[dxa$Sex=="Women"]<-as.numeric(1-p1F$.pred)
dxa$pred[dxa$Sex=="Men"]<-as.numeric(1-p1M$.pred)
dxa$AnthropoAge[dxa$Sex=="Women"]<-(log(-sW*log(1-dxa$pred[dxa$Sex=="Women"]))-b0W)/b1W
dxa$AnthropoAge[dxa$Sex=="Men"]<-(log(-sM*log(1-dxa$pred[dxa$Sex=="Men"]))-b0M)/b1M

#S-AnthropoAge
p1F1<-predict(gomp1aF1, newdata = dxa %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M1<-predict(gomp1aM1, newdata = dxa %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

dxa$pred[dxa$Sex=="Women"]<-as.numeric(1-p1F1$.pred)
dxa$pred[dxa$Sex=="Men"]<-as.numeric(1-p1M1$.pred)
dxa$AnthropoAge2[dxa$Sex=="Women"]<-(log(-sW1*log(1-dxa$pred[dxa$Sex=="Women"]))-b0W1)/b1W1
dxa$AnthropoAge2[dxa$Sex=="Men"]<-(log(-sM1*log(1-dxa$pred[dxa$Sex=="Men"]))-b0M1)/b1M1

#Accelerated 
m1<-lm(AnthropoAge~Age, data=dxa %>% filter(Sex=="Men"))
dxa$AnthropoAgeAccel[dxa$Sex=="Men"]<-m1$residuals
m1<-lm(AnthropoAge~Age, data=dxa %>% filter(Sex=="Women"))
dxa$AnthropoAgeAccel[dxa$Sex=="Women"]<-m1$residuals

m1<-lm(AnthropoAge2~Age, data=dxa %>% filter(Sex=="Men"))
dxa$AnthropoAgeAccel2[dxa$Sex=="Men"]<-m1$residuals
m1<-lm(AnthropoAge2~Age, data=dxa %>% filter(Sex=="Women"))
dxa$AnthropoAgeAccel2[dxa$Sex=="Women"]<-m1$residuals

m1<-lm(PhenoAge~Age, data=dxa)
dxa$PhenoAgeAccel<-m1$residuals

with(dxa, table("BMI"=BMI>30, "AnthropoAgeAccel"=AnthropoAgeAccel>0)) %>% print() %>% mcnemar.test()
with(dxa, table("S-A"=AnthropoAgeAccel2>0, "A"=AnthropoAgeAccel>0)) %>% print() %>% mcnemar.test()



### Data management: PhenoAge components ###
pheno<-NHANES1%>%filter(!is.infinite(PhenoAge)&!is.na(PhenoAge))%>%
  dplyr::select(SEQN,Age,PhenoAge,mortstat,permth_int,Sex,Ethnicity,BMI,Thigh_circumference,Arm_circumference,
                Triceps_skinfold,Subscapular_skinfold,Leg_length,Arm_length,Height,Weight,ICE,
                Glucose,Albumin,RDW,MCV,LymP,WBC,Cr,CRP,ALP) %>% filter(!duplicated(SEQN)) %>% drop_na() %>%
  mutate("Ethnicity"=factor(Ethnicity, labels = c("White", "Black", "Mexican-American", "Hispanic", "Other")),
       "tr_weight"=log(Weight),
       "tr_armc"=log(Arm_circumference),
       "tr_ice"=log(ICE),
       "tr_imc"=log(BMI),
       "tr_tric"=sqrt(Triceps_skinfold),
       "tr_subs"=sqrt(Subscapular_skinfold),
       "tr_arm"=log(Arm_length),
       "tr_thigh"=log(Thigh_circumference))

#AnthropoAge
p1F<-predict(gomp1aF, newdata = pheno %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M<-predict(gomp1aM, newdata = pheno %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

pheno$pred[pheno$Sex=="Women"]<-as.numeric(1-p1F$.pred)
pheno$pred[pheno$Sex=="Men"]<-as.numeric(1-p1M$.pred)
pheno$AnthropoAge[pheno$Sex=="Women"]<-(log(-sW*log(1-pheno$pred[pheno$Sex=="Women"]))-b0W)/b1W
pheno$AnthropoAge[pheno$Sex=="Men"]<-(log(-sM*log(1-pheno$pred[pheno$Sex=="Men"]))-b0M)/b1M

#S-AnthropoAge
p1F1<-predict(gomp1aF1, newdata = pheno %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M1<-predict(gomp1aM1, newdata = pheno %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))

pheno$pred[pheno$Sex=="Women"]<-as.numeric(1-p1F1$.pred)
pheno$pred[pheno$Sex=="Men"]<-as.numeric(1-p1M1$.pred)
pheno$AnthropoAge2[pheno$Sex=="Women"]<-(log(-sW1*log(1-pheno$pred[pheno$Sex=="Women"]))-b0W1)/b1W1
pheno$AnthropoAge2[pheno$Sex=="Men"]<-(log(-sM1*log(1-pheno$pred[pheno$Sex=="Men"]))-b0M1)/b1M1

#Accelerated
m1<-lm(AnthropoAge~Age, data=pheno %>% filter(Sex=="Men"))
pheno$AnthropoAgeAccel[pheno$Sex=="Men"]<-m1$residuals
m1<-lm(AnthropoAge~Age, data=pheno %>% filter(Sex=="Women"))
pheno$AnthropoAgeAccel[pheno$Sex=="Women"]<-m1$residuals

m1<-lm(AnthropoAge2~Age, data=pheno %>% filter(Sex=="Men"))
pheno$AnthropoAgeAccel2[pheno$Sex=="Men"]<-m1$residuals
m1<-lm(AnthropoAge2~Age, data=pheno %>% filter(Sex=="Women"))
pheno$AnthropoAgeAccel2[pheno$Sex=="Women"]<-m1$residuals



### Spiderplot: AnthropoAge components ###
#*p<0.05, **p<0.001, ***p<0.0001
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Age~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Women"), wilcox.test(BMI~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(ICE~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Weight~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Arm_length~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Triceps_skinfold~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Subscapular_skinfold~AnthropoAgeAccel>0)))$p.value #**
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Thigh_circumference~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Arm_circumference~AnthropoAgeAccel>0)))$p.value #**
(with(dxa%>%filter(Sex=="Women"), wilcox.test(AnthropoAge~AnthropoAgeAccel>0)))$p.value #***

data1 <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
             dplyr::select(Age,BMI,ICE,Weight,Arm_length,Triceps_skinfold,Subscapular_skinfold,
                           Thigh_circumference,Arm_circumference,AnthropoAge,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("CA","BMI***","WHtR***","Weight***","Arm\nlength***","Tricipital\nskinfold***",
                           "Subscapular\nskinfold**","Thigh\ncircumference","Arm\ncircumference**","AnthropoAge***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6a<-ggplotify::as.ggplot(~radarchart(data1, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#*p<0.05, **p<0.001, ***p<0.0001
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Age~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(BMI~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(ICE~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Weight~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Arm_length~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Triceps_skinfold~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Subscapular_skinfold~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Thigh_circumference~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Arm_circumference~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(AnthropoAge~AnthropoAgeAccel>0)))$p.value #***

data2 <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
             dplyr::select(Age,BMI,ICE,Weight,Arm_length,Triceps_skinfold,Subscapular_skinfold,
                           Thigh_circumference,Arm_circumference,AnthropoAge,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("CA","BMI***","WHtR","Weight***","Arm length***","Tricipital\nskinfold***",
                           "Subscapular\nskinfold***","Thigh\ncircumference***","Arm\ncircumference***","AnthropoAge***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6d<-ggplotify::as.ggplot(~radarchart(data2, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))



### Spiderplot: DXA ###
#*p<0.05, **p<0.001, ***p<0.0001
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_fmi~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_lmi~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_ratio~AnthropoAgeAccel>0)))$p.value #*
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_fap~AnthropoAgeAccel>0)))$p.value #*
(with(dxa%>%filter(Sex=="Women"), wilcox.test(afmi~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Women"), wilcox.test(almi~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(trunk_fmi~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(trunk_ratio~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(lumbar_bmd~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_bmd~AnthropoAgeAccel>0)))$p.value #NS

data3 <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
             dplyr::select(total_fmi,total_lmi,total_ratio,total_fap,afmi,almi,trunk_fmi,
                           trunk_ratio,lumbar_bmd,total_bmd,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("Fat mass\nindex***","Lean mass\nindex***","Total fat-lean\nratio*","Fat\npercentage*","AFMI",
                           "ALMI***","Trunk fat\nmass index***","Trunk\nratio***","Lumbar BMD","Total BMD")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6b<-ggplotify::as.ggplot(~radarchart(data3, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#*p<0.05, **p<0.001, ***p<0.0001
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_fmi~AnthropoAgeAccel>0)))$p.value #**
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_lmi~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_ratio~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_fap~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(afmi~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(almi~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(trunk_fmi~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(trunk_ratio~AnthropoAgeAccel>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(lumbar_bmd~AnthropoAgeAccel>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_bmd~AnthropoAgeAccel>0)))$p.value #***

data4 <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
             dplyr::select(total_fmi,total_lmi,total_ratio,total_fap,afmi,almi,trunk_fmi,
                           trunk_ratio,lumbar_bmd,total_bmd,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("Fat mass\nindex**","Lean mass\nindex***","Total fat-lean\nratio","Fat\npercentage","AFMI***",
                           "ALMI***","Trunk fat\nmass index","Trunk\nratio","Lumbar BMD***","Total BMD***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6e<-ggplotify::as.ggplot(~radarchart(data4, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))



### Spiderplots: PhenoAge components ###
#*p<0.05, **p<0.001, ***p<0.0001
(with(pheno%>%filter(Sex=="Women"), wilcox.test(PhenoAge~AnthropoAgeAccel>0)))$p.value #*
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Glucose~AnthropoAgeAccel>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Cr~AnthropoAgeAccel>0)))$p.value #*
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Albumin~AnthropoAgeAccel>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(ALP~AnthropoAgeAccel>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(WBC~AnthropoAgeAccel>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(LymP~AnthropoAgeAccel>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(MCV~AnthropoAgeAccel>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Women"), wilcox.test(RDW~AnthropoAgeAccel>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(CRP~AnthropoAgeAccel>0)))$p.value #***

data5 <- ((accel1 <- pheno %>% filter(Sex=="Women") %>% 
             dplyr::select(PhenoAge,Glucose,Cr,Albumin,ALP,WBC,LymP,MCV,RDW,CRP,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("PhenoAge*","Glucose","Creatinine*","Albumin***","ALP***","WBC***",
                           "Lymphocyte\nPercentage***","MCV","RDW***","CRP***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6c<-ggplotify::as.ggplot(~radarchart(data5, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#*p<0.05, **p<0.001, ***p<0.0001
(with(pheno%>%filter(Sex=="Men"), wilcox.test(PhenoAge~AnthropoAgeAccel>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Glucose~AnthropoAgeAccel>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Cr~AnthropoAgeAccel>0)))$p.value #***
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Albumin~AnthropoAgeAccel>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(ALP~AnthropoAgeAccel>0)))$p.value #***
(with(pheno%>%filter(Sex=="Men"), wilcox.test(WBC~AnthropoAgeAccel>0)))$p.value #***
(with(pheno%>%filter(Sex=="Men"), wilcox.test(LymP~AnthropoAgeAccel>0)))$p.value #*
(with(pheno%>%filter(Sex=="Men"), wilcox.test(MCV~AnthropoAgeAccel>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(RDW~AnthropoAgeAccel>0)))$p.value #**
(with(pheno%>%filter(Sex=="Men"), wilcox.test(CRP~AnthropoAgeAccel>0)))$p.value #**

data6 <- ((accel1 <- pheno %>% filter(Sex=="Men") %>% 
             dplyr::select(PhenoAge,Glucose,Cr,Albumin,ALP,WBC,LymP,MCV,RDW,CRP,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("PhenoAge","Glucose","Creatinine***","Albumin","ALP***","WBC***",
                           "Lymphocyte\nPercentage*","MCV","RDW**","CRP**")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6f<-ggplotify::as.ggplot(~radarchart(data6, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


### Joint ###
dxa$accel_comb<-factor((dxa$PhenoAgeAccel>0)+2*(dxa$AnthropoAgeAccel>0))
g1<-nhanes0 %>% mutate(accel1=factor(accel1, labels = c("Physiological Aging", "Accelerated AnthropoAge"))) %>%
  ggplot(aes(x=accel1, y=Age, color=accel1))+scale_color_manual(values=colors_in)+geom_point(size=10)+
  theme(legend.position = "top",legend.title=element_text(size=20*.90), 
        legend.text=element_text(size=22),legend.key.size = unit(2,"cm"))+labs(color="Group")
legend<-get_legend(g1)


fig6<-ggarrange(ggarrange(annotate_figure(ggarrange(f6a,labels = "a", font.label = list(size=25)), left=text_grob("Women", face = "bold", size=25, rot = 90),top=text_grob("Anthropometry", face = "bold", size=25)),
                           annotate_figure(ggarrange(f6b,labels = "b", font.label = list(size=25)), top=text_grob("DXA", face = "bold", size=25)),
                           annotate_figure(ggarrange(f6c,labels = "c", font.label = list(size=25)), top=text_grob("PhenoAge", face = "bold", size=25)), ncol=3, nrow=1),
                 annotate_figure(ggarrange(f6d,f6e,f6f,labels = letters[4:6], font.label = list(size=25), nrow=1, ncol=3), left=text_grob("Men", face = "bold", size=25, rot = 90)), nrow=2, ncol=1)+
  geom_subview(x=0.5, y=1-0.985, subview=legend)

ggsave(file = "Figure4.jpg", 
       fig6,
       bg = "transparent",
       width = 62.5*1.15, 
       height = 45*1.05,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)



#### S-AnthropoAgeAccel>0 Spiderplots ####
### Spiderplot: AnthropoAge components ###
#*p<0.05, **p<0.001, ***p<0.0001
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Age~AnthropoAgeAccel2>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Women"), wilcox.test(BMI~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(ICE~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Weight~AnthropoAgeAccel2>0)))$p.value #**
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Arm_length~AnthropoAgeAccel2>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Triceps_skinfold~AnthropoAgeAccel2>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Subscapular_skinfold~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Thigh_circumference~AnthropoAgeAccel2>0)))$p.value #*
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Arm_circumference~AnthropoAgeAccel2>0)))$p.value #**
(with(dxa%>%filter(Sex=="Women"), wilcox.test(AnthropoAge~AnthropoAgeAccel2>0)))$p.value #***

data1 <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
             dplyr::select(Age,BMI,ICE,Weight,Arm_length,Triceps_skinfold,Subscapular_skinfold,
                           Thigh_circumference,Arm_circumference,AnthropoAge,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("CA","BMI***","WHtR***","Weight**","Arm\nlength","Tricipital\nskinfold",
                           "Subscapular\nskinfold***","Thigh\ncircumference*","Arm\ncircumference**","AnthropoAge***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4a<-ggplotify::as.ggplot(~radarchart(data1, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#*p<0.05, **p<0.001, ***p<0.0001
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Age~AnthropoAgeAccel2>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(BMI~AnthropoAgeAccel2>0)))$p.value #**
(with(dxa%>%filter(Sex=="Men"), wilcox.test(ICE~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Weight~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Arm_length~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Triceps_skinfold~AnthropoAgeAccel2>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Subscapular_skinfold~AnthropoAgeAccel2>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Thigh_circumference~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Arm_circumference~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(AnthropoAge~AnthropoAgeAccel2>0)))$p.value #*

data2 <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
             dplyr::select(Age,BMI,ICE,Weight,Arm_length,Triceps_skinfold,Subscapular_skinfold,
                           Thigh_circumference,Arm_circumference,AnthropoAge,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("CA","BMI**","WHtR***","Weight***","Arm length***","Tricipital\nskinfold",
                           "Subscapular\nskinfold","Thigh\ncircumference***","Arm\ncircumference***","AnthropoAge*")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4d<-ggplotify::as.ggplot(~radarchart(data2, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))



### Spiderplot: DXA ###
#*p<0.05, **p<0.001, ***p<0.0001
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_fmi~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_lmi~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_ratio~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_fap~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(afmi~AnthropoAgeAccel2>0)))$p.value #*
(with(dxa%>%filter(Sex=="Women"), wilcox.test(almi~AnthropoAgeAccel2>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Women"), wilcox.test(trunk_fmi~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(trunk_ratio~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(lumbar_bmd~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_bmd~AnthropoAgeAccel2>0)))$p.value #***

data3 <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
             dplyr::select(total_fmi,total_lmi,total_ratio,total_fap,afmi,almi,trunk_fmi,
                           trunk_ratio,lumbar_bmd,total_bmd,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("Fat mass\nindex***","Lean mass\nindex***","Total fat-lean\nratio***","Fat\npercentage***","AFMI*",
                           "ALMI","Trunk fat\nmass index***","Trunk\nratio***","Lumbar BMD***","Total BMD***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4b<-ggplotify::as.ggplot(~radarchart(data3, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#*p<0.05, **p<0.001, ***p<0.0001
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_fmi~AnthropoAgeAccel2>0)))$p.value #*
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_lmi~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_ratio~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_fap~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(afmi~AnthropoAgeAccel2>0)))$p.value #NS
(with(dxa%>%filter(Sex=="Men"), wilcox.test(almi~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(trunk_fmi~AnthropoAgeAccel2>0)))$p.value #**
(with(dxa%>%filter(Sex=="Men"), wilcox.test(trunk_ratio~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(lumbar_bmd~AnthropoAgeAccel2>0)))$p.value #***
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_bmd~AnthropoAgeAccel2>0)))$p.value #***

data4 <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
             dplyr::select(total_fmi,total_lmi,total_ratio,total_fap,afmi,almi,trunk_fmi,
                           trunk_ratio,lumbar_bmd,total_bmd,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("Fat mass\nindex*","Lean mass\nindex***","Total fat-lean\nratio***","Fat\npercentage***","AFMI",
                           "ALMI***","Trunk fat\nmass index**","Trunk\nratio***","Lumbar BMD***","Total BMD***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4e<-ggplotify::as.ggplot(~radarchart(data4, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))



### Spiderplots: PhenoAge components ###
#*p<0.05, **p<0.001, ***p<0.0001
(with(pheno%>%filter(Sex=="Women"), wilcox.test(PhenoAge~AnthropoAgeAccel2>0)))$p.value #**
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Glucose~AnthropoAgeAccel2>0)))$p.value #*
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Cr~AnthropoAgeAccel2>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Albumin~AnthropoAgeAccel2>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(ALP~AnthropoAgeAccel2>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(WBC~AnthropoAgeAccel2>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(LymP~AnthropoAgeAccel2>0)))$p.value #***
(with(pheno%>%filter(Sex=="Women"), wilcox.test(MCV~AnthropoAgeAccel2>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Women"), wilcox.test(RDW~AnthropoAgeAccel2>0)))$p.value #**
(with(pheno%>%filter(Sex=="Women"), wilcox.test(CRP~AnthropoAgeAccel2>0)))$p.value #***

data5 <- ((accel1 <- pheno %>% filter(Sex=="Women") %>% 
             dplyr::select(PhenoAge,Glucose,Cr,Albumin,ALP,WBC,LymP,MCV,RDW,CRP,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("PhenoAge**","Glucose*","Creatinine","Albumin***","ALP***","WBC***",
                           "Lymphocyte\nPercentage***","MCV","RDW**","CRP***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4c<-ggplotify::as.ggplot(~radarchart(data5, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#*p<0.05, **p<0.001, ***p<0.0001
(with(pheno%>%filter(Sex=="Men"), wilcox.test(PhenoAge~AnthropoAgeAccel2>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Glucose~AnthropoAgeAccel2>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Cr~AnthropoAgeAccel2>0)))$p.value #***
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Albumin~AnthropoAgeAccel2>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(ALP~AnthropoAgeAccel2>0)))$p.value #***
(with(pheno%>%filter(Sex=="Men"), wilcox.test(WBC~AnthropoAgeAccel2>0)))$p.value #***
(with(pheno%>%filter(Sex=="Men"), wilcox.test(LymP~AnthropoAgeAccel2>0)))$p.value #*
(with(pheno%>%filter(Sex=="Men"), wilcox.test(MCV~AnthropoAgeAccel2>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(RDW~AnthropoAgeAccel2>0)))$p.value #NS
(with(pheno%>%filter(Sex=="Men"), wilcox.test(CRP~AnthropoAgeAccel2>0)))$p.value #***

data6 <- ((accel1 <- pheno %>% filter(Sex=="Men") %>% 
             dplyr::select(PhenoAge,Glucose,Cr,Albumin,ALP,WBC,LymP,MCV,RDW,CRP,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
            `colnames<-`(c("PhenoAge","Glucose","Creatinine***","Albumin","ALP***","WBC***",
                           "Lymphocyte\nPercentage*","MCV","RDW","CRP***")) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4f<-ggplotify::as.ggplot(~radarchart(data6, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


### Joint ###
supp4<-ggarrange(ggarrange(annotate_figure(ggarrange(s4a,labels = "a", font.label = list(size=25)), left=text_grob("Women", face = "bold", size=25, rot = 90),top=text_grob("Anthropometry", face = "bold", size=25)),
                 annotate_figure(ggarrange(s4b,labels = "b", font.label = list(size=25)), top=text_grob("DXA", face = "bold", size=25)),
                 annotate_figure(ggarrange(s4c,labels = "c", font.label = list(size=25)), top=text_grob("PhenoAge", face = "bold", size=25)), ncol=3, nrow=1),
                 annotate_figure(ggarrange(s4d,s4e,s4f,labels = letters[4:6], font.label = list(size=25), nrow=1, ncol=3), left=text_grob("Men", face = "bold", size=25, rot = 90)), nrow=2, ncol=1)+
                   geom_subview(x=0.5, y=1-0.985, subview=legend)
                 
ggsave(file = "SuppFig4.jpg", 
       supp4,
       bg = "transparent",
       width = 62.5*1.15, 
       height = 45*1.05,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)
