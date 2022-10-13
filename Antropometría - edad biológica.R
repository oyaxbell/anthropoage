# AnthropoAge, a novel approach to integrate body composition into the estimation of biological age 
# Data Analysis: Carlos A. Fermin-Martinez, Alejandro Marquez-Salinas & Omar Yaxmehen Bello-Chavolla
# Latest version of Analysis October 10th, 2022
# For any question regarding analysis contact Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

####---- Database management ----####
pacman::p_load(haven, tidyverse, ggpubr, lmtest, nortest, gtools, data.table, caret, glmnet, survival, flextable, blandr, BlandAltmanLeh, corrplot,
               rms, bestNormalize, flexsurv, pROC, timeROC, fmsb, factoextra, gridExtra,  nhanesA, wesanderson,forestmodel, ggedit,dummy,
               FactoMineR, fpc, NbClust, ggimage, glmnet, ggsci, survminer, cluster, ggplotify, UpSetR, nortest, viridis, officer, magrittr)

#Extra functions 
conc95 <- function(x){
  y <- (summary(x)$concordance[2])*1.96
  c <- summary(x)$concordance[1]
  c_low <- c-y; c_up <- c+y
  `names<-`(c(c,c_low,c_up),c("Concordance","Lower 95-CI","Upper 95-CI"))}

hr95 <- function(x){
  a <- summary(x)$conf.int[1,1]
  a_low <- summary(x)$conf.int[1,3]
  a_up <- summary(x)$conf.int[1,4]
  `names<-`(c(a,a_low,a_up),c("HR","Lower 95-CI","Upper 95-CI"))}

ci_quick<-function(x){
  a=round(x,3)
  paste0("HR=",a[1]," (",a[2],"-",a[3],")")}

p_aster <- function(x){
  y=x<c(0.05,0.01,0.001)
  z=sum(y)
  if(z==3){"***"}
  else if(z==2){"**"}
  else if(z==1){"*"}
  else if(z==0){""}}


#Please save as WINDOWS-1252
#setwd("C:/Users/investigacion/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/AnthropoAge")
#setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/AnthropoAge")
#setwd("~/UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ALEJANDRO MARQUEZ SALINAS - Antropometría")
#setwd("C:/Users/facmed/UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ALEJANDRO MARQUEZ SALINAS - Antropometría")
#setwd("C:/Users/facmed/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/Antropometría")

setwd("/Users/carlosfermin/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/Antropometría")
options(warn = 0)

### Loading databases ###
#NHANES-III
NHANES3<-fread("nhanes3.csv", na = c("", "N/A", "NA","na", "#N/A", "88888", "8888", "888888", "888","5555","999","9999", "99999", "999999", "9998"))
mortalidad3<-fread("nhanes3_mort.csv"); mortalidad3$SEQN<-mortalidad3$seqn
#NHANES-IV
NHANES<-fread("nhanes4_dxa.csv",na = c("", "N/A", "NA","na", "#N/A", "88888", "8888", "888888", "888","5555","999","9999", "99999", "999999"))
mortalidad<-fread("nhanes_mortalidad.csv"); mortalidad$SEQN<-mortalidad$seqn
source("functions_predict.R")

### NHANES-III management ###
NHANES0 <- NHANES3 %>% filter(HSAGEIR.x>=20)%>% 
  mutate("BMI"=BMPBMI, "Thigh_circumference"=BMPTHICI, "Arm_circumference"=BMPARMC, "Weight"= BMPWT,
         "Height"=BMPHT, "Waist"=BMPWAIST, "Triceps_skinfold"=BMPTRI, "Subscapular_skinfold"=BMPSUB,
         "Leg_length"=BMPLEG, "Arm_length"=BMPARML, "ICE"=BMPWAIST/BMPHT, "Diabetes"=HAD1,
         "Hypertension"=HAE2, "Asthma"=HAC1E, "Arthritis"=HAC1A, "Heart_failure"=HAC1C,
         "Heart_attack"=HAF10, "Stroke"=HAC1D, "Emphysema"=HAC1G, "Bronchitis"=HAC1F, "Malignancy"=HAC1O,
         "Sex"=factor(HSSEX.x, levels= c(1,2),labels= c("Men", "Women")), "Age"=HSAGEIR.x, "Ethnicity"=DMARETHN.x)
#Race/Ethnicity: 1=White, 2=Black, 3=Mexican-American, 4=Other

#PhenoAge
NHANES0$PhenoAge = 141.5 + ((log(-0.00553*log(1-(1-exp((-1.51714*exp(-19.907-0.0336*NHANES0$AMP+0.0095*NHANES0$CEPSI + 0.1953*NHANES0$G1PSI+
                                                                       0.0954*log(NHANES0$CRP)-0.0120*NHANES0$LMPPCNT+0.0268*NHANES0$MVPSI+
                                                                       0.3306*NHANES0$RWP+0.00188*NHANES0$AP+0.0554*NHANES0$WCP+
                                                                       0.0804*NHANES0$HSAGEIR.x))/(0.0076927)))))))/(0.09165)
#Cause-specific mortality
NHANES0 <- merge(NHANES0,mortalidad3,by="SEQN")
NHANES0$ucod_leading[is.na(NHANES0$ucod_leading)] <- 11
d1<-dummy::dummy(NHANES0 %>% transmute(as.factor(ucod_leading)))
colnames(d1)<-c("Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                "Nephone_diseases","Other_causes","Alive"); NHANES0<-cbind(NHANES0, d1)

nhanes0<-NHANES0 %>% dplyr::select(SEQN, Age, Sex, Ethnicity, mortstat,permth_int,PhenoAge,BMI, Thigh_circumference, Arm_circumference,
                                   Triceps_skinfold,Subscapular_skinfold,Leg_length, Arm_length, Waist, Height, Weight, ICE,
                                   "Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                                   "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                                   "Nephone_diseases","Other_causes","Alive", Diabetes, Hypertension, Asthma, Arthritis, Heart_failure, 
                                   Heart_attack, Stroke, Emphysema, Bronchitis, Malignancy); nhanes0$id <- rep(1,nrow(nhanes0))

### NHANES-IV management ### 
NHANES1 <- NHANES %>% filter(RIDAGEYR>=20) %>% 
  mutate("DXDTOFAT_N" = (DXDTOFAT/1000)/((BMXWT/100)^2), "DXXTRFAT_N" = (DXXTRFAT/1000)/((BMXWT/100)^2), "DXXHEFAT_N" = (DXXHEFAT/1000)/((BMXWT/100)^2),
         "DXXLAFAT_N" = (DXXLAFAT/1000)/((BMXWT/100)^2), "DXXRAFAT_N" = (DXXRAFAT/1000)/((BMXWT/100)^2), "DXXLLFAT_N" = (DXXLLFAT/1000)/((BMXWT/100)^2),
         "DXXRLFAT_N" = (DXXRLFAT/1000)/((BMXWT/100)^2), "DXXTRFAT_DXDTOFAT" = (DXXTRFAT)/(DXDTOFAT), "Weight" = BMXWT, "Waist" = BMXWAIST, "Height" = BMXHT,
         "BMI" = BMXBMI, "Calf_circumference" = BMXCALF, "Arm_circumference" = BMXARMC, "Thigh_circumference" = BMXTHICR, "Triceps_skinfold" = BMXTRI,
         "Subscapular_skinfold" = BMXSUB, "Leg_length" = BMXLEG, "ICE" = BMXWAIST/BMXHT, "Arm_length"= BMXARML, "METSIR" = (log(Glucose*2+TG)*BMXBMI)/log(HDL),
         "EXTSUP_FAT" = DXXLAFAT_N + DXXRAFAT_N, "EXTINF_FAT" = DXXLLFAT_N + DXXRLFAT_N, "DXDHELE_N" = (DXDHELE/1000), "DXDLALE_M" = (DXDLALE/1000), "DXDRALE_N" = (DXDRALE/1000),
         "DXDRLLE_N" = (DXDRLLE/1000), "DXDLALE_N" = (DXDLALE/1000), "DXDLLLE_N" = (DXDLLLE/1000), "DXDTOLE_N" = (DXDTOLE/1000), "DXDTRLE_N" = (DXDTRLE/1000)) %>%    
  mutate("METS_VF" = 4.466 + 0.011*(log(METSIR)^3)+ 3.239*(log(ICE)^3)-0.319*(2-RIAGENDR) + 0.594*(log(RIDAGEYR)), "EXTSUP_FAT_N" = (EXTSUP_FAT/1000)/((Height/100)^2),
         "EXTINF_FAT_N" = (EXTINF_FAT/1000)/((Height/100)^2), "EXTSUP_FAT_DXDTOFAT" = (EXTSUP_FAT/1000)/(DXDTOFAT/1000), "EXTINF_FAT_DXDTOFAT" = (EXTINF_FAT/1000)/(DXDTOFAT/1000)) %>% 
  mutate("Sex"=factor(RIAGENDR, levels= c(1,2),labels= c("Men", "Women")), "Age"=RIDAGEYR, "Diabetes"=DIQ010, "Hypertension"=BPQ020, "Asthma"=MCQ010,
         "Arthritis"=MCQ160A, "Heart_failure"=MCQ160B, "Heart_attack"=MCQ160E, "Stroke"=MCQ160F, "Emphysema"=MCQ160G, "Bronchitis"=MCQ160K, "Malignancy"=MCQ220)

#Recode race/ethnicity
NHANES1$Ethnicity[NHANES1$RIDRETH1==5] <- 4
NHANES1$Ethnicity[NHANES1$RIDRETH1==1] <- 3; NHANES1$Ethnicity[NHANES1$RIDRETH1==2] <- 4
NHANES1$Ethnicity[NHANES1$RIDRETH1==3] <- 1; NHANES1$Ethnicity[NHANES1$RIDRETH1==4] <- 2

#PhenoAge
NHANES1$PhenoAge <- 141.5 + (((log(-0.00553* log(1-(1-exp((-1.51714*exp(-19.907 - 0.0336*NHANES1$Albumin + 0.0095*(NHANES1$Creatinine*88.4) + 
                                                                          0.1953*(NHANES1$Glucose*0.0555) + 0.0954*log(NHANES1$CRP) -
                                                                          0.0120*NHANES1$LymphP + 0.0268*NHANES1$MCV + 0.3306*NHANES1$RDW + 
                                                                          0.00188*NHANES1$ALP + 0.0554*NHANES1$WBC +
                                                                          0.0804*NHANES1$Age))/(0.0076927)))))))/(0.09165))
#Cause-specific mortality
NHANES1 <- merge(NHANES1,mortalidad,by="SEQN")
NHANES1$ucod_leading[is.na(NHANES1$ucod_leading)] <- 11
d1<-dummy::dummy(NHANES1 %>% transmute(as.factor(ucod_leading)))
colnames(d1)<-c("Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                "Nephone_diseases","Other_causes","Alive"); NHANES1<-cbind(NHANES1,d1)

nhanes1<-NHANES1 %>% dplyr::select(SEQN, Age, Sex, Ethnicity, mortstat,permth_int,PhenoAge,BMI, Thigh_circumference, Arm_circumference,
                                   Triceps_skinfold,Subscapular_skinfold,Leg_length, Arm_length, Waist, Height, Weight, ICE,
                                   "Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                                   "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                                   "Nephone_diseases","Other_causes","Alive", Diabetes, Hypertension, Asthma, Arthritis, Heart_failure, 
                                   Heart_attack, Stroke, Emphysema, Bronchitis, Malignancy); nhanes1$id<-rep(2, nrow(nhanes1))


### NHANES III and IV ###
nhanes0 <- nhanes0 %>% filter(!duplicated(SEQN)); nhanes1 <- nhanes1 %>% filter(!duplicated(SEQN))
nhanes <- rbind(nhanes0, nhanes1)

#Comorbidities
nhanes$Diabetes<-ifelse(nhanes$Diabetes==1, 1, 0); nhanes$Asthma<-ifelse(nhanes$Asthma==1, 1, 0)
nhanes$Arthritis<-ifelse(nhanes$Arthritis==1, 1, 0); nhanes$Heart_failure<-ifelse(nhanes$Heart_failure==1, 1, 0)
nhanes$Heart_attack<-ifelse(nhanes$Heart_attack==1, 1, 0); nhanes$Emphysema<-ifelse(nhanes$Emphysema==1, 1, 0)
nhanes$Bronchitis<-ifelse(nhanes$Bronchitis==1, 1, 0); nhanes$Malignancy<-ifelse(nhanes$Malignancy==1, 1, 0)
nhanes$Stroke<-ifelse(nhanes$Stroke==1, 1, 0); nhanes$Hypertension<-ifelse(nhanes$Hypertension==1, 1, 0)

#Number of comorbidities
nhanes$num_comorb<-nhanes$Diabetes+nhanes$Asthma+nhanes$Arthritis+nhanes$Heart_failure+nhanes$Heart_attack+nhanes$Emphysema+nhanes$Bronchitis+nhanes$Malignancy+nhanes$Stroke+nhanes$Hypertension

#Race/Ethnicity
nhanes$Ethnicity<-factor(nhanes$Ethnicity, labels = c("White", "Black", "Mexican-American", "Other"))

####---- Missing values management ----####
#Step 1: Remove duplicated subjects
nrow(nhanes); table(nhanes$id)
nhanes_0<-nhanes %>% filter(!duplicated(SEQN))
nrow(nhanes_0); table(nhanes_0$id) #n=23651; 13866 from NHANES-III, 9785 from NHANES-IV

#Step 2: Remove missing values from ALL variables except PhenoAge
nhanes0<-nhanes_0[,-c(7)] %>% drop_na()
nhanes0<-(merge(nhanes0,nhanes_0[,c(1,7)],by = "SEQN"))

#Step 3: Replace all infinite values from PhenoAge with NA
nhanes0$PhenoAge[is.infinite(nhanes0$PhenoAge)]<-NA 
nrow(nhanes0); table(nhanes0$id) #n=18794; 11774 from NHANES-III, 7020 from NHANES-IV

#IF we removed PhenoAge missing values
nrow(nhanes0 %>% na.omit); table((nhanes0 %>% na.omit)$id) #n=17450; 10823 from NHANES-III, 6627 from NHANES-IV

#Step 4: Divide database according to sex
nhanes_men<-nhanes0 %>% filter(Sex=="Men"); nhanes_women<-nhanes0 %>% filter(Sex=="Women")
nrow(nhanes_men); table(nhanes_men$id) #9289 men; 5728 from NHANES-III, 3561 from NHANES-IV
nrow(nhanes_women); table(nhanes_women$id) #9505 women; 6046 from NHANES-III, 3459 from NHANES-IV


####---- Variable transformations ----####
#Transformed anthropometric measures
nhanes0$tr_weight<-log(nhanes0$Weight)
nhanes0$tr_imc<-log(nhanes0$BMI)
nhanes0$tr_ice<-nhanes0$ICE**(1/3)
nhanes0$tr_subs<-(nhanes0$Subscapular_skinfold)**(1/3)
nhanes0$tr_tric<-(nhanes0$Triceps_skinfold)**(1/3)
nhanes0$tr_height<-(nhanes0$Height)**(1/3)
nhanes0$tr_leg<-nhanes0$Leg_length #No transformation
nhanes0$tr_arm<-log(nhanes0$Arm_length)
nhanes0$tr_thigh<-log(nhanes0$Thigh_circumference)
nhanes0$tr_armc<-sqrt(nhanes0$Arm_circumference)

#Supplementary figure 1
ttab_s <- nhanes0 %>% select(
  Weight, BMI, ICE, Subscapular_skinfold, Triceps_skinfold, Height, Leg_length, Arm_length, Thigh_circumference, Arm_circumference) %>%
  apply(2,ad.test) %>% sapply(extract, "statistic") %>% as.numeric %>% cbind(nhanes0 %>% select(
    tr_weight, tr_imc, tr_ice, tr_subs, tr_tric, tr_height, tr_leg, tr_arm, tr_thigh, tr_armc) %>% apply(2,ad.test) %>%
      sapply(extract, "statistic") %>% as.numeric()) %>% apply(2,round,3); options(scipen = -1000); ttab_pA <- nhanes0 %>% select(
        Weight, BMI, ICE, Subscapular_skinfold, Triceps_skinfold, Height, Leg_length, Arm_length, Thigh_circumference, Arm_circumference) %>%
  apply(2,ad.test) %>% sapply(extract, "p.value") %>% as.character() %>% strsplit(split="e-") %>%
  sapply(as.numeric) %>% round(3); options(scipen = 1000); ttab_pA <- paste0(ttab_pA[1,], "e-", ttab_pA[2,]); options(
    scipen = -1000); ttab_pB <- nhanes0 %>% select(
      tr_weight, tr_imc, tr_ice, tr_subs, tr_tric, tr_height, tr_leg, tr_arm, tr_thigh, tr_armc) %>%
  apply(2,ad.test) %>% sapply(extract, "p.value") %>% as.character() %>% strsplit(split="e-") %>%
  sapply(as.numeric) %>% round(3); options(scipen = 1000); ttab_pB <- paste0(ttab_pB[1,], "e-", ttab_pB[2,]); ttab_n <- c(
    "Weight (kg)", "Body Mass Index (kg/m2)", "Waist-to-Height ratio", "Subscapular skinfold (cm)", "Triceps skinfold (cm)",
    "Height (cm)", "Leg length (cm)", "Arm length (cm)", "Thigh circumference (cm)", "Arm circumference (cm)"); ttab_t <-c(
      "Logarithmic", "Logarithmic", "Cubic root", "Cubic root", "Cubic root", "Cubic root", "None", "Logarithmic", "Logarithmic", "Square root"
      ); tr_tab <- data.frame(ttab_n, round((ttab_s[,1]), 3), ttab_pA, ttab_t, round((ttab_s[,2]), 3), ttab_pB) %>% `colnames<-`(c(
        "Variable", "A-statistic (prior)", "p-value (prior)", "Transformation", "A-statistic (transformed)", "p-value (transformed)")
        ); supptab2<-align(flextable(tr_tab,cwidth = c(2,1.5,1.5,2,1)),align = "center",part = "all"); options(scipen = 10)

suppFZ <- ggarrange(nrow=2, ncol=5, common.legend = T, legend = "bottom", 
                    (ggplot(nhanes0, aes(x=tr_weight, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Weight (log)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_imc, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Body Mass Index (log)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_ice, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Waist-to-Height ratio (cubic root)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_subs, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Subscapular skinfold (cubic root)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_tric, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Triceps skinfold (cubic root)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_height, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Height (cubic root)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_leg, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Leg length")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_arm, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Arm length (log)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_thigh, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Thigh circumference (log)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top")),
                    (ggplot(nhanes0, aes(x=tr_armc, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+
                       theme_pubclean()+xlab("Arm circumference (square root)")+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "top"))) %>% 
  annotate_figure(top = text_grob("Distribution of transformed anthropometric measurements", face = "bold", size = 20)) %>% 
  ggarrange(" ", ncol=1, nrow=2, heights = c(1,0.5)) + annotation_custom(tableGrob(
    tr_tab, rows=NULL, theme=ttheme_minimal(base_size = 13, padding = unit(c(23, 3.5), "mm"))), xmin=0, xmax=1, ymin=0.12, ymax=0.22)

ggsave(suppFZ,filename = "SuppFig1.jpg", width = 38.1, height = 23.8, units=c("cm"), dpi = 300, limitsize = FALSE)


####---- Descriptive statistics ----####
## Flowchart diagram
#Overall NHANES population
N1 <- NHANES3 %>% filter(!duplicated(SEQN)) %>% nrow + NHANES %>% filter(!duplicated(SEQN)) %>% nrow; N1
N1.III <- NHANES3 %>% filter(!duplicated(SEQN)) %>% nrow; N1.III; N1.IV <- NHANES %>% filter(!duplicated(SEQN)) %>% nrow; N1.IV
#Filtered by age >=20 years
NHANES3 %>% filter(!duplicated(SEQN) & HSAGEIR.x>=20) %>% nrow + NHANES %>% filter(!duplicated(SEQN) & RIDAGEYR>=20) %>% nrow
NHANES3 %>% filter(!duplicated(SEQN) & HSAGEIR.x>=20) %>% nrow; NHANES %>% filter(!duplicated(SEQN) & RIDAGEYR>=20) %>% nrow
#After removing missing data from anthropometric and mortality data
nrow(nhanes0); table(nhanes0$id)
#After removing missing data from laboratory data and PhenoAge
nrow(nhanes0 %>% na.omit); table(nhanes0 %>% na.omit %>% select(id))

#Anthropometry: Men vs Women
((nhanes0[,c(7:17)]) %>% apply(2,function(x){wilcox.test(x~nhanes0$Sex)}) %>% sapply(extract,"p.value") %>%
    as.numeric %>% `names<-`(names(nhanes0[,c(7:17)])))[c(10,8,7,11,2,9,12,6:3)-1]


##Supplementary table 2
# Sex
Male0 <- table(nhanes0$Sex)[2]
pMale0 <- round(((table(nhanes0$Sex)%>%prop.table())[2])*100,1)
Male1 <- table((nhanes0%>% filter(id==1))$Sex)[2]
pMale1 <- round(((table((nhanes0%>% filter(id==1))$Sex)%>%prop.table())[2])*100,1)
Male2 <- table((nhanes0%>% filter(id==2))$Sex)[2]
pMale2 <- round(((table((nhanes0%>% filter(id==2))$Sex)%>%prop.table())[2])*100,1)
p1 <- format.pval(prop.test(x=c(Male1,Male2), n=c(nrow((nhanes0%>% filter(id==1))),nrow((nhanes0%>% filter(id==2)))))$p.value, eps = .001, digits = 3) 

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

#Ethnicity
#White
WEth0 <- table(nhanes0$Ethnicity)[1]
WpEth0 <- round(((table(nhanes0$Ethnicity)%>%prop.table())[1])*100,1)
WEth1 <- table((nhanes0%>% filter(id==1))$Ethnicity)[1]
WpEth1 <- round(((table((nhanes0%>% filter(id==1))$Ethnicity)%>%prop.table())[1])*100,1)
WEth2 <- table((nhanes0%>% filter(id==2))$Ethnicity)[1]
WpEth2 <- round(((table((nhanes0%>% filter(id==2))$Ethnicity)%>%prop.table())[1])*100,1)
#Black
BEth0 <- table(nhanes0$Ethnicity)[2]
BpEth0 <- round(((table(nhanes0$Ethnicity)%>%prop.table())[2])*100,1)
BEth1 <- table((nhanes0%>% filter(id==1))$Ethnicity)[2]
BpEth1 <- round(((table((nhanes0%>% filter(id==1))$Ethnicity)%>%prop.table())[2])*100,1)
BEth2 <- table((nhanes0%>% filter(id==2))$Ethnicity)[2]
BpEth2 <- round(((table((nhanes0%>% filter(id==2))$Ethnicity)%>%prop.table())[2])*100,1)
#Mexican-American
MEth0 <- table(nhanes0$Ethnicity)[3]
MpEth0 <- round(((table(nhanes0$Ethnicity)%>%prop.table())[3])*100,1)
MEth1 <- table((nhanes0%>% filter(id==1))$Ethnicity)[3]
MpEth1 <- round(((table((nhanes0%>% filter(id==1))$Ethnicity)%>%prop.table())[3])*100,1)
MEth2 <- table((nhanes0%>% filter(id==2))$Ethnicity)[3]
MpEth2 <- round(((table((nhanes0%>% filter(id==2))$Ethnicity)%>%prop.table())[3])*100,1)
pEthn <- ((with(nhanes0, table(Ethnicity,id)) %>% prop.test)$p.value) %>% 
  format.pval(eps = .001, digits = 3)

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
p3 <- format.pval(wilcox.test(nhanes0$PhenoAge~nhanes0$id)$p.value, eps = .001, digits = 3) 

# >= 1 comorb
nhanes0$comorb<-ifelse(nhanes0$num_comorb>=1, 1, 0)
Comorb0 <- table(nhanes0$comorb)[2]
pComorb0 <- round(((table(nhanes0$comorb)%>%prop.table())[2])*100,1)
Comorb1 <- table((nhanes0%>% filter(id==1))$comorb)[2]
pComorb1 <- round(((table((nhanes0%>% filter(id==1))$comorb)%>%prop.table())[2])*100,1)
Comorb2 <- table((nhanes0%>% filter(id==2))$comorb)[2]
pComorb2 <- round(((table((nhanes0%>% filter(id==2))$comorb)%>%prop.table())[2])*100,1)
p4 <- format.pval(prop.test(x=c(Comorb1,Comorb2), n=c(nrow((nhanes0%>% filter(id==1))),nrow((nhanes0%>% filter(id==2)))))$p.value, eps = .001, digits = 3) 

# Comorb list
nhanestab3 <- nhanes0 %>% filter(id==1); nhanestab4 <- nhanes0 %>% filter(id==2)
comorb_list0 <- paste0(
  (nhanes0[,c(29:38)] %>% apply(2,sum)) %>% sort(decreasing = T), " (",
  ((((nhanes0[,c(29:38)] %>% apply(2,sum)) / nrow(nhanes0)) %>% sort(decreasing = T) %>% round(3)) * 100), ")")
comorb_list1 <- paste0(
  (nhanestab3[,c(29:38)] %>% apply(2,sum)) %>% sort(decreasing = T), " (",
  ((((nhanestab3[,c(29:38)] %>% apply(2,sum)) / nrow(nhanestab3)) %>% sort(decreasing = T) %>% round(3)) * 100), ")")
comorb_list2 <- paste0(
  (nhanestab4[,c(29:38)] %>% apply(2,sum)) %>% sort(decreasing = T), " (",
  ((((nhanestab4[,c(29:38)] %>% apply(2,sum)) / nrow(nhanestab4)) %>% sort(decreasing = T) %>% round(3)) * 100), ")")
p_comlist <- nhanes0[,c(29:38)] %>% lapply(table, nhanes0$id) %>% lapply(prop.test) %>%
  sapply(extract, "p.value") %>% as.numeric %>% format.pval(eps = .001, digits = 3) 

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
p6<- format.pval(wilcox.test(nhanes0$permth_int~nhanes0$id)$p.value, eps = .001, digits = 3)

## Anthropometric measurements ##
#Height Women
hei_w0 <- round(summary(nhanes_women$Height)[3],2)
hei_w0.1 <- round(summary(nhanes_women$Height)[2],2)
hei_w0.3 <- round(summary(nhanes_women$Height)[5],2)
hei_w1 <- round(summary((nhanes_women%>% filter(id==1))$Height)[3],2)
hei_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Height)[2],2)
hei_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Height)[5],2)
hei_w2 <- round(summary((nhanes_women%>% filter(id==2))$Height)[3],2)
hei_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Height)[2],2)
hei_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Height)[5],2)
p_hei1<- format.pval(wilcox.test(nhanes_women$Height~nhanes_women$id)$p.value, eps = .001, digits = 3)
#Height Men
hei_m0 <- round(summary(nhanes_men$Height)[3],2)
hei_m0.1 <- round(summary(nhanes_men$Height)[2],2)
hei_m0.3 <- round(summary(nhanes_men$Height)[5],2)
hei_m1 <- round(summary((nhanes_men%>% filter(id==1))$Height)[3],2)
hei_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Height)[2],2)
hei_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Height)[5],2)
hei_m2 <- round(summary((nhanes_men%>% filter(id==2))$Height)[3],2)
hei_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Height)[2],2)
hei_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Height)[5],2)
p_hei2<- format.pval(wilcox.test(nhanes_men$Height~nhanes_men$id)$p.value, eps = .001, digits = 3)

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
p15<- format.pval(wilcox.test(nhanes_women$Arm_length~nhanes_women$id)$p.value, eps = .001, digits = 3)
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
p16<- format.pval(wilcox.test(nhanes_men$Arm_length~nhanes_men$id)$p.value, eps = .001, digits = 3)

#Leg length Women
leg_w0 <- round(summary(nhanes_women$Leg_length)[3],2)
leg_w0.1 <- round(summary(nhanes_women$Leg_length)[2],2)
leg_w0.3 <- round(summary(nhanes_women$Leg_length)[5],2)
leg_w1 <- round(summary((nhanes_women%>% filter(id==1))$Leg_length)[3],2)
leg_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Leg_length)[2],2)
leg_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Leg_length)[5],2)
leg_w2 <- round(summary((nhanes_women%>% filter(id==2))$Leg_length)[3],2)
leg_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Leg_length)[2],2)
leg_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Leg_length)[5],2)
p_leg1<- format.pval(wilcox.test(nhanes_women$Leg_length~nhanes_women$id)$p.value, eps = .001, digits = 3)
#Leg length Men
leg_m0 <- round(summary(nhanes_men$Leg_length)[3],2)
leg_m0.1 <- round(summary(nhanes_men$Leg_length)[2],2)
leg_m0.3 <- round(summary(nhanes_men$Leg_length)[5],2)
leg_m1 <- round(summary((nhanes_men%>% filter(id==1))$Leg_length)[3],2)
leg_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Leg_length)[2],2)
leg_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Leg_length)[5],2)
leg_m2 <- round(summary((nhanes_men%>% filter(id==2))$Leg_length)[3],2)
leg_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Leg_length)[2],2)
leg_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Leg_length)[5],2)
p_leg2<- format.pval(wilcox.test(nhanes_men$Leg_length~nhanes_men$id)$p.value, eps = .001, digits = 3)

#Weight Women
wei_w0 <- round(summary(nhanes_women$Weight)[3],2)
wei_w0.1 <- round(summary(nhanes_women$Weight)[2],2)
wei_w0.3 <- round(summary(nhanes_women$Weight)[5],2)
wei_w1 <- round(summary((nhanes_women%>% filter(id==1))$Weight)[3],2)
wei_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Weight)[2],2)
wei_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Weight)[5],2)
wei_w2 <- round(summary((nhanes_women%>% filter(id==2))$Weight)[3],2)
wei_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Weight)[2],2)
wei_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Weight)[5],2)
p_wei1<- format.pval(wilcox.test(nhanes_women$Weight~nhanes_women$id)$p.value, eps = .001, digits = 3)
#Weight Men
wei_m0 <- round(summary(nhanes_men$Weight)[3],2)
wei_m0.1 <- round(summary(nhanes_men$Weight)[2],2)
wei_m0.3 <- round(summary(nhanes_men$Weight)[5],2)
wei_m1 <- round(summary((nhanes_men%>% filter(id==1))$Weight)[3],2)
wei_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Weight)[2],2)
wei_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Weight)[5],2)
wei_m2 <- round(summary((nhanes_men%>% filter(id==2))$Weight)[3],2)
wei_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Weight)[2],2)
wei_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Weight)[5],2)
p_wei2<- format.pval(wilcox.test(nhanes_men$Weight~nhanes_men$id)$p.value, eps = .001, digits = 3)

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
p7<- format.pval(wilcox.test(nhanes_women$BMI~nhanes_women$id)$p.value, eps = .001, digits = 3)
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
p8<- format.pval(wilcox.test(nhanes_men$BMI~nhanes_men$id)$p.value, eps = .001, digits = 3)

#Waist Women
cin_w0 <- round(summary(nhanes_women$Waist)[3],2)
cin_w0.1 <- round(summary(nhanes_women$Waist)[2],2)
cin_w0.3 <- round(summary(nhanes_women$Waist)[5],2)
cin_w1 <- round(summary((nhanes_women%>% filter(id==1))$Waist)[3],2)
cin_w1.1 <- round(summary((nhanes_women%>% filter(id==1))$Waist)[2],2)
cin_w1.3 <- round(summary((nhanes_women%>% filter(id==1))$Waist)[5],2)
cin_w2 <- round(summary((nhanes_women%>% filter(id==2))$Waist)[3],2)
cin_w2.1 <- round(summary((nhanes_women%>% filter(id==2))$Waist)[2],2)
cin_w2.3 <- round(summary((nhanes_women%>% filter(id==2))$Waist)[5],2)
p_cin1<- format.pval(wilcox.test(nhanes_women$Waist~nhanes_women$id)$p.value, eps = .001, digits = 3)
#Waist Men
cin_m0 <- round(summary(nhanes_men$Waist)[3],2)
cin_m0.1 <- round(summary(nhanes_men$Waist)[2],2)
cin_m0.3 <- round(summary(nhanes_men$Waist)[5],2)
cin_m1 <- round(summary((nhanes_men%>% filter(id==1))$Waist)[3],2)
cin_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Waist)[2],2)
cin_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Waist)[5],2)
cin_m2 <- round(summary((nhanes_men%>% filter(id==2))$Waist)[3],2)
cin_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Waist)[2],2)
cin_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Waist)[5],2)
p_cin2<- format.pval(wilcox.test(nhanes_men$Waist~nhanes_men$id)$p.value, eps = .001, digits = 3)

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
p9<- format.pval(wilcox.test(nhanes_women$ICE~nhanes_women$id)$p.value, eps = .001, digits = 3)
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
p10<- format.pval(wilcox.test(nhanes_men$ICE~nhanes_men$id)$p.value, eps = .001, digits = 3)

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
p19<- format.pval(wilcox.test(nhanes_women$Subscapular_skinfold~nhanes_women$id)$p.value, eps = .001, digits = 3)
#Subscapular skinfold Men
subs_m0 <- round(summary(nhanes_men$Subscapular_skinfold)[3],2)
subs_m0.1 <- round(summary(nhanes_men$Subscapular_skinfold)[2],2)
subs_m0.3 <- round(summary(nhanes_men$Subscapular_skinfold)[5],2)
subs_m1 <- round(summary((nhanes_men%>% filter(id==1))$Subscapular_skinfold)[3],2)
subs_m1.1 <- round(summary((nhanes_men%>% filter(id==1))$Subscapular_skinfold)[2],2)
subs_m1.3 <- round(summary((nhanes_men%>% filter(id==1))$Subscapular_skinfold)[5],2)
subs_m2 <- round(summary((nhanes_men%>% filter(id==2))$Subscapular_skinfold)[3],2)
subs_m2.1 <- round(summary((nhanes_men%>% filter(id==2))$Subscapular_skinfold)[2],2)
subs_m2.3 <- round(summary((nhanes_men%>% filter(id==2))$Subscapular_skinfold)[5],2)
p20<- format.pval(wilcox.test(nhanes_men$Subscapular_skinfold~nhanes_men$id)$p.value, eps = .001, digits = 3)

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
p17<- format.pval(wilcox.test(nhanes_women$Triceps_skinfold~nhanes_women$id)$p.value, eps = .001, digits = 3)
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
p18<- format.pval(wilcox.test(nhanes_men$Triceps_skinfold~nhanes_men$id)$p.value, eps = .001, digits = 3)

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
p13<- format.pval(wilcox.test(nhanes_women$Arm_circumference~nhanes_women$id)$p.value, eps = .001, digits = 3)
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
p14<- format.pval(wilcox.test(nhanes_men$Arm_circumference~nhanes_men$id)$p.value, eps = .001, digits = 3)

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
p11<- format.pval(wilcox.test(nhanes_women$Thigh_circumference~nhanes_women$id)$p.value, eps = .001, digits = 3)
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
p12<- format.pval(wilcox.test(nhanes_men$Thigh_circumference~nhanes_men$id)$p.value, eps = .001, digits = 3)

#Table
t0_N <- c("Female (%)","Age (years)", "Non-Hispanic White (%)", "Non-Hispanic Black (%)", "Mexican-American (%)", "PhenoAge (years)",
          ">=1 comorbidity (%)", "Arterial hypertension", "Arthritis", "Asthma", "Diabetes Mellitus", "Bronchitis", "Malignancy", "Heart attack",
          "Heart failure", "Stroke", "Emphysema", "Mortality (%)", "Follow-up (months)", "Height Females (cm)", "Height Males (cm)",
          "Arm length Females (cm)", "Arm length Males (cm)", "Leg length Females (cm)", "Leg length Males (cm)",
          "Weight Females (kg)", "Weight Males (kg)", "BMI Females (kg/m2)", "BMI Males (kg/m2)",
          "Waist circumference Females (cm)", "Waist circumference Males (cm)", "WHtR Females", "WHtR Males",
          "Subscapular skinfold Females (cm)", "Subscapcular skinfold Males (mm)",
          "Triceps skinfold Females (mm)", "Triceps skinfold Males (mm)",
          "Arm circumference Females (cm)", "Arm circumference Males (cm)",
          "Thigh circumference Females (cm)", "Thigh circumference Males (cm)")

t0_O <- c(paste0(Male0," (",pMale0,")"), paste0(Age0," (",Age0.1,"-",Age0.3,")"),
          paste0(WEth0," (",WpEth0,")"), paste0(BEth0," (",BpEth0,")"), paste0(MEth0," (",MpEth0,")"),
          paste0(PhenoAge0," (",PhenoAge0.1,"-",PhenoAge0.3,")"), paste0(Comorb0," (",pComorb0,")"),
          comorb_list0, paste0(Mort0," (",pMort0,")"),paste0(fut0," (",fut0.1,"-",fut0.3,")"),
          paste0(hei_w0," (",hei_w0.1,"-",hei_w0.3,")"),paste0(hei_m0," (",hei_m0.1,"-",hei_m0.3,")"),
          paste0(arml_w0," (",arml_w0.1,"-",arml_w0.3,")"), paste0(arml_m0," (",arml_m0.1,"-",arml_m0.3,")"),
          paste0(leg_w0," (",leg_w0.1,"-",leg_w0.3,")"),paste0(leg_m0," (",leg_m0.1,"-",leg_m0.3,")"),
          paste0(wei_w0," (",wei_w0.1,"-",wei_w0.3,")"),paste0(wei_m0," (",wei_m0.1,"-",wei_m0.3,")"),
          paste0(bmi_w0," (",bmi_w0.1,"-",bmi_w0.3,")"),paste0(bmi_m0," (",bmi_m0.1,"-",bmi_m0.3,")"),
          paste0(cin_w0," (",cin_w0.1,"-",cin_w0.3,")"),paste0(cin_m0," (",cin_m0.1,"-",cin_m0.3,")"),
          paste0(ice_w0," (",ice_w0.1,"-",ice_w0.3,")"),paste0(ice_m0," (",ice_m0.1,"-",ice_m0.3,")"),
          paste0(subs_w0," (",subs_w0.1,"-",subs_w0.3,")"), paste0(subs_m0," (",subs_m0.1,"-",subs_m0.3,")"),
          paste0(tric_w0," (",tric_w0.1,"-",tric_w0.3,")"), paste0(tric_m0," (",tric_m0.1,"-",tric_m0.3,")"),
          paste0(armc_w0," (",armc_w0.1,"-",armc_w0.3,")"), paste0(armc_m0," (",armc_m0.1,"-",armc_m0.3,")"),
          paste0(thigh_w0," (",thigh_w0.1,"-",thigh_w0.3,")"), paste0(thigh_m0," (",thigh_m0.1,"-",thigh_m0.3,")"))

t0_3 <- c(paste0(Male1," (",pMale1,")"), paste0(Age1," (",Age1.1,"-",Age1.3,")"),
          paste0(WEth1," (",WpEth1,")"), paste0(BEth1," (",BpEth1,")"), paste0(MEth1," (",MpEth1,")"),
          paste0(PhenoAge1," (",PhenoAge1.1,"-",PhenoAge1.3,")"), paste0(Comorb1," (",pComorb1,")"),
          comorb_list1, paste0(Mort1," (",pMort1,")"),paste0(fut1," (",fut1.1,"-",fut1.3,")"),
          paste0(hei_w1," (",hei_w1.1,"-",hei_w1.3,")"),paste0(hei_m1," (",hei_m1.1,"-",hei_m1.3,")"),
          paste0(arml_w1," (",arml_w1.1,"-",arml_w1.3,")"), paste0(arml_m1," (",arml_m1.1,"-",arml_m1.3,")"),
          paste0(leg_w1," (",leg_w1.1,"-",leg_w1.3,")"),paste0(leg_m1," (",leg_m1.1,"-",leg_m1.3,")"),
          paste0(wei_w1," (",wei_w1.1,"-",wei_w1.3,")"),paste0(wei_m1," (",wei_m1.1,"-",wei_m1.3,")"),
          paste0(bmi_w1," (",bmi_w1.1,"-",bmi_w1.3,")"),paste0(bmi_m1," (",bmi_m1.1,"-",bmi_m1.3,")"),
          paste0(cin_w1," (",cin_w1.1,"-",cin_w1.3,")"),paste0(cin_m1," (",cin_m1.1,"-",cin_m1.3,")"),
          paste0(ice_w1," (",ice_w1.1,"-",ice_w1.3,")"),paste0(ice_m1," (",ice_m1.1,"-",ice_m1.3,")"),
          paste0(subs_w1," (",subs_w1.1,"-",subs_w1.3,")"), paste0(subs_m1," (",subs_m1.1,"-",subs_m1.3,")"),
          paste0(tric_w1," (",tric_w1.1,"-",tric_w1.3,")"), paste0(tric_m1," (",tric_m1.1,"-",tric_m1.3,")"),
          paste0(armc_w1," (",armc_w1.1,"-",armc_w1.3,")"), paste0(armc_m1," (",armc_m1.1,"-",armc_m1.3,")"),
          paste0(thigh_w1," (",thigh_w1.1,"-",thigh_w1.3,")"), paste0(thigh_m1," (",thigh_m1.1,"-",thigh_m1.3,")"))

t0_4 <- c(paste0(Male2," (",pMale2,")"), paste0(Age2," (",Age2.1,"-",Age2.3,")"),
          paste0(WEth2," (",WpEth2,")"), paste0(BEth2," (",BpEth2,")"), paste0(MEth2," (",MpEth2,")"),
          paste0(PhenoAge2," (",PhenoAge2.1,"-",PhenoAge2.3,")"), paste0(Comorb2," (",pComorb2,")"),
          comorb_list2, paste0(Mort2," (",pMort2,")"),paste0(fut2," (",fut2.1,"-",fut2.3,")"),
          paste0(hei_w2," (",hei_w2.1,"-",hei_w2.3,")"),paste0(hei_m2," (",hei_m2.1,"-",hei_m2.3,")"),
          paste0(arml_w2," (",arml_w2.1,"-",arml_w2.3,")"), paste0(arml_m2," (",arml_m2.1,"-",arml_m2.3,")"),
          paste0(leg_w2," (",leg_w2.1,"-",leg_w2.3,")"),paste0(leg_m2," (",leg_m2.1,"-",leg_m2.3,")"),
          paste0(wei_w2," (",wei_w2.1,"-",wei_w2.3,")"),paste0(wei_m2," (",wei_m2.1,"-",wei_m2.3,")"),
          paste0(bmi_w2," (",bmi_w2.1,"-",bmi_w2.3,")"),paste0(bmi_m2," (",bmi_m2.1,"-",bmi_m2.3,")"),
          paste0(cin_w2," (",cin_w2.1,"-",cin_w2.3,")"),paste0(cin_m2," (",cin_m2.1,"-",cin_m2.3,")"),
          paste0(ice_w2," (",ice_w2.1,"-",ice_w2.3,")"),paste0(ice_m2," (",ice_m2.1,"-",ice_m2.3,")"),
          paste0(subs_w2," (",subs_w2.1,"-",subs_w2.3,")"), paste0(subs_m2," (",subs_m2.1,"-",subs_m2.3,")"),
          paste0(tric_w2," (",tric_w2.1,"-",tric_w2.3,")"), paste0(tric_m2," (",tric_m2.1,"-",tric_m2.3,")"),
          paste0(armc_w2," (",armc_w2.1,"-",armc_w2.3,")"), paste0(armc_m2," (",armc_m2.1,"-",armc_m2.3,")"),
          paste0(thigh_w2," (",thigh_w2.1,"-",thigh_w2.3,")"), paste0(thigh_m2," (",thigh_m2.1,"-",thigh_m2.3,")"))

t0_P <- c(p1,p2,"","",pEthn,p3,p4,p_comlist,p5,p6, p_hei1,p_hei2,p15,p16,p_leg1,p_leg2,
          p_wei1,p_wei2,p7,p8, p_cin1,p_cin2,p9,p10, p19,p20,p17,p18, p13,p14,p11,p12)

tab0<-data.frame("Names"= t0_N, "Overall"= t0_O, "NIII"= t0_3, "NIV"= t0_4, "P"= t0_P) %>% 
  `names<-`(c("Characteristics","Overall (n=18,794)","NHANES-III (n=11,774)","NHANES-IV (n=7,020)","P-value")) %>% 
  flextable(cwidth = c(2,1.5,1.5,2,1)) %>% align(align = "center",part = "all")
save_as_docx(tab0,path="tabla1.docx")

####---- Models for anthropometric measurements ---####
### Models in men ###
gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+poly(BMI,2)+num_comorb+shape(Ethnicity),dist="gompertz",data=nhanes_men)
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
gomp1<-flexsurvreg(Surv(permth_int,mortstat)~poly(Age,3)+num_comorb+poly(BMI,2)+shape(Ethnicity),dist="gompertz",data=nhanes_women)
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


### Run Figure 2 ###
nhanes_fin<-rbind(nhanes_women, nhanes_men)
F_names <- c("Height (cm)", "Arm length (cm)","Leg length (cm)",
             "Weight (Kg)", "Body-mass index (kg/m2)",
             "Waist-to-height ratio", "Subscapular skinfold (cm)","Triceps skinfold (cm)",
             "Thigh circumference (cm)", "Arm circumference (cm)")

#Height
F_lims <- c(140,193); f1A <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=Height, y=risk9, col=Sex, fill=Sex))+geom_smooth(alpha=0.3,fullrange=F)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=Height, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[1])+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "none")))

#Arm_length
F_lims <- c(25.5,45); f1B <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=Arm_length, y=risk6, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=Arm_length, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[2])+scale_y_continuous(name="Density", breaks = c(0,.18))+theme(legend.position = "none")))

#Leg_length
F_lims <- c(28,51); f1C <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=Leg_length, y=risk8, col=Sex, fill=Sex))+geom_smooth(alpha=0.3,fullrange=F)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=Leg_length, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[3])+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "none")))

#Weight
F_lims <- c(35,130); f1D <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=Weight, y=risk10, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=Weight, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[4])+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "none")))

#BMI
F_lims <- c(16,50); f1E <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=BMI, y=risk1, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=BMI, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[5])+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "none")))

#ICE
F_lims <- c(0.35,0.8); f1F <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=ICE, y=risk2, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=ICE, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[6])+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "none")))

#Subscapular_skinfold
F_lims <- c(0,45); f1G <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=Subscapular_skinfold, y=risk3, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=Subscapular_skinfold, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[7])+scale_y_continuous(name="Density", n.breaks = 2)+theme(legend.position = "none")))

#Triceps_skinfold
F_lims <- c(0,43); f1H <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=Triceps_skinfold, y=risk4, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=Triceps_skinfold, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[8])+scale_y_continuous(name="Density", breaks = c(0,.08))+theme(legend.position = "none")))

#Arm_circumference
F_lims <- c(21,45); f1I <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=Arm_circumference, y=risk7, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=Arm_circumference, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[10])+scale_y_continuous(name="Density", breaks = c(0,.1))+theme(legend.position = "none")))
  
#Thigh_circumference
F_lims <- c(36,71); f1J <- ggarrange(
  ncol = 1, nrow = 2, heights = c(1,0.35),
  (ggplot(nhanes_fin, aes(x=Thigh_circumference, y=risk5, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     scale_x_continuous(name="", limits = NULL)+scale_y_continuous(name="10-year mortality risk", limits=NULL)+theme(legend.position = "none")+
     scale_fill_manual(values=c("#994455","#6699CC")) + coord_cartesian(ylim = c(0,0.8), xlim = F_lims)),
  (ggplot(nhanes_fin, aes(x=Thigh_circumference, fill=Sex))+geom_density(alpha=0.75,size=0.3)+scale_fill_manual(values=c("#994455","#6699CC"))+theme_pubclean()+
     coord_cartesian(xlim = F_lims)+xlab(F_names[9])+scale_y_continuous(name="Density", breaks = c(0,.07))+theme(legend.position = "none")))


F_legend <- get_legend(ggplot(nhanes_fin, aes(x=Height, y=risk9, col=Sex, fill=Sex))+geom_smooth(alpha=0.3)+theme_pubclean()+
                         scale_color_manual(values=c("#994455","#6699CC"))+scale_fill_manual(values=c("#994455","#6699CC"))+
                         theme(legend.key.size = unit(0.8, 'cm'), legend.text = element_text(size=11.5), legend.title = element_text(size=14)))
fig1<-ggarrange(f1A, f1B, f1C, f1D, f1E, f1F, f1G, f1H, f1I, f1J, ncol=5, nrow=2, labels=letters[1:10], common.legend=T, legend="bottom", legend.grob=F_legend)

ggsave(fig1,filename = "Figure1.jpg", width = 43.75*0.9, height = 24*0.8, units=c("cm"), dpi = 300, limitsize = FALSE)

#Submission
ggsave(fig1,filename = "Submission/Figures/Figure_2.pdf",
       width = 39.375, height = 19.2, units=c("cm"), dpi = 600, limitsize = FALSE)

#Graphical abstract
(ggplot(nhanes_fin, aes(x=Arm_circumference, y=risk7, col=Sex, fill=Sex))+
    geom_smooth(alpha=0.3)+scale_color_manual(values=c("#994455","#6699CC"))+
    theme_void() + coord_cartesian(ylim = c(0,0.8), xlim = c(21,45))) + theme(legend.position = "none")


####---- Development of AnthropoAge ----####
### Training and validation ###
train_nhanes <- nhanes0[nhanes0$id==1,]
test_nhanes <- nhanes0[nhanes0$id==2,]

###Gompertz models###
options(scipen=10)
gomp1aM<-flexsurvreg(Surv(permth_int,mortstat) ~ Age + poly(tr_ice,2) +
                       tr_armc + poly(tr_thigh,2) + shape(Ethnicity),dist="gompertz",
                     data=train_nhanes%>%filter(Sex=="Men")); gomp1aM; BIC(gomp1aM)

gomp1aF<-flexsurvreg(Surv(permth_int,mortstat)~Age + tr_weight + tr_ice + 
                        tr_subs + tr_tric + poly(tr_thigh,2) + shape(Ethnicity), dist="gompertz",
                     data=train_nhanes%>%filter(Sex=="Women")); gomp1aF;BIC(gomp1aF)

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


### Supplementary table 3 ###
coef1M<-round(coef(gomp1aM),4)
confint1M<-round(confint(gomp1aM)[,1],4)
confint2M<-round(confint(gomp1aM)[,2],4)
bic1<-round(BIC(gomp1aM),1)

coef1F<-round(coef(gomp1aF),4)
confint1F<-round(confint(gomp1aF)[,1],4)
confint2F<-round(confint(gomp1aF)[,2],4)
bic2<-round(BIC(gomp1aF),1)

tab2 <- data.frame(
  "Model"= c("AnthropoAge Males", paste0("BIC ",bic1), rep(" ",9),
             "AnthropoAge Females", paste0("BIC ",bic2), rep(" ",10)),
  "Parameter"=c("Shape", "Rate", "Chronological Age", "WHtR (OP,1)", "WHtR (OP,2)",
                "Arm circumference", "Thigh circumference (OP,1)", "Thigh circumference (OP,2)",
                "Shape (Non-Hispanic Black)", "Shape (Mexican-American)", "Shape (Other race/ethnicity)",
                "Shape", "Rate", "Chronological Age", "Weight", "WHtR",
                "Subscapular skinfold", "Triceps skinfold", "Thigh circumference (OP,1)", "Thigh circumference (OP,2)",
                "Shape (Non-Hispanic Black)", "Shape (Mexican-American)", "Shape (Other race/ethnicity)"),
  "B-coefficient"=c(coef1M,coef1F), "Lower 95%CI"=c(confint1M,confint1F), "Upper 95%CI"=c(confint2M,confint1F))
rownames(tab2)<-NULL; tab2<- tab2 %>% `names<-`(c("Model","Parameter","B-coefficient","Lower 95%CI","Upper 95%CI")) %>% 
  flextable(cwidth = c(1.5,2.5,1,1,1)) %>% align(align = "center",part = "all")
save_as_docx(tab2,path="tabla2.docx")


### Training dataset ###
p1F<-predict(gomp1aF, newdata =train_nhanes %>%filter(Sex=="Women"),type="survival", ci=F, times = c(120))
p1M<-predict(gomp1aM, newdata =train_nhanes %>%filter(Sex=="Men"),type="survival", ci=F, times = c(120))
train_nhanes$pred[train_nhanes$Sex=="Women"]<-as.numeric(1-p1F$.pred)
train_nhanes$pred[train_nhanes$Sex=="Men"]<-as.numeric(1-p1M$.pred)
train_nhanes$AnthropoAge[train_nhanes$Sex=="Women"]<-(log(-sW*log(1-train_nhanes$pred[train_nhanes$Sex=="Women"]))-b0W)/b1W
train_nhanes$AnthropoAge[train_nhanes$Sex=="Men"]<-(log(-sM*log(1-train_nhanes$pred[train_nhanes$Sex=="Men"]))-b0M)/b1M
# Accelerated metrics #
m1F<-lm(AnthropoAge~Age, data=train_nhanes %>% filter(Sex=="Women"))
train_nhanes$AnthropoAgeAccel[train_nhanes$Sex=="Women"]<-m1F$residuals
m1M<-lm(AnthropoAge~Age, data=train_nhanes %>% filter(Sex=="Men"))
train_nhanes$AnthropoAgeAccel[train_nhanes$Sex=="Men"]<-m1M$residuals
m1Ph<-lm(PhenoAge~Age, data=train_nhanes %>% filter(!is.infinite(PhenoAge)))
train_nhanes$PhenoAgeAccel[!is.na(train_nhanes$PhenoAge)]<-m1Ph$residuals
# Sex differences #
tapply(train_nhanes$AnthropoAge, train_nhanes$Sex, quantile)
wilcox.test(train_nhanes$AnthropoAge~as.numeric(train_nhanes$Sex)) #AnthropoAge
tapply(train_nhanes$AnthropoAgeAccel, train_nhanes$Sex, quantile)
wilcox.test(train_nhanes$AnthropoAgeAccel~as.numeric(train_nhanes$Sex)) #AnthropoAgeAccel
tapply(train_nhanes$PhenoAge, train_nhanes$Sex, quantile, na.rm=T)
wilcox.test(train_nhanes$PhenoAge~as.numeric(train_nhanes$Sex)) #PhenoAge
tapply(train_nhanes$PhenoAgeAccel, train_nhanes$Sex, quantile, na.rm=T)
wilcox.test(train_nhanes$PhenoAgeAccel~as.numeric(train_nhanes$Sex)) #PhenoAgeAccel


### Testing dataset ###
p1F<-predict(gomp1aF, newdata = test_nhanes %>% filter(Sex=="Women") ,type="survival", ci=F, times = c(120))
p1M<-predict(gomp1aM, newdata = test_nhanes %>% filter(Sex=="Men"),type="survival", ci=F, times = c(120))
test_nhanes$pred[test_nhanes$Sex=="Women"]<-as.numeric(1-p1F$.pred)
test_nhanes$pred[test_nhanes$Sex=="Men"]<-as.numeric(1-p1M$.pred)
test_nhanes$AnthropoAge[test_nhanes$Sex=="Women"]<-(log(-sW*log(1-test_nhanes$pred[test_nhanes$Sex=="Women"]))-b0W)/b1W
test_nhanes$AnthropoAge[test_nhanes$Sex=="Men"]<-(log(-sM*log(1-test_nhanes$pred[test_nhanes$Sex=="Men"]))-b0M)/b1M
# Accelerated metrics #
m1.1F<-lm(AnthropoAge~Age, data=test_nhanes %>% filter(Sex=="Women"))
test_nhanes$AnthropoAgeAccel[test_nhanes$Sex=="Women"]<-m1.1F$residuals
m1.1M<-lm(AnthropoAge~Age, data=test_nhanes %>% filter(Sex=="Men"))
test_nhanes$AnthropoAgeAccel[test_nhanes$Sex=="Men"]<-m1.1M$residuals
m1.1Ph<-lm(PhenoAge~Age, data=test_nhanes %>% filter(!is.infinite(PhenoAge)))
test_nhanes$PhenoAgeAccel[!is.na(test_nhanes$PhenoAge)]<-m1.1Ph$residuals
# Sex differences #
tapply(test_nhanes$AnthropoAge, test_nhanes$Sex, quantile) %>% lapply(round,1)
wilcox.test(test_nhanes$AnthropoAge~as.numeric(test_nhanes$Sex)) #AnthropoAge
tapply(test_nhanes$AnthropoAgeAccel, test_nhanes$Sex, quantile) %>% lapply(round,2)
wilcox.test(test_nhanes$AnthropoAgeAccel~as.numeric(test_nhanes$Sex)) #AnthropoAgeAccel
tapply(test_nhanes$PhenoAge, test_nhanes$Sex, quantile, na.rm=T) %>% lapply(round,1)
wilcox.test(test_nhanes$PhenoAge~as.numeric(test_nhanes$Sex)) #PhenoAge
tapply(test_nhanes$PhenoAgeAccel, test_nhanes$Sex, quantile, na.rm=T) %>% lapply(round,2)
wilcox.test(test_nhanes$PhenoAgeAccel~as.numeric(test_nhanes$Sex)) #PhenoAgeAccel

### Save models as RDA files to develop ShinyApp ###
# English version #
save(gomp1aM, file="Shiny App/anthropoage/Models/1_CAnthropo_M.rda")
save(gomp1aF, file="Shiny App/anthropoage/Models/2_CAnthropo_F.rda")
save(gomp1bM, file="Shiny App/anthropoage/Models/3_CAge_M.rda")
save(gomp1bF, file="Shiny App/anthropoage/Models/4_CAge_F.rda")
save(m1M, file="Shiny App/anthropoage/Models/5_CAccel_M.rda")
save(m1F, file="Shiny App/anthropoage/Models/6_CAccel_F.rda")
save(m1Ph, file="Shiny App/anthropoage/Models/13_PhenoAgeAccel.rda")
# Spanish version #
save(gomp1aM, file="Shiny App/anthropoage_es/Models/1_CAnthropo_M.rda")
save(gomp1aF, file="Shiny App/anthropoage_es/Models/2_CAnthropo_F.rda")
save(gomp1bM, file="Shiny App/anthropoage_es/Models/3_CAge_M.rda")
save(gomp1bF, file="Shiny App/anthropoage_es/Models/4_CAge_F.rda")
save(m1M, file="Shiny App/anthropoage_es/Models/5_CAccel_M.rda")
save(m1F, file="Shiny App/anthropoage_es/Models/6_CAccel_F.rda")
save(m1Ph, file="Shiny App/anthropoage_es/Models/13_PhenoAgeAccel.rda")


### Figure 2 ###
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

fig2<-ggarrange(annotate_figure(ggarrange(f1,f2,f5,f6, labels=letters[c(1,2,5,6)], common.legend = TRUE), top=text_grob("Training cohort", face = "bold", size=12)),
                annotate_figure(ggarrange(f3,f4,f7,f8, labels=letters[c(3,4,7,8)], common.legend = TRUE), top=text_grob("Validation cohort", face = "bold", size=12)), nrow=1, ncol=2)

ggsave(file = "Figure2.jpg", fig2, bg = "transparent", width = 27,  height = 15, units=c("cm"), dpi = 500, limitsize = FALSE)
#Submission
ggsave(fig2, file = "Submission/Figures/Figure_3.pdf", bg = "transparent", width = 27,  height = 15, units=c("cm"), dpi = 600, limitsize = FALSE)
#Graphical abstract
ggplot(test_nhanes, aes(x=AnthropoAge, y=Age, col=Sex))+theme_void()+geom_jitter(size=0.25,alpha=0.6)+
  geom_smooth(method="lm")+scale_color_manual(values=c("#994455","#8DCBE4"))+
  theme(legend.position = "none")


### Correlation between anthropometric variables ###
rphA<-nhanes0 %>% filter(Sex=="Women") %>% select(
  Height, Arm_length, Leg_length, Weight, BMI, ICE, Subscapular_skinfold, Triceps_skinfold, Arm_circumference, Thigh_circumference
  ) %>% as.matrix %>% rcorr(type = "spearman"); rphB<-nhanes0 %>% filter(Sex=="Men") %>% select(
  Height, Arm_length, Leg_length, Weight, BMI, ICE, Subscapular_skinfold, Triceps_skinfold, Arm_circumference, Thigh_circumference
    ) %>% as.matrix %>% rcorr(type = "spearman"); F_names <- c(
          "Height", "Arm Length", "Leg length", "Weight", "BMI", "WHtR",
          "Subscapular\nskinfold", "Triceps\nskinfold", "Arm\ncircumference", "Thigh\ncircumference"
          ); rownames(rphA$r) <- F_names; colnames(rphA$r) <- F_names; rownames(rphB$r) <- F_names; colnames(rphB$r) <- F_names

corrplot_full<-function(x,y){print(x); print(y)}
FigSX_A <- as.ggplot(
  ~corrplot_full(corrplot(rphA$r,method="number",type="lower",add=F,p.mat=rphA$p,sig.level=.05/50,tl.pos="tl", cl.pos="r"),
                 corrplot(rphA$r, method="circle", type="upper", add=T, p.mat=rphA$p, sig.level=.05/50, tl.pos = "n"))) %>% 
  annotate_figure(top = text_grob("Correlations in women", face = "bold", size = 15, vjust = 1.5))
FigSX_B <- as.ggplot(
  ~corrplot_full(corrplot(rphB$r,method="number",type="lower",add=F,p.mat=rphB$p,sig.level=.05/50,tl.pos="tl", cl.pos="r"),
                 corrplot(rphB$r, method="circle", type="upper", add=T, p.mat=rphB$p, sig.level=.05/50, tl.pos = "n"))) %>% 
  annotate_figure(top = text_grob("Correlations in men", face = "bold", size = 15, vjust = 1.5))

ggarrange(FigSX_A,FigSX_B, labels = letters[1:2], nrow=1, ncol=2) %>%
  ggsave(file = "SuppFig2.jpg", bg = "transparent", width = 37.1475,  height = 19.52625, units=c("cm"), dpi = 300, limitsize = FALSE)


### Multicolinearity within each model ###
coxph(Surv(permth_int,mortstat)~Age+poly(tr_thigh,2)+tr_armc+
        poly(tr_ice,2)+strata(Ethnicity), data=train_nhanes%>%filter(Sex=="Men")) %>% vif()
coxph(Surv(permth_int,mortstat)~Age+tr_weight+poly(tr_thigh,2)+
        tr_ice+tr_tric+tr_subs+strata(Ethnicity), data=train_nhanes%>%filter(Sex=="Women")) %>% vif()

coxph(Surv(permth_int,mortstat)~Age+poly(tr_imc,2)+tr_ice+strata(Ethnicity),
            data=train_nhanes%>%filter(Sex=="Men")) %>% vif
coxph(Surv(permth_int,mortstat)~Age+poly(tr_imc,2)+tr_ice+strata(Ethnicity),
            data=train_nhanes%>%filter(Sex=="Women")) %>% vif

####---- Simplified AnthropoAge (S-AnthropoAge) ----####
train_nhanes1 <- nhanes0[nhanes0$id==1,]
test_nhanes1 <- nhanes0[nhanes0$id==2,]

#Gompertz w/anthropometric variables
gomp1aM1<-flexsurvreg(Surv(permth_int,mortstat)~Age+poly(tr_imc,2)+tr_ice+shape(Ethnicity),
                      dist="gompertz",data=train_nhanes1%>%filter(Sex=="Men")); gomp1aM1; BIC(gomp1aM1)
gomp1aF1<-flexsurvreg(Surv(permth_int,mortstat)~Age+poly(tr_imc,2)+tr_ice+shape(Ethnicity),
                      dist="gompertz",data=train_nhanes1%>%filter(Sex=="Women")); gomp1aF1; BIC(gomp1aF1)


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
m2F<-lm(AnthropoAge2~Age, data=train_nhanes1 %>% filter(Sex=="Women"))
train_nhanes1$AnthropoAgeAccel2[train_nhanes1$Sex=="Women"]<-m2F$residuals
m2M<-lm(AnthropoAge2~Age, data=train_nhanes1 %>% filter(Sex=="Men"))
train_nhanes1$AnthropoAgeAccel2[train_nhanes1$Sex=="Men"]<-m2M$residuals
#S-AnthropoAgeAccel (test)
m2.1F<-lm(AnthropoAge2~Age, data=test_nhanes1 %>% filter(Sex=="Women"))
test_nhanes1$AnthropoAgeAccel2[test_nhanes1$Sex=="Women"]<-m2.1F$residuals
m2.1M<-lm(AnthropoAge2~Age, data=test_nhanes1 %>% filter(Sex=="Men"))
test_nhanes1$AnthropoAgeAccel2[test_nhanes1$Sex=="Men"]<-m2.1M$residuals

### Save models for development of ShinyApp ###
# English version #
save(gomp1aM1, file="Shiny App/anthropoage/Models/7_SAnthropo_M.rda")
save(gomp1aF1, file="Shiny App/anthropoage/Models/8_SAnthropo_F.rda")
save(gomp1bM1, file="Shiny App/anthropoage/Models/9_SAge_M.rda")
save(gomp1bF1, file="Shiny App/anthropoage/Models/10_SAge_F.rda")
save(m2M, file="Shiny App/anthropoage/Models/11_SAccel_M.rda")
save(m2F, file="Shiny App/anthropoage/Models/12_SAccel_F.rda")
# Spanish version #
save(gomp1aM1, file="Shiny App/anthropoage_es/Models/7_SAnthropo_M.rda")
save(gomp1aF1, file="Shiny App/anthropoage_es/Models/8_SAnthropo_F.rda")
save(gomp1bM1, file="Shiny App/anthropoage_es/Models/9_SAge_M.rda")
save(gomp1bF1, file="Shiny App/anthropoage_es/Models/10_SAge_F.rda")
save(m2M, file="Shiny App/anthropoage_es/Models/11_SAccel_M.rda")
save(m2F, file="Shiny App/anthropoage_es/Models/12_SAccel_F.rda")


### Table 2B  ###
coef1M1<-round(coef(gomp1aM1),4)
confint1M1<-round(confint(gomp1aM1)[,1],4)
confint2M1<-round(confint(gomp1aM1)[,2],4)
bic11<-round(BIC(gomp1aM1),1)

coef1F1<-round(coef(gomp1aF1),4)
confint1F1<-round(confint(gomp1aF1)[,1],4)
confint2F1<-round(confint(gomp1aF1)[,2],4)
bic21<-round(BIC(gomp1aF1),1)

tab2.1 <- data.frame(
  "Model"=c("S-AnthropoAge Males", paste0("BIC ",bic11), c(rep(" ",7)),
            "S-AnthropoAge Females", paste0("BIC ",bic21), rep(" ",7)),
  "Parameter"=c("Shape", "Rate", "Chronological Age", "BMI (OP,1)","BMI (OP,2)", "WHtR",
                "Shape (Non-Hispanic Black)", "Shape (Mexican American)", "Shape (Other race/ethnicity)",
                "Shape", "Rate", "Chronological Age", "BMI (OP,1)","BMI (OP,2)", "WHtR",
                "Shape (Non-Hispanic Black)", "Shape (Mexican American)", "Shape (Other race/ethnicity)"),
  "B-coefficient"=c(coef1M1,coef1F1), "Lower 95%CI"=c(confint1M1,confint1F1), "Upper 95%CI"=c(confint2M1,confint1F1)
  ); rownames(tab2.1)<-NULL

tab2.1 %>% `names<-`(c("Model","Parameter","B-coefficient","Lower 95%CI","Upper 95%CI")) %>% 
  flextable(cwidth = c(1.5,2.5,1,1,1)) %>% align(align = "center",part = "all") %>% 
  save_as_docx(path="tabla2B.docx")


####---- Bland-Altman ----####
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


#Bland-Altman AnthropoAge vs S-AnthropoAge
stats1<-blandr::blandr.statistics(nhanes0$AnthropoAge, nhanes0$AnthropoAge2, sig.level = 0.95)
s<-round(stats1$bias,3); sU<-round(stats1$biasUpperCI,3); sL<-round(stats1$biasLowerCI,3)
ann1 <- paste0("Bias = ",s," (",sL," - ",sU,")")

antro <- nhanes0%>%dplyr::select(AnthropoAge, AnthropoAge2)
ICC <- irr::icc(antro, model = "twoway", type = "agreement", unit = "single")
i1 <- round(ICC$value,4); i2 <- format.pval(ICC$p.value,eps="0.001")
ann2 <- paste0("ICC = ", i1, ", 95%CI ", (ICC$lbound %>% round(4)), "-", (ICC$ubound %>% round(4)), "," ," p", i2)

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
ggsave(file = "SuppFig3.jpg", supp2, bg = "transparent", width = 35, height = 15, units=c("cm"), dpi = 300, limitsize = FALSE)



### Race/Ethnicity ###
train_nhanes <- nhanes0 %>% filter(id==1); test_nhanes <- nhanes0 %>% filter(id==2)

s3a <- ggplot(test_nhanes, aes(y=PhenoAge, x=Age, col=Ethnicity, linetype=Ethnicity)) + geom_jitter(alpha=0.2, color="gray") + geom_smooth(method="lm") +
  theme_classic() + ylab("PhenoAge (years)") + xlab("CA (years)") + theme(legend.position = "bottom") +
  scale_color_manual(values=c("#6699CC","#994455","#A98203","#364B9A")) + scale_linetype_manual(values=c("solid","11","dashed","longdash")) +
  theme(axis.title = element_text(size=13), axis.text = element_text(size=9.5)) + theme(legend.position = "none")

s3b <- ggplot(test_nhanes, aes(x=Ethnicity, y=PhenoAgeAccel, fill=Ethnicity)) + geom_boxplot() +
  theme_classic() + ylab("PhenoAgeAccel (years)") + xlab("") + theme(legend.position = "none") +
  scale_fill_manual(values=c("#6699CC","#994455","#A98203","#364B9A")) + ggpubr::stat_compare_means(size=4.5,label.x=2, label.y=50) + 
  theme(axis.title = element_text(size=13), axis.text = element_text(size=11)) + ylim(c(-20,54)) + ggbreak::scale_y_break(c(22,46),expand = T)

s3c <- ggplot(test_nhanes, aes(y=AnthropoAge, x=Age, col=Ethnicity, linetype=Ethnicity)) + geom_jitter(alpha=0.2, color="gray") + geom_smooth(method="lm") +
  theme_classic() + ylab("AnthropoAge (years)") + xlab("CA (years)") + theme(legend.position = "bottom") +
  scale_color_manual(values=c("#6699CC","#994455","#A98203","#364B9A")) + scale_linetype_manual(values=c("solid","11","dashed","longdash")) +
  theme(axis.title = element_text(size=13), axis.text = element_text(size=9.5)) + theme(legend.position = "none")

s3d <- ggplot(test_nhanes, aes(x=Ethnicity, y=AnthropoAgeAccel, fill=Ethnicity)) + geom_boxplot() +
  theme_classic() + ylab("AnthropoAgeAccel (years)") + xlab("") + theme(legend.position = "none") +
  scale_fill_manual(values=c("#6699CC","#994455","#A98203","#364B9A")) + ggpubr::stat_compare_means(size=4.5,label.x=2, label.y=50) + 
  theme(axis.title = element_text(size=13), axis.text = element_text(size=11)) + ylim(c(-20,54)) + ggbreak::scale_y_break(c(22,46),expand = T)

s3e <- ggplot(test_nhanes, aes(y=AnthropoAge2, x=Age, col=Ethnicity, linetype=Ethnicity)) + geom_jitter(alpha=0.2, color="gray") + geom_smooth(method="lm") +
  theme_classic() + ylab("S-AnthropoAge (years)") + xlab("CA (years)") + theme(legend.position = "bottom") +
  scale_color_manual(values=c("#6699CC","#994455","#A98203","#364B9A")) + scale_linetype_manual(values=c("solid","11","dashed","longdash")) +
  theme(axis.title = element_text(size=13), axis.text = element_text(size=9.5)) +
  theme(legend.position="bottom", legend.key.size = unit(8,"mm"), legend.text = element_text(size=12), legend.title = element_text(size=15))

s3f <- ggplot(test_nhanes, aes(x=Ethnicity, y=AnthropoAgeAccel2, fill=Ethnicity)) + geom_boxplot() +
  theme_classic() + ylab("S-AnthropoAgeAccel (years)") + xlab("") + theme(legend.position = "none") +
  scale_fill_manual(values=c("#6699CC","#994455","#A98203","#364B9A")) + ggpubr::stat_compare_means(size=4.5,label.x=2, label.y=50) + 
  theme(axis.title = element_text(size=13), axis.text = element_text(size=11)) + ylim(c(-20,54)) + ggbreak::scale_y_break(c(22,46),expand = T)

supp3a <- ggarrange(s3a,s3c,s3e, labels = letters[1:3], common.legend = TRUE, legend = "bottom", legend.grob = get_legend(s3e), ncol = 3)
supp3b <- ggarrange(print(s3b),print(s3d),print(s3f), labels = letters[4:6], common.legend = TRUE, legend = "none", ncol = 3)
supp3 <- ggarrange(supp3a, supp3b, ncol = 1)
ggsave(file = "SuppFig4.jpg", supp3, bg = "transparent", width = 45*0.8, height = 30*0.8, units=c("cm"), dpi = 300, limitsize = FALSE)

with(test_nhanes, tapply(PhenoAgeAccel, Ethnicity, median)) %>% round(2)
with(test_nhanes, tapply(AnthropoAgeAccel, Ethnicity, median)) %>% round(2)
with(test_nhanes, tapply(AnthropoAgeAccel2, Ethnicity, median)) %>% round(2)
table(test_nhanes$Ethnicity)

##TEST##
with(test_nhanes, tapply(PhenoAgeAccel, Ethnicity, quantile)) %>% lapply(round,2)
with(test_nhanes, kruskal.test(PhenoAgeAccel ~ Ethnicity)) #PhenoAge
with(test_nhanes, tapply(AnthropoAgeAccel, Ethnicity, quantile)) %>% lapply(round,2)
with(test_nhanes, kruskal.test(AnthropoAgeAccel ~ Ethnicity)) #AnthropoAge
with(test_nhanes, tapply(AnthropoAgeAccel2, Ethnicity, quantile)) %>% lapply(round,2)
with(test_nhanes, kruskal.test(AnthropoAgeAccel2 ~ Ethnicity)) #S-AnthropoAge


####---- ROC curves ----####
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

#roc.test(r1, r2, method = "boot")
#roc.test(r3, r2, method = "boot") #Anthropo, Pheno, Age
#roc.test(r4, r2, method = "boot"); roc.test(r5, r2, method = "boot")
#roc.test(r6, r2, method = "boot"); roc.test(r7, r2, method = "boot")
#roc.test(r8, r1, method = "boot"); roc.test(r8, r2, method = "boot") #S-Anthropo, Pheno, Age

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

fig3a <- ggroc(roc.list, size=1.5*0.65, aes("color","linetype"))+
  theme(legend.position="none")+ theme_minimal() + ggtitle("Training cohort") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid", size=0.2)+
  xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)") +
  labs(col="Scores", linetype="Scores")+ theme_pubclean()+
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size =17*0.65, face = 'plain'), axis.title.y = element_text(size=17*0.65, face = 'plain'),
        axis.text.x = element_text(size = 14*0.65, face ='plain'), axis.text.y = element_text(size=14*0.65, face ='plain'),
        legend.text = element_text(size = 17*0.65))+
  annotation_custom(tableGrob(table1, theme=ttheme_default(base_size = 9.75, padding = unit(c(2.5,2.5),"mm") )),
                    xmin=-0.5, xmax=0, ymin=0.2, ymax=0.4)+
  scale_color_manual(values=c("#994455","#6699CC","#EECC66","#364B9A","#F67B4E","#EE99AA","#C2E4EF"))+
  theme(legend.key.size = unit(12*0.65,"mm"))+scale_linetype_manual(values=c("solid","11","dashed","longdash","dotted","twodash","solid"))+
  theme(text=element_text(size=17*0.65),axis.text=element_text(size=15*0.65),legend.text=element_text(size=17*0.65),
        legend.position = "bottom", plot.title = element_text(size=18*0.65, face="bold", hjust=0.5))


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

#roc.test(r2, r1, method = "boot")
#roc.test(r3, r2, method = "boot")
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

fig3b <- ggroc(roc.list, size=1.5*0.65, aes("color","linetype"))+
  theme(legend.position="none")+ theme_minimal() + ggtitle("Validation cohort") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid", size=0.2)+
  xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)") +
  labs(col="Scores", linetype="Scores")+ theme_pubclean()+
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size =17*0.65, face = 'plain'), axis.title.y = element_text(size=17*0.65, face = 'plain'),
        axis.text.x = element_text(size = 14*0.65, face ='plain'), axis.text.y = element_text(size=14*0.65, face ='plain'),
        legend.text = element_text(size = 17*0.65))+
  annotation_custom(tableGrob(table1, theme=ttheme_default(base_size = 9.75, padding = unit(c(2.5,2.5),"mm") )),
                    xmin=-0.5, xmax=0, ymin=0.2, ymax=0.4)+
  scale_color_manual(values=c("#994455","#6699CC","#EECC66","#364B9A","#F67B4E","#EE99AA","#C2E4EF"))+
  scale_linetype_manual(values=c("solid","11","dashed","longdash","dotted","twodash","solid"))+
  theme(legend.key.size = unit(12*0.65,"mm"))+
  theme(text=element_text(size=17*0.65),axis.text=element_text(size=15*0.65),legend.text=element_text(size=17*0.65),
        legend.position = "bottom", plot.title = element_text(size=18*0.65, face="bold", hjust=0.5))

f3.1 <- ggarrange(fig3a," ",fig3b, nrow=1, ncol = 3,labels = c("a","","b"), widths=c(10,1,10),
          font.label = list("size"=17.5*0.65), legend.grob = get_legend(fig3a), legend = "bottom") 


####---- ROC by sex, comorbidities, ethnicity and age group ----####
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


### ROC curves by comorbidities ###
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


### ROC curves by age group ###
test_nhanes <- test_nhanes %>% mutate(age_cat = Age %>% cut(c(-Inf,65,80,Inf))) %>%
  mutate(age_cat=factor(age_cat, labels = c("Middle-age", "Old", "Very Old")))

RC_Antropo <- OptimalCutpoints::optimal.cutpoints(X = AnthropoAge~mortstat, tag.healthy = 0, categorical.cov = "age_cat",
                                                  methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                                  conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()

RC_Antropo2 <- OptimalCutpoints::optimal.cutpoints(X = AnthropoAge2~mortstat, tag.healthy = 0, categorical.cov = "age_cat",
                                                   methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                                   conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()

RC_Pheno <- OptimalCutpoints::optimal.cutpoints(X = PhenoAge~mortstat, tag.healthy = 0, categorical.cov = "age_cat",
                                                methods = "Youden", data = test_nhanes, pop.prev = NULL, ci.fit=FALSE,
                                                conf.level = 0.95, control = OptimalCutpoints::control.cutpoints()) %>% summary()


Age1A <- RC_Antropo$p.table[["Middle-age"]][["AUC_CI"]]
Age1B <- RC_Antropo$p.table[["Old"]][["AUC_CI"]]
Age1C <- RC_Antropo$p.table[["Very Old"]][["AUC_CI"]]

Age2A <- RC_Antropo2$p.table[["Middle-age"]][["AUC_CI"]]
Age2B <- RC_Antropo2$p.table[["Old"]][["AUC_CI"]]
Age2C <- RC_Antropo2$p.table[["Very Old"]][["AUC_CI"]]

Age3A <- RC_Pheno$p.table[["Middle-age"]][["AUC_CI"]]
Age3B <- RC_Pheno$p.table[["Old"]][["AUC_CI"]]
Age3C <- RC_Pheno$p.table[["Very Old"]][["AUC_CI"]]


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
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2*0.65, size=1.25*0.65, position=position_dodge(.9))+
  geom_col(position = "dodge", color="black", size=1.25*0.65, alpha=0.6) + 
  coord_cartesian(ylim=c(0.75,0.90)) + theme_pubclean() + 
  scale_fill_manual(values=c("#994455","#6699CC","#EECC66"), labels=c("PhenoAge", "AnthropoAge", "S-AnthropoAge"))+
  scale_x_discrete(labels=c("Male","Female")) + ggtitle("AUROC by sex\n(validation cohort)")+theme(legend.key.size = unit(12*0.65,"mm"))+
  theme(text=element_text(size=17*0.65),axis.text=element_text(size=15*0.65),legend.text=element_text(size=17*0.65),
        legend.position = "bottom", plot.title = element_text(size=18*0.65, face="bold", hjust=0.5))


fig3d <- data.frame(rep(c("1A","2SA","3P"),3), c(rep("A0",3),rep("B1",3),rep("C>=2",3)),
                    c(ROC_Com3A$ci[2],ROC_Com1A$ci[2],ROC_Com2A$ci[2],ROC_Com3B$ci[2],ROC_Com1B$ci[2],ROC_Com2B$ci[2],
                      ROC_Com3C$ci[2],ROC_Com1C$ci[2],ROC_Com2C$ci[2]),
                    c(ROC_Com3A$ci[1],ROC_Com1A$ci[1],ROC_Com2A$ci[1],ROC_Com3B$ci[1],ROC_Com1B$ci[1],ROC_Com2B$ci[1],
                      ROC_Com3C$ci[1],ROC_Com1C$ci[1],ROC_Com2C$ci[1]),
                    c(ROC_Com3A$ci[3],ROC_Com1A$ci[3],ROC_Com2A$ci[3],ROC_Com3B$ci[3],ROC_Com1B$ci[3],ROC_Com2B$ci[3],
                      ROC_Com3C$ci[3],ROC_Com1C$ci[3],ROC_Com2C$ci[3])) %>% 
  `names<-`(c("Variable","Comorbitidies","AUROC","Low","Up")) %>% 
  ggplot(aes(x=Comorbitidies, y=AUROC, fill=Variable)) +
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2*0.65, size=1.25*0.65, position=position_dodge(.9))+
  geom_col(position = "dodge", color="black", size=1.25*0.65, alpha=0.6) + 
  coord_cartesian(ylim=c(0.75,0.90)) + theme_pubclean() + scale_x_discrete(labels=c("0","1", ">=2")) + 
  scale_fill_manual(values=c("#994455","#6699CC","#EECC66"), labels=c("PhenoAge", "AnthropoAge", "S-AnthropoAge"))+
  ggtitle("AUROC by number of comorbidities\n(validation cohort)")+theme(legend.key.size = unit(12*0.65,"mm"))+
  theme(text=element_text(size=17*0.65),axis.text=element_text(size=15*0.65),legend.text=element_text(size=17*0.65),
        legend.position = "bottom", plot.title = element_text(size=18*0.65, face="bold", hjust=0.5))
f3.2<-ggarrange(fig3c," ",fig3d, nrow=1, ncol = 3,labels = c("c","","d"), widths=c(10,1,10),
                font.label = list("size"=17.5*0.65), legend.grob = get_legend(fig3c), legend = "bottom")

f3 <- ggarrange(f3.1,"",f3.2, nrow=3, ncol=1, heights = c(12,0.3,7.7))
ggsave(file = "Figure3.jpg", f3, bg = "transparent", width = 30, height = 21, units=c("cm"),dpi = 300, limitsize = FALSE)
#Submission
ggsave(f3, file = "Submission/Figures/Figure_4.pdf", bg = "transparent",
       width = 30, height = 21, units=c("cm"), dpi = 600, limitsize = F)



##AUROC by age group -  AnthropoAge##
ROC_Age1A <- test_nhanes%>%filter(age_cat=="Middle-age") %>% pROC::roc(mortstat, AnthropoAge, ci=T)
ROC_Age1B <- test_nhanes%>%filter(age_cat=="Old") %>% pROC::roc(mortstat, AnthropoAge, ci=T)
ROC_Age1C <- test_nhanes%>%filter(age_cat=="Very Old") %>% pROC::roc(mortstat, AnthropoAge, ci=T)

ROC_Age2A <- test_nhanes%>%filter(age_cat=="Middle-age") %>% pROC::roc(mortstat, AnthropoAge2, ci=T)
ROC_Age2B <- test_nhanes%>%filter(age_cat=="Old") %>% pROC::roc(mortstat, AnthropoAge2, ci=T)
ROC_Age2C <- test_nhanes%>%filter(age_cat=="Very Old") %>% pROC::roc(mortstat, AnthropoAge2, ci=T)

ROC_Age3A <- test_nhanes%>%filter(age_cat=="Middle-age") %>% pROC::roc(mortstat, PhenoAge, ci=T)
ROC_Age3B <- test_nhanes%>%filter(age_cat=="Old") %>% pROC::roc(mortstat, PhenoAge, ci=T)
ROC_Age3C <- test_nhanes%>%filter(age_cat=="Very Old") %>% pROC::roc(mortstat, PhenoAge, ci=T)

fig3e <- data.frame(rep(c("1A","2SA","3P"),3), c(rep("AMiddle-age",3),rep("BOld",3),rep("CVery Old",3)),
                    c(ROC_Age3A$ci[2],ROC_Age1A$ci[2],ROC_Age2A$ci[2],ROC_Age3B$ci[2],ROC_Age1B$ci[2],ROC_Age2B$ci[2],
                      ROC_Age3C$ci[2],ROC_Age1C$ci[2],ROC_Age2C$ci[2]),
                    c(ROC_Age3A$ci[1],ROC_Age1A$ci[1],ROC_Age2A$ci[1],ROC_Age3B$ci[1],ROC_Age1B$ci[1],ROC_Age2B$ci[1],
                      ROC_Age3C$ci[1],ROC_Age1C$ci[1],ROC_Age2C$ci[1]),
                    c(ROC_Age3A$ci[3],ROC_Age1A$ci[3],ROC_Age2A$ci[3],ROC_Age3B$ci[3],ROC_Age1B$ci[3],ROC_Age2B$ci[3],
                      ROC_Age3C$ci[3],ROC_Age1C$ci[3],ROC_Age2C$ci[3])) %>% 
  `names<-`(c("Variable","AgeCat","AUROC","Low","Up")) %>% 
  ggplot(aes(x=AgeCat, y=AUROC, fill=Variable)) +
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2, size=1.25, position=position_dodge(.9))+
  geom_col(position = "dodge", color="black", size=1.25, alpha=0.6) + 
  coord_cartesian(ylim=c(0.5,0.9)) + theme_pubclean() + xlab("Age group") +
  scale_fill_manual(values=c("#994455","#6699CC","#EECC66"), labels=c("PhenoAge", "AnthropoAge", "S-AnthropoAge"))+
  scale_x_discrete(labels=c("Middle-age","Old", "Very Old")) + ggtitle("AUROC by age group\n(validation cohort)")+theme(legend.key.size = unit(12,"mm"))+
  theme(text=element_text(size=17),axis.text=element_text(size=15),legend.text=element_text(size=17), legend.position = "bottom",
        plot.title = element_text(size=18, face="bold", hjust=0.5))


##AUROC by age group -  AnthropoAgeAccel##
ROC_Age4A <- test_nhanes%>%filter(age_cat=="Middle-age") %>% pROC::roc(mortstat, AnthropoAgeAccel, ci=T)
ROC_Age4B <- test_nhanes%>%filter(age_cat=="Old") %>% pROC::roc(mortstat, AnthropoAgeAccel, ci=T)
ROC_Age4C <- test_nhanes%>%filter(age_cat=="Very Old") %>% pROC::roc(mortstat, AnthropoAgeAccel, ci=T)

ROC_Age5A <- test_nhanes%>%filter(age_cat=="Middle-age") %>% pROC::roc(mortstat, AnthropoAgeAccel2, ci=T)
ROC_Age5B <- test_nhanes%>%filter(age_cat=="Old") %>% pROC::roc(mortstat, AnthropoAgeAccel2, ci=T)
ROC_Age5C <- test_nhanes%>%filter(age_cat=="Very Old") %>% pROC::roc(mortstat, AnthropoAgeAccel2, ci=T)

ROC_Age6A <- test_nhanes%>%filter(age_cat=="Middle-age") %>% pROC::roc(mortstat, PhenoAgeAccel, ci=T)
ROC_Age6B <- test_nhanes%>%filter(age_cat=="Old") %>% pROC::roc(mortstat, PhenoAgeAccel, ci=T)
ROC_Age6C <- test_nhanes%>%filter(age_cat=="Very Old") %>% pROC::roc(mortstat, PhenoAgeAccel, ci=T)

fig3f <- data.frame(rep(c("1A","2SA","3P"),3), c(rep("AMiddle-age",3),rep("BOld",3),rep("CVery Old",3)),
                    c(ROC_Age6A$ci[2],ROC_Age4A$ci[2],ROC_Age5A$ci[2],ROC_Age6B$ci[2],ROC_Age4B$ci[2],ROC_Age5B$ci[2],
                      ROC_Age6C$ci[2],ROC_Age4C$ci[2],ROC_Age5C$ci[2]),
                    c(ROC_Age6A$ci[1],ROC_Age4A$ci[1],ROC_Age5A$ci[1],ROC_Age6B$ci[1],ROC_Age4B$ci[1],ROC_Age5B$ci[1],
                      ROC_Age6C$ci[1],ROC_Age4C$ci[1],ROC_Age5C$ci[1]),
                    c(ROC_Age6A$ci[3],ROC_Age4A$ci[3],ROC_Age5A$ci[3],ROC_Age6B$ci[3],ROC_Age4B$ci[3],ROC_Age5B$ci[3],
                      ROC_Age6C$ci[3],ROC_Age4C$ci[3],ROC_Age5C$ci[3])) %>% 
  `names<-`(c("Variable","AgeCat","AUROC","Low","Up")) %>% 
  ggplot(aes(x=AgeCat, y=AUROC, fill=Variable)) +
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2, size=1.25, position=position_dodge(.9))+
  geom_col(position = "dodge", color="black", size=1.25, alpha=0.6) + 
  coord_cartesian(ylim=c(0.5,0.9)) + theme_pubclean() + xlab("Age group") +
  scale_fill_manual(values=c("#994455","#6699CC","#EECC66"), labels=c("PhenoAgeAccel", "AnthropoAgeAccel", "S-AnthropoAgeAccel"))+
  scale_x_discrete(labels=c("Middle-age","Old", "Very Old")) + ggtitle("AUROC by age group\n(validation cohort)")+theme(legend.key.size = unit(12,"mm"))+
  theme(text=element_text(size=17),axis.text=element_text(size=15),legend.text=element_text(size=17), legend.position = "bottom",
        plot.title = element_text(size=18, face="bold", hjust=0.5))

table(test_nhanes$age_cat)


####---- Cause-specific mortality ----####
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

HR1.1<-hr95(m1_cv); Z1.1<-summary(m1_cv)$coefficients[1,4]
HR1.2<-hr95(m2_cv); Z1.2<-summary(m2_cv)$coefficients[1,4]
HR1.3<-hr95(m3_cv); Z1.3<-summary(m3_cv)$coefficients[1,4]


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

HR2.1<-hr95(m1_db); Z2.1<-summary(m1_db)$coefficients[1,4]
HR2.2<-hr95(m2_db); Z2.2<-summary(m2_db)$coefficients[1,4]
HR2.3<-hr95(m3_db); Z2.3<-summary(m3_db)$coefficients[1,4]


### Stroke mortality ###
nhanes0$cer_mort <- NULL
nhanes0$cer_mort[nhanes0$Cerebrovascular_diseases==1 & nhanes0$mortstat==1]<-1
nhanes0$cer_mort[nhanes0$Cerebrovascular_diseases==0 & nhanes0$mortstat==1]<-2
nhanes0$cer_mort[nhanes0$Cerebrovascular_diseases==0 & nhanes0$mortstat==0]<-0
nhanes0$cer_mort<-factor(nhanes0$cer_mort, labels = c("Censored", "Stroke mortality", "Other causes"))
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

HR3.1<-hr95(m1_cer); Z3.1<-summary(m1_cer)$coefficients[1,4]
HR3.2<-hr95(m2_cer); Z3.2<-summary(m2_cer)$coefficients[1,4]
HR3.3<-hr95(m3_cer); Z3.3<-summary(m3_cer)$coefficients[1,4]


### Cancer-related mortality ###
nhanes0$can_mort <- NULL
nhanes0$can_mort[nhanes0$Malignant_neoplasms==1 & nhanes0$mortstat==1]<-1
nhanes0$can_mort[nhanes0$Malignant_neoplasms==0 & nhanes0$mortstat==1]<-2
nhanes0$can_mort[nhanes0$Malignant_neoplasms==0 & nhanes0$mortstat==0]<-0
nhanes0$can_mort<-factor(nhanes0$can_mort, labels = c("Censored", "Cancer mortality", "Other causes"))
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

HR4.1<-hr95(m1_can); Z4.1<-summary(m1_can)$coefficients[1,4]
HR4.2<-hr95(m2_can); Z4.2<-summary(m2_can)$coefficients[1,4]
HR4.3<-hr95(m3_can); Z4.3<-summary(m3_can)$coefficients[1,4]


### Influenza/Pneumonia mortality ###
nhanes0$inf <- NULL
nhanes0$inf[nhanes0$Influenza_or_pneumonia==1 & nhanes0$mortstat==1]<-1
nhanes0$inf[nhanes0$Influenza_or_pneumonia==0 & nhanes0$mortstat==1]<-2
nhanes0$inf[nhanes0$Influenza_or_pneumonia==0 & nhanes0$mortstat==0]<-0
nhanes0$inf<-factor(nhanes0$inf, labels = c("Censored", "Influenza mortality", "Other causes"))
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

HR5.1<-hr95(m1_inf); Z5.1<-summary(m1_inf)$coefficients[1,4]
HR5.2<-hr95(m2_inf); Z5.2<-summary(m2_inf)$coefficients[1,4]
HR5.3<-hr95(m3_inf); Z5.3<-summary(m3_inf)$coefficients[1,4]


### CKD mortality ###
nhanes0$ckd <- NULL
nhanes0$ckd[nhanes0$Nephone_diseases==1 & nhanes0$mortstat==1]<-1
nhanes0$ckd[nhanes0$Nephone_diseases==0 & nhanes0$mortstat==1]<-2
nhanes0$ckd[nhanes0$Nephone_diseases==0 & nhanes0$mortstat==0]<-0
nhanes0$ckd<-factor(nhanes0$ckd, labels = c("Censored", "CKD mortality", "Other causes"))
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

HR6.1<-hr95(m1_ckd); Z6.1<-summary(m1_ckd)$coefficients[1,4]
HR6.2<-hr95(m2_ckd); Z6.2<-summary(m2_ckd)$coefficients[1,4]
HR6.3<-hr95(m3_ckd); Z6.3<-summary(m3_ckd)$coefficients[1,4]


### Alzheimer mortality ###
nhanes0$alz <- NULL
nhanes0$alz[nhanes0$Alzheimer_disease==1 & nhanes0$mortstat==1]<-1
nhanes0$alz[nhanes0$Alzheimer_disease==0 & nhanes0$mortstat==1]<-2
nhanes0$alz[nhanes0$Alzheimer_disease==0 & nhanes0$mortstat==0]<-0
nhanes0$alz<-factor(nhanes0$alz, labels = c("Censored", "AD  mortality", "Other causes"))
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

HR7.1<-hr95(m1_alz); Z7.1<-summary(m1_alz)$coefficients[1,4]
HR7.2<-hr95(m2_alz); Z7.2<-summary(m2_alz)$coefficients[1,4]
HR7.3<-hr95(m3_alz); Z7.3<-summary(m3_alz)$coefficients[1,4]


### COPD mortality ###
nhanes0$copd <- NULL
nhanes0$copd[nhanes0$Chronic_lower_respiratory_diseases==1 & nhanes0$mortstat==1]<-1
nhanes0$copd[nhanes0$Chronic_lower_respiratory_diseases==0 & nhanes0$mortstat==1]<-2
nhanes0$copd[nhanes0$Chronic_lower_respiratory_diseases==0 & nhanes0$mortstat==0]<-0
nhanes0$copd<-factor(nhanes0$copd, labels = c("Censored", "COPD mortality", "Other causes"))
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

HR8.1<-hr95(m1_copd); Z8.1<-summary(m1_copd)$coefficients[1,4]
HR8.2<-hr95(m2_copd); Z8.2<-summary(m2_copd)$coefficients[1,4]
HR8.3<-hr95(m3_copd); Z8.3<-summary(m3_copd)$coefficients[1,4]

### All-cause mortality (Gompertz) ###
m1_all<-flexsurvreg(Surv(permth_int, mortstat)~ AnthropoAge+Sex+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")
m2_all<-flexsurvreg(Surv(permth_int, mortstat)~ AnthropoAge2+Sex+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")
m3_all<-flexsurvreg(Surv(permth_int, mortstat)~ PhenoAge+Sex+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")

HR9.1 <- (m1_all)$res[3,1:3] %>% exp; Z9.1 <- m1_all$res.t[3,4]
HR9.2 <- (m2_all)$res[3,1:3] %>% exp; Z9.2 <- m2_all$res.t[3,4]
HR9.3 <- (m3_all)$res[3,1:3] %>% exp; Z9.3 <- m3_all$res.t[3,4]

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
causes<-c("Cardiovascular", "Diabetes Mellitus", "Stroke", "Cancer", "Influenza/Pneumonia", "Nephritis/Nephrosis", "Alzheimer", "Chronic lower respiratory")

tab3<-data.frame("Cause-specific mortality"=c(causes),
                 "AnthropoAge \n c-statistic (95%CI)"=c(paste0(antro," (",antro_ci,")")),
                 "S-AnthropoAge \n c-statistic (95%CI)"=c(paste0(antro_s," (",antros_ci,")")),
                 "PhenoAge \n c-statistic (95%CI)"=c(paste0(pheno," (",pheno_ci,")")),
                 "AnthropoAge vs. \n S-AnthropoAge"=c(antro_b),
                 "AnthropoAge vs. \n PhenoAge"=c(antro_s_b),
                 "S-AnthropoAge vs. \n PhenoAge"=c(pheno_b))

tab3<-`names<-`(tab3,c("Cause-specific mortality",
                       "AnthropoAge \n c-statistic (95%CI)",
                       "S-AnthropoAge \n c-statistic (95%CI)",
                       "PhenoAge \n c-statistic (95%CI)",
                       "AnthropoAge vs. \n S-AnthropoAge",
                       "AnthropoAge vs. \n PhenoAge",
                       "S-AnthropoAge vs. \n PhenoAge"))

tab3<-align(flextable(tab3),align = "center",part = "all") %>% autofit()
doc <- read_docx() %>% body_add_flextable(value = tab3, split = TRUE) %>%
  body_end_section_landscape() %>% print( target = "tabla3.docx" )


### SUPPLEMENTARY TABLE 5 ###
HR_A1 <- (as.list(paste0("HR",1:9,".1")) %>% lapply(get)) %>% lapply(ci_quick) %>% as.character()
HR_A2 <- (as.list(paste0("HR",1:9,".2")) %>% lapply(get)) %>% lapply(ci_quick) %>% as.character()
HR_PH <- (as.list(paste0("HR",1:9,".3")) %>% lapply(get)) %>% lapply(ci_quick) %>% as.character()
Z_A1 <- paste0("SE=",(as.list(paste0("Z",1:9,".1")) %>% lapply(get)) %>% lapply(round,4) %>% as.character())
Z_A2 <- paste0("SE=",(as.list(paste0("Z",1:9,".2")) %>% lapply(get)) %>% lapply(round,4) %>% as.character())
Z_PH <- paste0("SE=",(as.list(paste0("Z",1:9,".3")) %>% lapply(get)) %>% lapply(round,4) %>% as.character())

supptab5 <- data.frame(c(causes, "All-cause mortality"),paste0(HR_A1,",\n",Z_A1),paste0(HR_A2,",\n",Z_A2),
                       paste0(HR_PH,",\n",Z_PH)) %>% `names<-`(c("Mortality","AnthropoAge","S-AnthropoAge","PhenoAge")) %>% 
  flextable() %>% align(align = "center",part = "all") %>% autofit()

doc <- read_docx() %>% body_add_flextable(value = supptab5, split = TRUE) %>%
  body_end_section_landscape() %>% print( target = "tabla_supp5.docx" )



##BA-Sex interaction
mort_c <- list(fgdata_cv, fgdata_db, fgdata_cer, fgdata_can, fgdata_inf, fgdata_ckd, fgdata_alz, fgdata_copd)
mrep <- function(x){coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge*Sex+Ethnicity+num_comorb, weight=fgwt, data=x)}

c(exp(coef(mrep(mort_c[[1]]))[7]),exp(confint(mrep(mort_c[[1]]))[7,])) %>% round(3) #CVD #*
c(exp(coef(mrep(mort_c[[2]]))[7]),exp(confint(mrep(mort_c[[2]]))[7,])) %>% round(3) #DM2 #*
c(exp(coef(mrep(mort_c[[3]]))[7]),exp(confint(mrep(mort_c[[3]]))[7,])) %>% round(3) #STR #*
c(exp(coef(mrep(mort_c[[4]]))[7]),exp(confint(mrep(mort_c[[4]]))[7,])) %>% round(3) #CAN #*
c(exp(coef(mrep(mort_c[[5]]))[7]),exp(confint(mrep(mort_c[[5]]))[7,])) %>% round(3) #INF
c(exp(coef(mrep(mort_c[[6]]))[7]),exp(confint(mrep(mort_c[[6]]))[7,])) %>% round(3) #CKD
c(exp(coef(mrep(mort_c[[7]]))[7]),exp(confint(mrep(mort_c[[7]]))[7,])) %>% round(3) #ALZ
c(exp(coef(mrep(mort_c[[8]]))[7]),exp(confint(mrep(mort_c[[8]]))[7,])) %>% round(3) #COP

#Cardiovascular mortality
#m1_cv2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cv)
#m2_cv2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cv)
#m3_cv2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cv)
#c(exp(coef(m1_cv2)[7]),exp(confint(m1_cv2)[7,])) %>% round(3); c(exp(coef(m2_cv2)[7]),exp(confint(m2_cv2)[7,])) %>% round(3); c(exp(coef(m3_cv2)[7]),exp(confint(m3_cv2)[7,])) %>% round(3)

#Diabetes mortality
#m1_dm2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_db)
#m2_dm2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_db)
#m3_dm2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_db)
#c(exp(coef(m1_dm2)[7]),exp(confint(m1_dm2)[7,])) %>% round(3); c(exp(coef(m2_dm2)[7]),exp(confint(m2_dm2)[7,])) %>% round(3); c(exp(coef(m3_dm2)[7]),exp(confint(m3_dm2)[7,])) %>% round(3)

#Stroke mortality
#m1_ce2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cer)
#m2_ce2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cer)
#m3_ce2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_cer)
#c(exp(coef(m1_ce2)[7]),exp(confint(m1_ce2)[7,])) %>% round(3); c(exp(coef(m2_ce2)[7]),exp(confint(m2_ce2)[7,])) %>% round(3); c(exp(coef(m3_ce2)[7]),exp(confint(m3_ce2)[7,])) %>% round(3)

#Cancer mortality
#m1_ca2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_can)
#m2_ca2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ AnthropoAge2*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_can)
#m3_ca2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PhenoAge*Sex+Ethnicity+num_comorb,weight=fgwt, data=fgdata_can)
#c(exp(coef(m1_ca2)[7]),exp(confint(m1_ca2)[7,])) %>% round(3); c(exp(coef(m2_ca2)[7]),exp(confint(m2_ca2)[7,])) %>% round(3); c(exp(coef(m3_ca2)[7]),exp(confint(m3_ca2)[7,])) %>% round(3)



####---- Multidimensional aging ----####
nhanes0$accel1<-ifelse(nhanes0$AnthropoAgeAccel>0, 1, 0)
nhanes0$accel2<-ifelse(nhanes0$AnthropoAgeAccel2>0, 1, 0)
nhanes0$accel3<-ifelse(nhanes0$PhenoAgeAccel>0, 1, 0)

### Kaplan-Meier
colors_border<-c("#6699CC","#994455")
colors_in<- rgb(col2rgb(colors_border)[1,],col2rgb(colors_border)[2,],
                col2rgb(colors_border)[3,],max=255,alpha=(50)*(255/100))

mod2_kma<-survfit(Surv(permth_int, mortstat) ~ accel1, data = nhanes0)
f3a<-ggsurvplot(mod2_kma, fontsize = 4.5*0.75, data = nhanes0, palette = colors_border,conf.int = T, risk.table = T,pval = TRUE, xlab="Time (Months)",
                ylab="Survival probability", title="AnthropoAgeAccel", pval.coord = c(0, 0.8), pval.size=5*0.7, legend.labs = c("Phys-aging","Accel-aging"),
                ylim= c(0.7,1.0), xlim=c(0, 120), break.y.by= c(0.1), break.x.by= c(40), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15*0.8),
                         text = element_text(hjust=0.5, family = "sans", size=15*0.75),
                         legend.text = element_text(hjust=0.5, size=11))))
f3a$table <- f3a$table + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5, color="black", margin = margin(r=20*0.7), face="bold", size=12*0.75, angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size=15*0.75)) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
f3a$plot <- f3a$plot + labs(fill=NULL, color=NULL) + theme(legend.justification = c(0.5,0.5), legend.key.size = unit(5,"mm"))
fig3a<-ggarrange(f3a$plot, f3a$table, heights = c(2, 0.7), ncol = 1, nrow = 2)

mod2_kmb<-survfit(Surv(permth_int, mortstat) ~ accel2, data = nhanes0)
f3b<-ggsurvplot(mod2_kmb, fontsize = 4.5*0.75, data = nhanes0, palette = colors_border,conf.int = T, risk.table = T,pval = TRUE, xlab="Time (Months)",
                ylab="Survival probability", title="S-AnthropoAgeAccel", pval.coord = c(0, 0.8), pval.size=5*0.7, legend.labs = c("Phys-aging","Accel-aging"),
                ylim= c(0.7,1.0), xlim=c(0, 120), break.y.by= c(0.1), break.x.by= c(40), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15*0.8),
                         text = element_text(hjust=0.5, family = "sans", size=15*0.75),
                         legend.text = element_text(hjust=0.5, size=11))))
f3b$table <- f3b$table + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5, color="black", margin = margin(r=20*0.7), face="bold", size=12*0.75, angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size=15*0.75)) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
f3b$plot <- f3b$plot + labs(fill=NULL, color=NULL) + theme(legend.justification = c(0.5,0.5), legend.key.size = unit(5,"mm"))
fig3b<-ggarrange(f3b$plot, f3b$table, heights = c(2, 0.7), ncol = 1, nrow = 2)


mod2_kmc<-survfit(Surv(permth_int, mortstat) ~ accel3, data = nhanes0)
f3c<-ggsurvplot(mod2_kmc, fontsize = 4.5*0.75, data = nhanes0, palette = colors_border,conf.int = T, risk.table = T,pval = TRUE, xlab="Time (Months)",
                ylab="Survival probability", title="PhenoAgeAccel", pval.coord = c(0, 0.8), pval.size=5*0.7, legend.labs = c("Phys-aging","Accel-aging"),
                ylim= c(0.7,1.0), xlim=c(0, 120), break.y.by= c(0.1), break.x.by= c(40), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15*0.8),
                         text = element_text(hjust=0.5, family = "sans", size=15*0.75),
                         legend.text = element_text(hjust=0.5, size=11))))
f3c$table <- f3c$table + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5, color="black", margin = margin(r=20*0.7), face="bold", size=12*0.75, angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size=15*0.75)) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
f3c$plot <- f3c$plot + labs(fill=NULL, color=NULL) + theme(legend.justification = c(0.5,0.5), legend.key.size = unit(5,"mm"))
fig3c<-ggarrange(f3c$plot, f3c$table, heights = c(2, 0.7), ncol = 1, nrow = 2)


#### Combined acceleration mortality ###
## Categorical interaction graph
nhanes0$accel_comb<-factor( 2*(nhanes0$PhenoAgeAccel>0) + (nhanes0$AnthropoAgeAccel>0)); table(nhanes0$accel_comb)
mod2_kmd<-survfit(Surv(permth_int, mortstat) ~ factor(accel_comb), data = nhanes0)

colors_border2<-c("#80C4E7","#C29300","#AA4D10","#911F2C")
f3d<-ggsurvplot(mod2_kmd, data = nhanes0, fontsize = 4.5*0.75, palette = colors_border2,conf.int = T, risk.table = T,pval = TRUE, xlab="Time (Months)",
                ylab="Survival probability", title="Multidomain aging", pval.coord = c(0, 0.8), pval.size=5*0.7,
                legend.labs = c("Phys-aging", "Accel-Anthropo", "Accel-Pheno", "Multi-accel"),
                ylim= c(0.7,1.0), xlim=c(0, 120), break.y.by= c(0.1), break.x.by= c(40), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15*0.8),
                         text = element_text(hjust=0.5, family = "sans", size=15*0.75),
                         legend.text = element_text(hjust=0.5, size=11))))
f3d$table <- f3d$table + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5, color="black", margin = margin(r=20*0.7), face="bold", size=12*0.75, angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size=15*0.75)) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
f3d$plot <- f3d$plot + labs(fill=NULL, color=NULL) + theme(legend.justification = c(0.5,0.5), legend.key.size = unit(5,"mm")) +
  scale_color_manual(values=colors_border2, labels=c("Phys-aging", "Accel\nAnthropoAge", "Accel\nPhenoAge", "Multi-accel")) +
  scale_fill_manual(values=colors_border2, labels=c("Phys-aging", "Accel\nAnthropoAge", "Accel\nPhenoAge", "Multi-accel"))
fig3d<-ggarrange(f3d$plot, f3d$table, heights = c(2, 0.7), ncol = 1, nrow = 2)

fig5<-ggarrange(fig3a, fig3b, "", "", fig3c, fig3d, labels=c("a","b","","","c","d"), heights=c(10,1,10), ncol=2, nrow=3)
ggsave(file="Figure5.jpg", fig5, bg="transparent", width=28, height=18.5, units=c("cm"), dpi=500, limitsize = FALSE)
#Submission
ggsave(fig5, file="Submission/Figures/Figure_5.pdf", bg="transparent",
       width=28, height=18.5, units=c("cm"), dpi=600, limitsize = FALSE)
#Graphical abstract
f3d$plot + theme_void() + theme(plot.title=element_blank(), legend.position = "none", text = element_blank())



19*(29/18.5)

#Categorical interaction AnthropoAgeAccel
m_fin<-flexsurvreg(Surv(permth_int, mortstat)~accel_comb+Sex+Age+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")
round(cbind("coef"=exp(coef(m_fin))[3:8], exp(confint(m_fin))[3:8,]),2)
#Categorical interaction S-AnthropoAgeAccel
nhanes0$accel_comb2<-factor(2*(nhanes0$PhenoAgeAccel>0)+(nhanes0$AnthropoAgeAccel2>0)); table(nhanes0$accel_comb2)
m_fin2<-flexsurvreg(Surv(permth_int, mortstat)~accel_comb2+Sex+Age+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")
round(cbind("coef"=exp(coef(m_fin2))[3:8], exp(confint(m_fin2))[3:8,]),2)
#Adjustment
m_fin3<-flexsurvreg(Surv(permth_int, mortstat)~AnthropoAgeAccel+PhenoAgeAccel+Sex+Age+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")
round(cbind("coef"=exp(coef(m_fin3))[3:7], exp(confint(m_fin3))[3:7,]),2)
#CONTINUOUS interaction AnthropoAgeAccel
m_fin4<-flexsurvreg(Surv(permth_int, mortstat)~AnthropoAgeAccel*PhenoAgeAccel+Age+Sex+num_comorb+shape(Ethnicity), data=nhanes0, dist = "gompertz")
round(cbind("coef"=exp(coef(m_fin4))[3:8], exp(confint(m_fin4))[3:8,]),5)


####---- Comorbidity Assessment ----####
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
remove(nhanes0d)
supp5<-ggarrange(s5a, s5b, labels=letters[1:2])
ggsave(file = "SuppFig5.jpg", supp5, bg = "transparent", width = 37,height = 20, units=c("cm"), dpi = 500, limitsize = FALSE)


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
ggsave(file = "SuppFig6.jpg", supp6, bg = "transparent", width = 37,height = 20, units=c("cm"), dpi = 500, limitsize = FALSE)



####---- Spiderplots Databases ----####
index_del <- nhanes$SEQN[duplicated(nhanes$SEQN)]

### Data management: Anthropometry and DXA ###
dxa<-NHANES1 %>% filter(!is.infinite(PhenoAge)&!is.na(PhenoAge)) %>% filter() %>% 
  dplyr::select(SEQN,Age,mortstat,permth_int,Sex,Ethnicity,SEQN, PhenoAge, BMI,ICE,Weight,Height,Waist,
                Leg_length,Arm_length,Triceps_skinfold,Subscapular_skinfold,Thigh_circumference,Arm_circumference,
                DXXHEFAT,DXDHELE,DXXLAFAT,DXDLALE,DXXLLFAT,DXDLLLE,DXXRAFAT,DXDRALE,DXXRLFAT,DXDRLLE,DXXLSBMD,
                DXXTRFAT,DXDTRLE,DXDTOFAT,DXDTOLE,DXDTOBMD) %>% filter(!duplicated(SEQN)) %>% na.omit() %>% 
  mutate("Ethnicity"=factor(Ethnicity, labels = c("White", "Black", "Mexican-American", "Other")),
         "tr_height"=(Height)**(1/3),
         "tr_leg"=Leg_length,
         "tr_arm"=log(Arm_length),
         "tr_imc"=log(BMI),
         "tr_weight"=log(Weight),
         "tr_ice"=(ICE)**(1/3),
         "tr_tric"=(Triceps_skinfold)**(1/3),
         "tr_subs"=(Subscapular_skinfold)**(1/3),
         "tr_armc"=sqrt(Arm_circumference),
         "tr_thigh"=log(Thigh_circumference),
         "arm_fmi"=(DXXLAFAT/1000+DXXRAFAT/1000)/((Height/100)^2), #fmi=fat mass index
         "arm_lmi"=(DXDLALE/1000+DXDRALE/1000)/((Height/100)^2), #lmi=lean mass index
         "leg_fmi"=(DXXLLFAT/1000+DXXRLFAT/1000)/((Height/100)^2),
         "leg_lmi"=(DXDLLLE/1000+DXDRLLE/1000)/((Height/100)^2),
         "afmi"=(DXXLAFAT/1000+DXXRAFAT/1000+DXXLLFAT/1000+DXXRLFAT/1000)/((Height/100)^2), #
         "almi"=(DXDLALE/1000+DXDRALE/1000+DXDLLLE/1000+DXDRLLE/1000)/((Height/100)^2), #
         "head_fmi"=(DXXHEFAT/1000)/((Height/100)^2),
         "head_lmi"=(DXDHELE/1000)/((Height/100)^2),
         "trunk_fmi"=(DXXTRFAT/1000)/((Height/100)^2), #
         "trunk_lmi"=(DXDTRLE/1000)/((Height/100)^2),
         "total_fmi"=(DXDTOFAT/1000)/((Height/100)^2), #
         "total_lmi"=(DXDTOLE/1000)/((Height/100)^2), #
         "total_fap"=(DXDTOFAT/1000)/Weight,
         "total_lep"=(DXDTOLE/1000)/Weight,
         "app_ratio"=(DXXLAFAT+DXXRAFAT+DXXLLFAT+DXXRLFAT)/(DXDLALE+DXDRALE+DXDLLLE+DXDRLLE),
         "trunk_ratio"=DXXTRFAT/DXDTRLE, #
         "total_ratio"=DXDTOFAT/DXDTOLE, #
         "lumbar_bmd"=DXXLSBMD, #
         "total_bmd"=DXDTOBMD, #
         "ttafm_ratio"=(DXXTRFAT)/(DXXLAFAT+DXXRAFAT+DXXLLFAT+DXXRLFAT) #Trunk-to-appendicular fat mass ratio
         ); nrow(dxa)
rownames(dxa) <- dxa$SEQN; dxa <- dxa[!(row.names(dxa) %in% index_del),]; nrow(dxa)

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

with(dxa, table("Pheno"=PhenoAgeAccel>0, "Anthro"=AnthropoAgeAccel>0)) %>% print() %>% mcnemar.test()
with(dxa, table("S-A"=AnthropoAgeAccel2>0, "A"=AnthropoAgeAccel>0)) %>% print() %>% mcnemar.test()


### Data management: PhenoAge components ###
pheno<-NHANES1%>%filter(!is.infinite(PhenoAge)&!is.na(PhenoAge))%>%
  dplyr::select(SEQN,Age,PhenoAge,mortstat,permth_int,Sex,Ethnicity,BMI,Thigh_circumference,Arm_circumference,Waist,
                Triceps_skinfold,Subscapular_skinfold,Leg_length,Arm_length,Height,Weight,ICE,
                Glucose,Albumin,RDW,MCV,LymphP,WBC,Creatinine,CRP,ALP) %>% filter(!duplicated(SEQN)) %>% drop_na() %>%
  mutate("Ethnicity"=factor(Ethnicity, labels = c("White", "Black", "Mexican-American", "Other")),
         "tr_weight"=log(Weight),
         "tr_armc"=sqrt(Arm_circumference),
         "tr_ice"=(ICE)**(1/3),
         "tr_imc"=log(BMI),
         "tr_tric"=(Triceps_skinfold)**(1/3),
         "tr_subs"=(Subscapular_skinfold)**(1/3),
         "tr_arm"=log(Arm_length),
         "tr_thigh"=log(Thigh_circumference)); nrow(pheno)
rownames(pheno) <- pheno$SEQN; pheno <- pheno[!(row.names(pheno) %in% index_del),]; nrow(pheno)

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

m1<-lm(PhenoAge~Age, data=pheno)
pheno$PhenoAgeAccel<-m1$residuals

####---- AnthropoAgeAccel>0 Spiderplots ----####
### Spiderplot: AnthropoAge components ###
#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Age~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a1
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Leg_length~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a2
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Weight~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a3
(with(dxa%>%filter(Sex=="Women"), wilcox.test(BMI~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a4
(with(dxa%>%filter(Sex=="Women"), wilcox.test(ICE~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a5
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Subscapular_skinfold~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a6
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Triceps_skinfold~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a7
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Arm_circumference~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a8
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Thigh_circumference~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a9
(with(dxa%>%filter(Sex=="Women"), wilcox.test(AnthropoAge~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("CA","Leg\nlength", "Weight", "BMI", "WHtR","Subscapular\nskinfold",
                  "Triceps\nskinfold", "Arm\ncircumference","Thigh\ncircumference","AnthropoAge"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data1 <- ((accel1 <- dxa %>% filter(Sex=="Women"&Age) %>% 
             dplyr::select(Age,Leg_length,Weight,BMI,ICE,Subscapular_skinfold,Triceps_skinfold,
                           Arm_circumference,Thigh_circumference,AnthropoAge,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6a<-ggplotify::as.ggplot(~radarchart(data1, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in, plwd=2.25, pty = ".",
                                      plty=1, title="", cex.main=1.5*0.4, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8*0.4, vlcex=0.9))+
  theme(plot.margin=unit(c(-1.5*0.4,-0.75*0.4,-2.5*0.4,-2.25*0.4),"cm"))





#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Age~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a1
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Leg_length~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a2
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Weight~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a3
(with(dxa%>%filter(Sex=="Men"), wilcox.test(BMI~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a4
(with(dxa%>%filter(Sex=="Men"), wilcox.test(ICE~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a5
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Subscapular_skinfold~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a6
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Triceps_skinfold~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a7
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Arm_circumference~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a8
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Thigh_circumference~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a9
(with(dxa%>%filter(Sex=="Men"), wilcox.test(AnthropoAge~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("CA","Leg\nlength", "Weight", "BMI", "WHtR","Subscapular\nskinfold",
                  "Triceps\nskinfold", "Arm\ncircumference","Thigh\ncircumference","AnthropoAge"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data2 <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
             dplyr::select(Age,Leg_length,Weight,BMI,ICE,Subscapular_skinfold,Triceps_skinfold,
                           Arm_circumference,Thigh_circumference,AnthropoAge,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6d<-ggplotify::as.ggplot(~radarchart(data2, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in, plwd=2.25, pty = ".",
                                      plty=1, title="", cex.main=1.5*0.4, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8*0.4, vlcex=0.9))+
  theme(plot.margin=unit(c(-1.5*0.4,-0.75*0.4,-2.5*0.4,-2.25*0.4),"cm"))

### Spiderplot: DXA ###
#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_fmi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a1
(with(dxa%>%filter(Sex=="Women"), wilcox.test(trunk_fmi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a2
(with(dxa%>%filter(Sex=="Women"), wilcox.test(afmi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a3
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_lmi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a4
(with(dxa%>%filter(Sex=="Women"), wilcox.test(almi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a5
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_ratio~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a6
(with(dxa%>%filter(Sex=="Women"), wilcox.test(trunk_ratio~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a7
(with(dxa%>%filter(Sex=="Women"), wilcox.test(ttafm_ratio~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a8
(with(dxa%>%filter(Sex=="Women"), wilcox.test(lumbar_bmd~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a9
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_bmd~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("Body\nFMI", "Trunk\nFMI", "AFMI", "Body\nLMI", "ALMI", "Body fat-lean ratio",
                  "Trunk fat-lean\nratio", "Trunk-Append\nfat ratio", "Lumbar\nBMD","Total BMD"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data3 <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
             dplyr::select(total_fmi, trunk_fmi, afmi, total_lmi, almi, total_ratio,
                           trunk_ratio, ttafm_ratio, lumbar_bmd, total_bmd,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6b<-ggplotify::as.ggplot(~radarchart(data3, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in, plwd=2.25, pty = ".",
                                       plty=1, title="", cex.main=1.5*0.4, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8*0.4, vlcex=0.9))+
  theme(plot.margin=unit(c(-1.5*0.4,-0.75*0.4,-2.5*0.4,-2.25*0.4),"cm"))

#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_fmi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a1
(with(dxa%>%filter(Sex=="Men"), wilcox.test(trunk_fmi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a2
(with(dxa%>%filter(Sex=="Men"), wilcox.test(afmi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a3
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_lmi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a4
(with(dxa%>%filter(Sex=="Men"), wilcox.test(almi~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a5
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_ratio~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a6
(with(dxa%>%filter(Sex=="Men"), wilcox.test(trunk_ratio~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a7
(with(dxa%>%filter(Sex=="Men"), wilcox.test(ttafm_ratio~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a8
(with(dxa%>%filter(Sex=="Men"), wilcox.test(lumbar_bmd~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a9
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_bmd~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("Body\nFMI", "Trunk\nFMI", "AFMI", "Body\nLMI", "ALMI", "Body fat-lean ratio",
                  "Trunk fat-lean\nratio", "Trunk-Append\nfat ratio", "Lumbar\nBMD","Total BMD"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data4 <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
             dplyr::select(total_fmi, trunk_fmi, afmi, total_lmi, almi, total_ratio,
                           trunk_ratio, ttafm_ratio, lumbar_bmd, total_bmd, AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6e<-ggplotify::as.ggplot(~radarchart(data4, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in, plwd=2.25, pty = ".",
                                       plty=1, title="", cex.main=1.5*0.4, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8*0.4, vlcex=0.9))+
  theme(plot.margin=unit(c(-1.5*0.4,-0.75*0.4,-2.5*0.4,-2.25*0.4),"cm"))


### Spiderplots: PhenoAge components ###
#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(pheno%>%filter(Sex=="Women"), wilcox.test(PhenoAge~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a1
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Glucose~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a2
(with(pheno%>%filter(Sex=="Women"), t.test(Creatinine~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a3
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Albumin~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a4
(with(pheno%>%filter(Sex=="Women"), wilcox.test(ALP~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a5
(with(pheno%>%filter(Sex=="Women"), wilcox.test(WBC~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a6
(with(pheno%>%filter(Sex=="Women"), wilcox.test(LymphP~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a7
(with(pheno%>%filter(Sex=="Women"), wilcox.test(MCV~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a8
(with(pheno%>%filter(Sex=="Women"), wilcox.test(RDW~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a9
(with(pheno%>%filter(Sex=="Women"), wilcox.test(CRP~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("PhenoAge","Glucose","Creatinine","Albumin","ALP","WBC",
                  "Lymphocyte\nPercentage","MCV","RDW","CRP"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data5 <- ((accel1 <- pheno %>% filter(Sex=="Women") %>% 
             dplyr::select(PhenoAge,Glucose,Creatinine,Albumin,ALP,WBC,LymphP,MCV,RDW,CRP,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6c<-ggplotify::as.ggplot(~radarchart(data5, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in, plwd=2.25, pty = ".",
                                      plty=1, title="", cex.main=1.5*0.4, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8*0.4, vlcex=0.9))+
  theme(plot.margin=unit(c(-1.5*0.4,-0.75*0.4,-2.5*0.4,-2.25*0.4),"cm"))


#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(pheno%>%filter(Sex=="Men"), wilcox.test(PhenoAge~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a1
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Glucose~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a2
(with(pheno%>%filter(Sex=="Men"), t.test(Creatinine~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a3
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Albumin~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a4
(with(pheno%>%filter(Sex=="Men"), wilcox.test(ALP~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a5
(with(pheno%>%filter(Sex=="Men"), wilcox.test(WBC~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a6
(with(pheno%>%filter(Sex=="Men"), wilcox.test(LymphP~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a7
(with(pheno%>%filter(Sex=="Men"), wilcox.test(MCV~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a8
(with(pheno%>%filter(Sex=="Men"), wilcox.test(RDW~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a9
(with(pheno%>%filter(Sex=="Men"), wilcox.test(CRP~AnthropoAgeAccel>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("PhenoAge","Glucose","Creatinine","Albumin","ALP","WBC",
                  "Lymphocyte\nPercentage","MCV","RDW","CRP"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data6 <- ((accel1 <- pheno %>% filter(Sex=="Men") %>% 
             dplyr::select(PhenoAge,Glucose,Creatinine,Albumin,ALP,WBC,LymphP,MCV,RDW,CRP,AnthropoAgeAccel) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

f6f<-ggplotify::as.ggplot(~radarchart(data6, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in, plwd=2.25, pty = ".",
                                      plty=1, title="", cex.main=1.5*0.4, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8*0.4, vlcex=0.9))+
  theme(plot.margin=unit(c(-1.5*0.4,-0.75*0.4,-2.5*0.4,-2.25*0.4),"cm"))



### Joint ###
dxa$accel_comb<-factor((dxa$PhenoAgeAccel>0)+2*(dxa$AnthropoAgeAccel>0))
g1<-nhanes0 %>% mutate(accel1=factor(accel1, labels = c("Physiological Aging", "Accelerated AnthropoAge"))) %>%
  ggplot(aes(x=accel1, y=Age, color=accel1))+scale_color_manual(values=colors_in)+geom_point(size=10*0.4)+
  theme(legend.position = "top",legend.title=element_text(size=18*0.4), legend.text=element_text(size=22*0.4),
        legend.key.size = unit(2*0.4,"cm"), legend.background = element_rect(fill="transparent"))+
  labs(color=NULL); legend<-get_legend(g1)

fig6<-ggarrange(ggarrange(
  annotate_figure(ggarrange(f6a,labels = "a", font.label = list(size=25*0.4)),
                  left=text_grob("Women", face = "bold", size=25*0.4, rot = 90),
                  top=text_grob("Anthropometry", face = "bold", size=25*0.4)),
  annotate_figure(ggarrange(f6b,labels = "b", font.label = list(size=25*0.4)),
                  top=text_grob("DXA", face = "bold", size=25*0.4)),
  annotate_figure(ggarrange(f6c,labels = "c", font.label = list(size=25*0.4)),
                  top=text_grob("PhenoAge", face = "bold", size=25*0.4)), ncol=3, nrow=1),
  annotate_figure(ggarrange(f6d,f6e,f6f,labels = letters[4:6], font.label = list(size=25*0.4), nrow=1, ncol=3),
                  left=text_grob("Men", face = "bold", size=25*0.4, rot = 90)), nrow=2, ncol=1)+
  geom_subview(x=0.5, y=1-0.985, subview=legend)

ggsave(file = "Figure4.jpg", fig6, bg = "transparent", width = 28.75, 
       height = 18.9, units=c("cm"), dpi = 500, limitsize = FALSE)
#Submission
ggsave(fig6, file="Submission/Figures/Figure_6.pdf", bg="transparent",
       width=28.75, height=18.9, units=c("cm"), dpi=600, limitsize = FALSE)



####---- S-AnthropoAgeAccel>0 Spiderplots ----####
### Spiderplot: AnthropoAge components ###
#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Age~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a1
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Leg_length~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a2
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Weight~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a3
(with(dxa%>%filter(Sex=="Women"), wilcox.test(BMI~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a4
(with(dxa%>%filter(Sex=="Women"), wilcox.test(ICE~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a5
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Subscapular_skinfold~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a6
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Triceps_skinfold~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a7
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Arm_circumference~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a8
(with(dxa%>%filter(Sex=="Women"), wilcox.test(Thigh_circumference~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a9
(with(dxa%>%filter(Sex=="Women"), wilcox.test(AnthropoAge2~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("CA","Leg\nlength", "Weight", "BMI", "WHtR","Subscapular\nskinfold",
                  "Triceps\nskinfold", "Arm\ncircumference","Thigh\ncircumference","S-AnthropoAge"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data1 <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
             dplyr::select(Age,Leg_length,Weight,BMI,ICE,Subscapular_skinfold,Triceps_skinfold,
                           Arm_circumference,Thigh_circumference,AnthropoAge2,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4a<-ggplotify::as.ggplot(~radarchart(data1, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Age~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a1
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Leg_length~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a2
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Weight~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a3
(with(dxa%>%filter(Sex=="Men"), wilcox.test(BMI~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a4
(with(dxa%>%filter(Sex=="Men"), wilcox.test(ICE~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a5
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Subscapular_skinfold~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a6
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Triceps_skinfold~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a7
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Arm_circumference~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a8
(with(dxa%>%filter(Sex=="Men"), wilcox.test(Thigh_circumference~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a9
(with(dxa%>%filter(Sex=="Men"), wilcox.test(AnthropoAge2~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("CA","Leg\nlength", "Weight", "BMI", "WHtR","Subscapular\nskinfold",
                  "Triceps\nskinfold", "Arm\ncircumference","Thigh\ncircumference","S-AnthropoAge"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data2 <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
             dplyr::select(Age,Leg_length,Weight,BMI,ICE,Subscapular_skinfold,Triceps_skinfold,
                           Arm_circumference,Thigh_circumference,AnthropoAge2,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4d<-ggplotify::as.ggplot(~radarchart(data2, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))



### Spiderplot: DXA ###
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_fmi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a1
(with(dxa%>%filter(Sex=="Women"), wilcox.test(trunk_fmi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a2
(with(dxa%>%filter(Sex=="Women"), wilcox.test(afmi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a3
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_lmi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a4
(with(dxa%>%filter(Sex=="Women"), wilcox.test(almi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a5
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_ratio~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a6
(with(dxa%>%filter(Sex=="Women"), wilcox.test(trunk_ratio~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a7
(with(dxa%>%filter(Sex=="Women"), wilcox.test(ttafm_ratio~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a8
(with(dxa%>%filter(Sex=="Women"), wilcox.test(lumbar_bmd~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a9
(with(dxa%>%filter(Sex=="Women"), wilcox.test(total_bmd~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("Body\nFMI", "Trunk\nFMI", "AFMI", "Body\nLMI", "ALMI", "Body fat-lean ratio",
                  "Trunk fat-lean\nratio", "Trunk-Append\nfat ratio", "Lumbar BMD","Total BMD"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data3 <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
             dplyr::select(total_fmi, trunk_fmi, afmi, total_lmi, almi, total_ratio,
                           trunk_ratio, ttafm_ratio, lumbar_bmd, total_bmd, AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4b<-ggplotify::as.ggplot(~radarchart(data3, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_fmi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a1
(with(dxa%>%filter(Sex=="Men"), wilcox.test(trunk_fmi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a2
(with(dxa%>%filter(Sex=="Men"), wilcox.test(afmi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a3
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_lmi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a4
(with(dxa%>%filter(Sex=="Men"), wilcox.test(almi~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a5
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_ratio~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a6
(with(dxa%>%filter(Sex=="Men"), wilcox.test(trunk_ratio~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a7
(with(dxa%>%filter(Sex=="Men"), wilcox.test(ttafm_ratio~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a8
(with(dxa%>%filter(Sex=="Men"), wilcox.test(lumbar_bmd~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a9
(with(dxa%>%filter(Sex=="Men"), wilcox.test(total_bmd~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("Body\nFMI", "Trunk\nFMI", "AFMI", "Body\nLMI", "ALMI", "Body fat-lean ratio",
                  "Trunk fat-lean\nratio", "Trunk-Append\nfat ratio", "Lumbar BMD","Total BMD"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data4 <- ((accel1 <- dxa %>% filter(Sex=="Men"&Age>=75) %>% 
             dplyr::select(total_fmi, trunk_fmi, afmi, total_lmi, almi, total_ratio,
                           trunk_ratio, ttafm_ratio, lumbar_bmd, total_bmd, AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4e<-ggplotify::as.ggplot(~radarchart(data4, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))



### Spiderplots: PhenoAge components ###
#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(pheno%>%filter(Sex=="Women"), wilcox.test(PhenoAge~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a1
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Glucose~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a2
(with(pheno%>%filter(Sex=="Women"), t.test(Creatinine~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a3
(with(pheno%>%filter(Sex=="Women"), wilcox.test(Albumin~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a4
(with(pheno%>%filter(Sex=="Women"), wilcox.test(ALP~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a5
(with(pheno%>%filter(Sex=="Women"), wilcox.test(WBC~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a6
(with(pheno%>%filter(Sex=="Women"), wilcox.test(LymphP~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a7
(with(pheno%>%filter(Sex=="Women"), wilcox.test(MCV~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a8
(with(pheno%>%filter(Sex=="Women"), wilcox.test(RDW~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a9
(with(pheno%>%filter(Sex=="Women"), wilcox.test(CRP~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("PhenoAge","Glucose","Creatinine","Albumin","ALP","WBC",
                  "Lymphocyte\nPercentage","MCV","RDW","CRP"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data5 <- ((accel1 <- pheno %>% filter(Sex=="Women") %>% 
             dplyr::select(PhenoAge,Glucose,Creatinine,Albumin,ALP,WBC,LymphP,MCV,RDW,CRP,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
            `rownames<-`(c("Phys-Aging", "Accel-Aging")) %>% rbind(rep(1,10), rep(-1,10)) %>% slice(3:4,1:2))

s4c<-ggplotify::as.ggplot(~radarchart(data5, axistype=1 , seg=2, pcol=colors_border, pfcol=colors_in,
                                      plwd=4 , plty=1, title="", cex.main=1.5, cglcol="grey30", cglty=2, axislabcol="gray35",
                                      caxislabels=seq(-1, 1, 1), cglwd=0.8, vlcex=1.65))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#1* p<0.05, 2*=p<0.01, 3*=p<0.001
(with(pheno%>%filter(Sex=="Men"), wilcox.test(PhenoAge~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a1
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Glucose~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a2
(with(pheno%>%filter(Sex=="Men"), t.test(Creatinine~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a3
(with(pheno%>%filter(Sex=="Men"), wilcox.test(Albumin~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a4
(with(pheno%>%filter(Sex=="Men"), wilcox.test(ALP~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a5
(with(pheno%>%filter(Sex=="Men"), wilcox.test(WBC~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a6
(with(pheno%>%filter(Sex=="Men"), wilcox.test(LymphP~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a7
(with(pheno%>%filter(Sex=="Men"), wilcox.test(MCV~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a8
(with(pheno%>%filter(Sex=="Men"), wilcox.test(RDW~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a9
(with(pheno%>%filter(Sex=="Men"), wilcox.test(CRP~AnthropoAgeAccel2>0)))$p.value %>% p_aster -> p_a10

spider_names <- c("PhenoAge","Glucose","Creatinine","Albumin","ALP","WBC",
                  "Lymphocyte\nPercentage","MCV","RDW","CRP"
) %>% paste0((as.list(paste0("p_a",1:10)) %>% lapply(get) %>% as.character()))

data6 <- ((accel1 <- pheno %>% filter(Sex=="Men") %>% 
             dplyr::select(PhenoAge,Glucose,Creatinine,Albumin,ALP,WBC,LymphP,MCV,RDW,CRP,AnthropoAgeAccel2) %>% 
             scale() %>% as.data.frame() %>%  mutate("accel"=as.numeric(AnthropoAgeAccel2>0))) %>%
            select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% `colnames<-`(spider_names) %>%
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

ggsave(file = "SuppFig7.jpg", supp4, bg = "transparent", width = 62.5*1.15,
       height = 45*1.05, units=c("cm"), dpi = 500, limitsize = FALSE)



####---- Multidimensional Spiderplots ----####
### Spiderplot: AnthropoAge components ###
### Spiderplots: Anthropometry ###
colors_border2<-c("#80C4E7","#C29300","#AA4D10","#911F2C")
data1F <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
              dplyr::select(Age,Leg_length,Weight,BMI,ICE,Subscapular_skinfold,Triceps_skinfold,
                            Arm_circumference,Thigh_circumference,AnthropoAge,AnthropoAgeAccel,PhenoAgeAccel) %>% scale() %>%
              as.data.frame() %>%  mutate("accelA"=as.numeric(AnthropoAgeAccel>0),"accelP"=as.numeric(PhenoAgeAccel>0)) %>% 
              mutate("accel"=(2*(accelP>0) + (accelA>0)))) %>%
             select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
             `colnames<-`(c("CA","Leg\nlength", "Weight", "BMI", "WHtR","Subscapular\nskinfold",
                            "Triceps\nskinfold", "Arm\ncircumference","Thigh\ncircumference","AnthropoAge")) %>%
             `rownames<-`(c("Phys-aging", "Accel-AnthropoAge", "Accel-PhenoAge", "Multi-accel")) %>%
             rbind(rep(0.7,10), rep(-0.7,10)) %>% slice(5:6,1:4))

data1M <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
              dplyr::select(Age,Leg_length,Weight,BMI,ICE,Subscapular_skinfold,Triceps_skinfold,
                            Arm_circumference,Thigh_circumference,AnthropoAge,AnthropoAgeAccel,PhenoAgeAccel) %>% scale() %>%
              as.data.frame() %>%  mutate("accelA"=as.numeric(AnthropoAgeAccel>0),"accelP"=as.numeric(PhenoAgeAccel>0)) %>% 
              mutate("accel"=(2*(accelP>0) + (accelA>0)))) %>%
             select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
             `colnames<-`(c("CA","Leg\nlength", "Weight", "BMI", "WHtR","Subscapular\nskinfold",
                            "Triceps\nskinfold", "Arm\ncircumference","Thigh\ncircumference","AnthropoAge")) %>%
             `rownames<-`(c("Phys-aging", "Accel-AnthropoAge", "Accel-PhenoAge", "Multi-accel")) %>%
             rbind(rep(0.7,10), rep(-0.7,10)) %>% slice(5:6,1:4))

fYa<-ggplotify::as.ggplot(~radarchart(data1F, axistype=1 , seg=2, pcol=colors_border2, pfcol=NULL,
                                      plwd=c(7.5) , plty= c("solid","F3","dashed","dotted"), title="", cex.main=1.5, cglcol="grey30", cglty=1, axislabcol="gray35",
                                      caxislabels=seq(-0.7, 0.7, 0.7), cglwd=0.8, vlcex=2, pty=NA))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))

fYd<-ggplotify::as.ggplot(~radarchart(data1M, axistype=1 , seg=2, pcol=colors_border2, pfcol=NULL,
                                      plwd=c(7.5) , plty= c("solid","F3","dashed","dotted"), title="", cex.main=1.5, cglcol="grey30", cglty=1, axislabcol="gray35",
                                      caxislabels=seq(-0.7, 0.7, 0.7), cglwd=0.8, vlcex=2, pty=NA))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


### Spiderplots: DXA ###
data2F <- ((accel1 <- dxa %>% filter(Sex=="Women") %>% 
              dplyr::select(total_fmi, trunk_fmi, afmi, total_lmi, almi, total_ratio, trunk_ratio,
                            ttafm_ratio, lumbar_bmd, total_bmd, AnthropoAgeAccel,PhenoAgeAccel) %>% 
              scale() %>% as.data.frame() %>%  mutate("accelA"=as.numeric(AnthropoAgeAccel>0),
                                                      "accelP"=as.numeric(PhenoAgeAccel>0)) %>% 
              mutate("accel"=(2*(accelP>0) + (accelA>0)))) %>%
             select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
             `colnames<-`(c("Body\nFMI", "Trunk\nFMI", "AFMI", "Body\nLMI", "ALMI", "Body fat-lean ratio",
                            "Trunk fat-lean\nratio", "Trunk-Append\nfat ratio", "Lumbar BMD","Total BMD")) %>%
             `rownames<-`(c("Phys-aging", "Accel-AnthropoAge", "Accel-PhenoAge", "Multi-accel")) %>%
             rbind(rep(0.7,10), rep(-0.7,10)) %>% slice(5:6,1:4))

data2M <- ((accel1 <- dxa %>% filter(Sex=="Men") %>% 
              dplyr::select(total_fmi, trunk_fmi, afmi, total_lmi, almi, total_ratio, trunk_ratio,
                            ttafm_ratio, lumbar_bmd, total_bmd, AnthropoAgeAccel,PhenoAgeAccel) %>% 
              scale() %>% as.data.frame() %>%  mutate("accelA"=as.numeric(AnthropoAgeAccel>0),
                                                      "accelP"=as.numeric(PhenoAgeAccel>0)) %>% 
              mutate("accel"=(2*(accelP>0) + (accelA>0)))) %>%
             select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
             `colnames<-`(c("Body\nFMI", "Trunk\nFMI", "AFMI", "Body\nLMI", "ALMI", "Body fat-lean ratio",
                            "Trunk fat-lean\nratio", "Trunk-Append\nfat ratio", "Lumbar BMD","Total BMD")) %>%
             `rownames<-`(c("Phys-aging", "Accel-AnthropoAge", "Accel-PhenoAge", "Multi-accel")) %>%
             rbind(rep(0.7,10), rep(-0.7,10)) %>% slice(5:6,1:4))

fYb<-ggplotify::as.ggplot(~radarchart(data2F, axistype=1 , seg=2, pcol=colors_border2, pfcol=NULL,
                                      plwd=c(7.5) , plty= c("solid","F3","dashed","dotted"), title="", cex.main=1.5, cglcol="grey30", cglty=1, axislabcol="gray35",
                                      caxislabels=seq(-0.7, 0.7, 0.7), cglwd=0.8, vlcex=2, pty=NA))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))

fYe<-ggplotify::as.ggplot(~radarchart(data2M, axistype=1 , seg=2, pcol=colors_border2, pfcol=NULL,
                                      plwd=c(7.5) , plty= c("solid","F3","dashed","dotted"), title="", cex.main=1.5, cglcol="grey30", cglty=1, axislabcol="gray35",
                                      caxislabels=seq(-0.7, 0.7, 0.7), cglwd=0.8, vlcex=2, pty=NA))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))



### Spiderplots: PhenoAge components ###
data3F <- ((accel1 <- pheno %>% filter(Sex=="Women") %>% 
              dplyr::select(PhenoAge,Glucose,Creatinine,Albumin,ALP,WBC,LymphP,MCV,RDW,CRP,AnthropoAgeAccel, PhenoAgeAccel) %>% 
              scale() %>% as.data.frame() %>%  mutate("accelA"=as.numeric(AnthropoAgeAccel>0),"accelP"=as.numeric(PhenoAgeAccel>0)) %>% 
              mutate("accel"=(2*(accelP>0) + (accelA>0)))) %>%
             select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
             `colnames<-`(c("PhenoAge","Glucose","Creatinine","Albumin","ALP","WBC",
                            "Lymphocyte\nPercentage","MCV","RDW","CRP")) %>%
             `rownames<-`(c("Phys-aging", "Accel-AnthropoAge", "Accel-PhenoAge", "Multi-accel")) %>%
             rbind(rep(0.7,10), rep(-0.7,10)) %>% slice(5:6,1:4))

data3M <- ((accel1 <- pheno %>% filter(Sex=="Men") %>% 
              dplyr::select(PhenoAge,Glucose,Creatinine,Albumin,ALP,WBC,LymphP,MCV,RDW,CRP,AnthropoAgeAccel, PhenoAgeAccel) %>% 
              scale() %>% as.data.frame() %>%  mutate("accelA"=as.numeric(AnthropoAgeAccel>0),"accelP"=as.numeric(PhenoAgeAccel>0)) %>% 
              mutate("accel"=(2*(accelP>0) + (accelA>0)))) %>%
             select(1:10) %>% aggregate(list(accel1$accel), median) %>% select(2:11) %>% 
             `colnames<-`(c("PhenoAge","Glucose","Creatinine","Albumin","ALP","WBC",
                            "Lymphocyte\nPercentage","MCV","RDW","CRP")) %>%
             `rownames<-`(c("Phys-aging", "Accel-AnthropoAge", "Accel-PhenoAge", "Multi-accel")) %>%
             rbind(rep(0.7,10), rep(-0.7,10)) %>% slice(5:6,1:4))

fYc<-ggplotify::as.ggplot(~radarchart(data3F, axistype=1 , seg=2, pcol=colors_border2, pfcol=NULL,
                                      plwd=c(7.5) , plty= c("solid","F3","dashed","dotted"), title="", cex.main=1.5, cglcol="grey30", cglty=1, axislabcol="gray35",
                                      caxislabels=seq(-0.7, 0.7, 0.7), cglwd=0.8, vlcex=2, pty=NA))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))

fYf<-ggplotify::as.ggplot(~radarchart(data3M, axistype=1 , seg=2, pcol=colors_border2, pfcol=NULL,
                                      plwd=c(7.5) , plty= c("solid","F3","dashed","dotted"), title="", cex.main=1.5, cglcol="grey30", cglty=1, axislabcol="gray35",
                                      caxislabels=seq(-0.7, 0.7, 0.7), cglwd=0.8, vlcex=2, pty=NA))+
  theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


GGG<-accel1 %>% mutate(accel=factor(accel)) %>% ggplot(aes(AnthropoAgeAccel,PhenoAgeAccel,color=accel,linetype=accel)) +
  geom_smooth(method="lm", se=F, size=2.5) + labs(color="",linetype="") +
  scale_color_manual(values=colors_border2, labels=c("Phys-aging", "Accel\nAnthropoAge", "Accel\nPhenoAge", "Multi-accel")) +
  scale_linetype_manual(values=c("solid","longdash","dashed","dotted"), labels=c("Phys-aging", "Accel\nAnthropoAge", "Accel\nPhenoAge", "Multi-accel"))+ 
  theme(legend.position = "top",legend.title=element_text(size=20*.90),
        legend.text=element_text(size=22),legend.key.size = unit(2,"cm")); legend2 <- get_legend(GGG)


figY<-ggarrange(ggarrange(annotate_figure(ggarrange(fYa,labels = "a", font.label = list(size=25)), left=text_grob("Women", face = "bold", size=25, rot = 90),top=text_grob("Anthropometry", face = "bold", size=25)),
                          annotate_figure(ggarrange(fYb,labels = "b", font.label = list(size=25)), top=text_grob("DXA", face = "bold", size=25)),
                          annotate_figure(ggarrange(fYc,labels = "c", font.label = list(size=25)), top=text_grob("PhenoAge", face = "bold", size=25)), ncol=3, nrow=1),
                annotate_figure(ggarrange(fYd,fYe,fYf,labels = letters[4:6], font.label = list(size=25), nrow=1, ncol=3), left=text_grob("Men", face = "bold", size=25, rot = 90)), nrow=2, ncol=1)+
  geom_subview(x=0.5, y=1-0.985, subview=legend2)

ggsave(file = "SuppFig8.jpg", figY, bg = "transparent", width = 62.5*1.15,
       height = 45*1.05, units=c("cm"), dpi = 500, limitsize = FALSE)
