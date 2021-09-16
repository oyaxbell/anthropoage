library(haven); library(ggplot2); library(ggpubr); library(tidyverse); 
library(gtools); library(data.table); library(nhanesA); library(haven); library(ggplot2); library(ggpubr); library(tidyverse); 
library(gtools); library(data.table); library(mice)

## Lista de tablas ##
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/AnthropoAge")


# DXA #
DXA1999_2000 <- nhanesDXA(2000)
DXA2001_2002 <- nhanesDXA(2002)
DXA2003_2004 <- nhanesDXA(2004)
DXA2005_2006 <- nhanesDXA(2006)
DXA2011_2012 <- nhanes("DXX_G")
DXA2013_2014 <- nhanes("DXX_H")
DXA2015_2016 <- nhanes("DXX_I")
DXA2017_2018 <-nhanes("DXX_J")


DXA_Antro <- smartbind(DXA1999_2000, DXA2001_2002, DXA2003_2004, DXA2005_2006, DXA2011_2012, DXA2013_2014, DXA2015_2016, DXA2017_2018)
view(DXA_Antro)
write_csv(DXA_Antro,"DXA_Antro.csv") # Aqu? pongan ustedes la localizaci?n que quieren.

# Laboratorios #
labs2000names  <- nhanesTables('LAB', 2000, namesonly=TRUE)
labs2000 <- lapply(labs2000names, nhanes)
names(labs2000) <- labs2000names

labs2002names  <- nhanesTables('LAB', 2002, namesonly=TRUE)
labs2002 <- lapply(labs2002names, nhanes)
names(labs2002) <- labs2002names

labs2004names  <- nhanesTables('LAB', 2004, namesonly=TRUE)
labs2004 <- lapply(labs2004names, nhanes)
names(labs2004) <- labs2004names

labs2006names  <- nhanesTables('LAB', 2006, namesonly=TRUE)
labs2006 <- lapply(labs2006names, nhanes)
names(labs2006) <- labs2006names

labs2008names  <- nhanesTables('LAB', 2008, namesonly=TRUE)
labs2008 <- lapply(labs2008names, nhanes)
names(labs2008) <- labs2008names

labs2010names  <- nhanesTables('LAB', 2010, namesonly=TRUE)
labs2010 <- lapply(labs2010names, nhanes)
names(labs2010) <- labs2010names

labs2012names  <- nhanesTables('LAB', 2012, namesonly=TRUE)
labs2012 <- lapply(labs2012names, nhanes)
names(labs2012) <- labs2012names

labs2014names  <- nhanesTables('LAB', 2014, namesonly=TRUE)
labs2014 <- lapply(labs2014names, nhanes)
names(labs2014) <- labs2014names

labs2016names  <- nhanesTables('LAB', 2016, namesonly=TRUE)
labs2016 <- lapply(labs2016names, nhanes)
names(labs2016) <- labs2016names

labs2018names  <- nhanesTables('LAB', 2018, namesonly=TRUE)
labs2018 <- lapply(labs2018names, nhanes)
names(labs2018) <- labs2018names

labs_Antro <- list(labs2000, labs2002, labs2004, labs2006, labs2008,labs2010,labs2012, labs2014, labs2016, labs2018)

# Demogr?ficos #
demo_A <- nhanes("DEMO")
demo_B <- nhanes("DEMO_B")
demo_C <- nhanes("DEMO_C")
demo_D <- nhanes("DEMO_D")
demo_E <- nhanes("DEMO_E")
demo_F <- nhanes("DEMO_F")
demo_G <- nhanes("DEMO_G")
demo_H <- nhanes("DEMO_H")
demo_I <- nhanes("DEMO_I")
demo_J <- nhanes("DEMO_J")
demo <- smartbind(demo_A, demo_B, demo_C, demo_D,demo_E,demo_F, demo_G, demo_H, demo_I, demo_J)
write_csv(demo, "demo_nhanes.csv") # Aqu? pongan ustedes la localizaci?n que quieren.


# Salud
health_A <- nhanes("MCQ")
health_B <- nhanes("MCQ_B")
health_C <- nhanes("MCQ_C")
health_D <- nhanes("MCQ_D")
health_E <- nhanes("MCQ_E")
health_F <- nhanes("MCQ_F")
health_G <- nhanes("MCQ_G")
health_H <- nhanes("MCQ_H")
health_I <- nhanes("MCQ_I")
health_J <- nhanes("MCQ_J")
health <- smartbind(health_A, health_B, health_C, health_D,health_E,health_F, health_G, health_H, health_I, health_J)
write_csv(health, "health_nhanes.csv") 

# Diabetes
diabetes_A <- nhanes("DIQ")
diabetes_B <- nhanes("DIQ_B")
diabetes_C <- nhanes("DIQ_C")
diabetes_D <- nhanes("DIQ_D")
diabetes_E <- nhanes("DIQ_E")
diabetes_F <- nhanes("DIQ_F")
diabetes_G <- nhanes("DIQ_G")
diabetes_H <- nhanes("DIQ_H")
diabetes_I <- nhanes("DIQ_I")
diabetes_J <- nhanes("DIQ_J")
diabetes <- smartbind(diabetes_A, diabetes_B, diabetes_C, diabetes_D,diabetes_E,diabetes_F, diabetes_G, diabetes_H, diabetes_I, diabetes_J)
write_csv(diabetes, "diabetes_nhanes.csv") 

# Antropometr?a #
BMX <- nhanes("BMX")
BMX_B <- nhanes("BMX_B")
BMX_C <- nhanes("BMX_C")
BMX_D <- nhanes("BMX_D")
BMX_E <- nhanes("BMX_E")
BMX_F <- nhanes("BMX_F")
BMX_G <- nhanes("BMX_G")
BMX_H <- nhanes("BMX_H")
BMX_I <- nhanes("BMX_I")
BMX_J <- nhanes("BMX_J")

antro <- smartbind(BMX, BMX_B, BMX_C, BMX_D,BMX_E,BMX_F, BMX_G, BMX_H, BMX_I, BMX_J)
write_csv(antro, "antro_nhanes.csv") # Aqu? pongan ustedes la localizaci?n que quieren.

#### Variables Edad Fenot?pica ####
# Variables DEMOGR?FICAS
# 2000 - https://wwwn.cdc.gov/nchs/nhanes/Search/variablelist.aspx?Component=Demographics&CycleBeginYear=1999
# 2002 - https://wwwn.cdc.gov/nchs/nhanes/Search/variablelist.aspx?Component=Demographics&CycleBeginYear=2001
# 2004 - https://wwwn.cdc.gov/nchs/nhanes/Search/variablelist.aspx?Component=Demographics&CycleBeginYear=2003
# 2006 - https://wwwn.cdc.gov/nchs/nhanes/Search/variablelist.aspx?Component=Demographics&CycleBeginYear=2005
# 2012 - https://wwwn.cdc.gov/nchs/nhanes/Search/variablelist.aspx?Component=Demographics&CycleBeginYear=2011
# 2014 - https://wwwn.cdc.gov/nchs/nhanes/Search/variablelist.aspx?Component=Demographics&CycleBeginYear=2013
# 2016 - https://wwwn.cdc.gov/nchs/nhanes/Search/variablelist.aspx?Component=Demographics&CycleBeginYear=2015
# 2018 - https://wwwn.cdc.gov/nchs/nhanes/Search/variablelist.aspx?Component=Demographics&CycleBeginYear=2017

# Variables de lABORATORIO
# 2000 - https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&CycleBeginYear=1999
# 2002 - https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&CycleBeginYear=2001
# 2004 - https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&CycleBeginYear=2003
# 2006 - https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&CycleBeginYear=2005
# 2012 - https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&CycleBeginYear=2011
# 2014 - https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&CycleBeginYear=2013
# 2016 - https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&CycleBeginYear=2015
# 2018 - https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&CycleBeginYear=2017
  
NHANES = data.frame(SEQN=c(demo_A$SEQN,demo_B$SEQN,demo_C$SEQN,demo_D$SEQN,demo_E$SEQN,demo_F$SEQN,demo_G$SEQN,demo_H$SEQN,demo_I$SEQN,demo_J$SEQN),
                    Sex=c(demo_A$RIAGENDR,demo_B$RIAGENDR,demo_C$RIAGENDR,demo_D$RIAGENDR,demo_E$RIAGENDR,demo_F$RIAGENDR,demo_G$RIAGENDR,demo_H$RIAGENDR,demo_I$RIAGENDR,demo_J$RIAGENDR),
                    Age=c(demo_A$RIDAGEYR,demo_B$RIDAGEYR,demo_C$RIDAGEYR,demo_D$RIDAGEYR,demo_E$RIDAGEYR,demo_F$RIDAGEYR,demo_G$RIDAGEYR,demo_H$RIDAGEYR,demo_I$RIDAGEYR,demo_J$RIDAGEYR),
                    Ethnicity=c(demo_A$RIDRETH1, demo_B$RIDRETH1,demo_C$RIDRETH1,demo_D$RIDRETH1,demo_E$RIDRETH1,demo_F$RIDRETH1,demo_G$RIDRETH1, demo_H$RIDRETH1,demo_I$RIDRETH1,demo_J$RIDRETH1))
write_csv(NHANES, "demo_short_nhanes.csv")
NHANES$SEQN<-as.character(NHANES$SEQN)


comorb=data.frame(SEQN=c(health_A$SEQN,health_B$SEQN,health_C$SEQN,health_D$SEQN,health_E$SEQN,health_F$SEQN,health_G$SEQN,health_H$SEQN,health_I$SEQN,health_J$SEQN),
                  Diabetes=c(diabetes_A$DIQ010,diabetes_B$DIQ010,diabetes_C$DIQ010,diabetes_D$DIQ010,diabetes_E$DIQ010,diabetes_F$DIQ010,diabetes_G$DIQ010,diabetes_H$DIQ010,diabetes_I$DIQ010,diabetes_J$DIQ010),
                  Asthma=c(health_A$MCQ010,health_B$MCQ010,health_C$MCQ010,health_D$MCQ010,health_E$MCQ010,health_F$MCQ010,health_G$MCQ010,health_H$MCQ010,health_I$MCQ010,health_J$MCQ010),
                 Overweight=c(health_A$MCQ080,health_B$MCQ080,health_C$MCQ080,health_D$MCQ080,health_E$MCQ080,health_F$MCQ080,health_G$MCQ080,health_H$MCQ080,health_I$MCQ080,health_J$MCQ080),
                 Arthritis=c(health_A$MCQ160F,health_B$MCQ160F,health_C$MCQ160F,health_D$MCQ160F,health_E$MCQ160F,health_F$MCQ160F,health_G$MCQ160F,health_H$MCQ160F,health_I$MCQ160F,health_J$MCQ160F),
                 Heart_Failure=c(health_A$MCQ160B,health_B$MCQ160B,health_C$MCQ160B,health_D$MCQ160B,health_E$MCQ160B,health_F$MCQ160B,health_G$MCQ160B,health_H$MCQ160B,health_I$MCQ160B,health_J$MCQ160B),
                 CAD=c(health_A$MCQ160C,health_B$MCQ160C,health_C$MCQ160C,health_D$MCQ160C,health_E$MCQ160C,health_F$MCQ160C,health_G$MCQ160C,health_H$MCQ160C,health_I$MCQ160C,health_J$MCQ160C),
                 Angina=c(health_A$MCQ160D,health_B$MCQ160D,health_C$MCQ160D,health_D$MCQ160D,health_E$MCQ160D,health_F$MCQ160D,health_G$MCQ160D,health_H$MCQ160D,health_I$MCQ160D,health_J$MCQ160D),
                 Heart_attack=c(health_A$MCQ160E,health_B$MCQ160E,health_C$MCQ160E,health_D$MCQ160E,health_E$MCQ160E,health_F$MCQ160E,health_G$MCQ160E,health_H$MCQ160E,health_I$MCQ160E,health_J$MCQ160E),
                 Stroke=c(health_A$MCQ160F,health_B$MCQ160F,health_C$MCQ160F,health_D$MCQ160F,health_E$MCQ160F,health_F$MCQ160F,health_G$MCQ160F,health_H$MCQ160F,health_I$MCQ160F,health_J$MCQ160F),
                 Emphysema=c(health_A$MCQ160G,health_B$MCQ160G,health_C$MCQ160G,health_D$MCQ160G,health_E$MCQ160G,health_F$MCQ160G,health_G$MCQ160G,health_H$MCQ160G,health_I$MCQ160G,health_J$MCQ160G),
                 Bronchitis=c(health_A$MCQ160K,health_B$MCQ160K,health_C$MCQ160K,health_D$MCQ160K,health_E$MCQ160K,health_F$MCQ160K,health_G$MCQ160K,health_H$MCQ160K,health_I$MCQ160K,health_J$MCQ160K),
                 Liver_Disease=c(health_A$MCQ160L,health_B$MCQ160L,health_C$MCQ160L,health_D$MCQ160L,health_E$MCQ160L,health_F$MCQ160L,health_G$MCQ160L,health_H$MCQ160L,health_I$MCQ160L,health_J$MCQ160L),
                 Malignancy=c(health_A$MCQ220,health_B$MCQ220,health_C$MCQ220,health_D$MCQ220,health_E$MCQ220,health_F$MCQ220,health_G$MCQ220,health_H$MCQ220,health_I$MCQ220,health_J$MCQ220))
                 
write_csv(NHANES, "comorb_short_nhanes.csv")
NHANES$SEQN<-as.character(NHANES$SEQN)


# Ethnicity2 ####
df2 = data.frame(SEQN=c(demo_G$SEQN, demo_H$SEQN,demo_I$SEQN,demo_J$SEQN),
                 Ethnicity2=c(demo_G$RIDRETH3, demo_H$RIDRETH3,demo_I$RIDRETH3,demo_J$RIDRETH3))

# Glucose ####
df2 = data.frame(SEQN=c(labs2000$LAB10AM$SEQN, labs2002$L10AM_B$SEQN, labs2004$L10AM_C$SEQN, labs2006$GLU_D$SEQN, 
                        labs2008$GLU_E$SEQN,labs2010$GLU_F$SEQN,labs2012$GLU_G$SEQN, labs2014$GLU_H$SEQN, labs2016$GLU_I$SEQN, labs2018$BIOPRO_J$SEQN),
                 Glucose=c(labs2000$LAB10AM$LBXGLU, labs2002$L10AM_B$LBXGLU, labs2004$L10AM_C$LBXGLU, labs2006$GLU_D$LBXGLU, 
                           labs2008$GLU_E$LBXGLU,labs2010$GLU_F$LBXGLU,labs2012$GLU_G$LBXGLU, labs2014$GLU_H$LBXGLU, labs2016$GLU_I$LBXGLU, labs2018$BIOPRO_J$LBXSGL))
df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")

# Alb?mina ####
df2 = data.frame(SEQN=c(labs2000$LAB18$SEQN, labs2002$L40_B$SEQN, labs2004$L40_C$SEQN, labs2006$BIOPRO_D$SEQN, 
                        labs2008$BIOPRO_E$SEQN,labs2010$BIOPRO_F$SEQN,labs2012$BIOPRO_G$SEQN, labs2014$BIOPRO_H$SEQN, labs2016$BIOPRO_I$SEQN, labs2018$BIOPRO_J$SEQN),
                 Albumin=c(labs2000$LAB18$LBXSAL, labs2002$L40_B$LBXSAL, labs2004$L40_C$LBXSAL, labs2006$BIOPRO_D$LBXSAL, 
                           labs2008$BIOPRO_E$LBXSAL,labs2010$BIOPRO_F$LBXSAL,labs2012$BIOPRO_G$LBXSAL, labs2014$BIOPRO_H$LBXSAL, labs2016$BIOPRO_I$LBXSAL, labs2018$BIOPRO_J$LBXSAL))
df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")

# Creatinina ####
df2 = data.frame(SEQN=c(labs2000$LAB18$SEQN, labs2002$L40_B$SEQN, labs2004$L40_C$SEQN, labs2006$BIOPRO_D$SEQN, 
                        labs2008$BIOPRO_E$SEQN,labs2010$BIOPRO_F$SEQN,labs2012$BIOPRO_G$SEQN, labs2014$BIOPRO_H$SEQN, labs2016$BIOPRO_I$SEQN, labs2018$BIOPRO_J$SEQN),
                 Cr=c(labs2000$LAB18$LBXSCR, labs2002$L40_B$LBDSCR, labs2004$L40_C$LBXSCR, labs2006$BIOPRO_D$LBXSCR, 
                      labs2008$BIOPRO_E$LBXSCR,labs2010$BIOPRO_F$LBXSCR,labs2012$BIOPRO_G$LBXSCR, labs2014$BIOPRO_H$LBXSCR, labs2016$BIOPRO_I$LBXSCR, labs2018$BIOPRO_J$LBXSCR))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Prote?na C reactiva ####

df2 = data.frame(SEQN=c(labs2000$LAB11$SEQN, labs2002$L11_B$SEQN, labs2004$L11_C$SEQN, labs2006$CRP_D$SEQN, labs2008$CRP_E$SEQN,
                        labs2010$CRP_F$SEQN,labs2016$HSCRP_I$SEQN, labs2018$HSCRP_J$SEQN),
                 CRP=c(labs2000$LAB11$LBXCRP, labs2002$L11_B$LBXCRP, labs2004$L11_C$LBXCRP, labs2006$CRP_D$LBXCRP, labs2008$CRP_E$LBXCRP,
                       labs2010$CRP_F$LBXCRP,labs2016$HSCRP_I$LBXHSCRP, labs2018$HSCRP_J$LBXHSCRP))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Porcentaje de linfocitos ####
df2 = data.frame(SEQN=c(labs2000$LAB25$SEQN, labs2002$L25_B$SEQN, labs2004$L25_C$SEQN, labs2006$CBC_D$SEQN, labs2008$CBC_E$SEQN,
                        labs2012$CBC_G$SEQN, labs2014$CBC_H$SEQN, labs2016$CBC_I$SEQN, labs2018$CBC_J$SEQN),
                 LymP=c(labs2000$LAB25$LBXLYPCT, labs2002$L25_B$LBXLYPCT, labs2004$L25_C$LBXLYPCT, labs2006$CBC_D$LBXLYPCT, labs2008$CBC_E$LBXLYPCT,
                        labs2012$CBC_G$LBXLYPCT, labs2014$CBC_H$LBXLYPCT, labs2016$CBC_I$LBXLYPCT, labs2018$CBC_J$LBXLYPCT))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Volumen celular medio ####

df2 = data.frame(SEQN=c(labs2000$LAB25$SEQN, labs2002$L25_B$SEQN, labs2004$L25_C$SEQN, labs2006$CBC_D$SEQN, labs2008$CBC_E$SEQN,
                        labs2012$CBC_G$SEQN, labs2014$CBC_H$SEQN, labs2016$CBC_I$SEQN, labs2018$CBC_J$SEQN),
                 MCV=c(labs2000$LAB25$LBXMCVSI, labs2002$L25_B$LBXMCVSI, labs2004$L25_C$LBXMCVSI, labs2006$CBC_D$LBXMCVSI,labs2008$CBC_E$LBXMCVSI, 
                       labs2012$CBC_G$LBXMCVSI, labs2014$CBC_H$LBXMCVSI, labs2016$CBC_I$LBXMCVSI, labs2018$CBC_J$LBXMCVSI))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Amplitud de distribuci?n eritrocitaria ####

df2 = data.frame(SEQN=c(labs2000$LAB25$SEQN, labs2002$L25_B$SEQN, labs2004$L25_C$SEQN, labs2006$CBC_D$SEQN, labs2008$CBC_E$SEQN,
                        labs2012$CBC_G$SEQN, labs2014$CBC_H$SEQN, labs2016$CBC_I$SEQN, labs2018$CBC_J$SEQN),
                 RDW=c(labs2000$LAB25$LBXRDW, labs2002$L25_B$LBXRDW, labs2004$L25_C$LBXRDW, labs2006$CBC_D$LBXRDW, labs2008$CBC_E$LBXRDW,
                       labs2012$CBC_G$LBXRDW, labs2014$CBC_H$LBXRDW, labs2016$CBC_I$LBXRDW, labs2018$CBC_J$LBXRDW))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Fosfatasa alcalina ####

df2 = data.frame(SEQN=c(labs2000$LAB18$SEQN, labs2002$L40_B$SEQN, labs2004$L40_C$SEQN, labs2006$BIOPRO_D$SEQN, labs2008$BIOPRO_E$SEQN,
                        labs2012$BIOPRO_G$SEQN, labs2014$BIOPRO_H$SEQN, labs2016$BIOPRO_I$SEQN, labs2018$BIOPRO_J$SEQN),
                 ALP=c(labs2000$LAB18$LBXSAPSI, labs2002$L40_B$LBDSAPSI, labs2004$L40_C$LBXSAPSI, labs2006$BIOPRO_D$LBXSAPSI, labs2008$BIOPRO_E$LBXSAPSI,
                       labs2012$BIOPRO_G$LBXSAPSI, labs2014$BIOPRO_H$LBXSAPSI, labs2016$BIOPRO_I$LBXSAPSI, labs2018$BIOPRO_J$LBXSAPSI))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Recuento de gl?bulos blancos 1000cells/uL ####

df2 = data.frame(SEQN=c(labs2000$LAB25$SEQN, labs2002$L25_B$SEQN, labs2004$L25_C$SEQN, labs2006$CBC_D$SEQN, labs2008$CBC_E$SEQN,
                        labs2012$CBC_G$SEQN, labs2014$CBC_H$SEQN, labs2016$CBC_I$SEQN, labs2018$CBC_J$SEQN),
                 WBC=c(labs2000$LAB25$LBXWBCSI, labs2002$L25_B$LBXWBCSI, labs2004$L25_C$LBXWBCSI, labs2006$CBC_D$LBXWBCSI, labs2008$CBC_E$LBXWBCSI,
                       labs2012$CBC_G$LBXWBCSI, labs2014$CBC_H$LBXWBCSI, labs2016$CBC_I$LBXWBCSI, labs2018$CBC_J$LBXWBCSI))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Insulina ####

df2 = data.frame(SEQN=c(labs2000$LAB10AM$SEQN, labs2002$L10AM_B$SEQN, labs2004$L10AM_C$SEQN, labs2006$GLU_D$SEQN, 
                        labs2008$GLU_E$SEQN,labs2010$GLU_F$SEQN,labs2012$GLU_G$SEQN, labs2014$INS_H$SEQN, labs2016$INS_I$SEQN,labs2018$INS_J$SEQN),
                 Insulin=c(labs2000$LAB10AM$LBXINSI, labs2002$L10AM_B$LBXINSI, labs2004$L10AM_C$LBDINSI, labs2006$GLU_D$LBDINSI, 
                           labs2008$GLU_E$LBDINSI,labs2010$GLU_F$LBDINSI,labs2012$GLU_G$LBDINSI, labs2014$INS_H$LBDINSI, labs2016$INS_I$LBDINSI,labs2018$INS_J$LBDINSI))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Colesterol total ####

df2 = data.frame(SEQN=c(labs2000$Lab13$SEQN, labs2002$l13_b$SEQN, labs2004$l13_c$SEQN, labs2006$TCHOL_D$SEQN, 
                        labs2008$TCHOL_E$SEQN,labs2010$TCHOL_F$SEQN,labs2012$TCHOL_G$SEQN, labs2014$TCHOL_H$SEQN, labs2016$TCHOL_I$SEQN, labs2018$TCHOL_J$SEQN),
                 TCholesterol=c(labs2000$Lab13$LBXTC, labs2002$l13_b$LBXTC, labs2004$l13_c$LBXTC, labs2006$TCHOL_D$LBXTC, 
                                labs2008$TCHOL_E$LBXTC,labs2010$TCHOL_F$LBXTC,labs2012$TCHOL_G$LBXTC, labs2014$TCHOL_H$LBXTC, labs2016$TCHOL_I$LBXTC, labs2018$TCHOL_J$LBXTC))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")

# HDL ####

df2 = data.frame(SEQN=c(labs2000$Lab13$SEQN, labs2002$l13_b$SEQN, labs2004$l13_c$SEQN, labs2006$HDL_D$SEQN, 
                        labs2008$HDL_E$SEQN,labs2010$HDL_F$SEQN,labs2012$HDL_G$SEQN, labs2014$HDL_H$SEQN, labs2016$HDL_I$SEQN, labs2018$HDL_J$SEQN),
                 HDL=c(labs2000$Lab13$LBDHDL, labs2002$l13_b$LBDHDL, labs2004$l13_c$LBXHDD, labs2006$HDL_D$LBDHDD, 
                       labs2008$HDL_E$LBDHDD,labs2010$HDL_F$LBDHDD,labs2012$HDL_G$LBDHDD, labs2014$HDL_H$LBDHDD, labs2016$HDL_I$LBDHDD, labs2018$HDL_J$LBDHDD))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Triglic?ridos ####

df2 = data.frame(SEQN=c(labs2000$LAB18$SEQN, labs2002$L40_B$SEQN, labs2004$L40_C$SEQN, labs2006$BIOPRO_D$SEQN, 
                        labs2008$BIOPRO_E$SEQN,labs2010$BIOPRO_F$SEQN,labs2012$BIOPRO_G$SEQN, labs2014$BIOPRO_H$SEQN, labs2016$BIOPRO_I$SEQN, labs2018$BIOPRO_J$SEQN),
                 TAG=c(labs2000$LAB18$LBXSTR, labs2002$L40_B$LBXSTR, labs2004$L40_C$LBXSTR, labs2006$BIOPRO_D$LBXSTR, 
                       labs2008$BIOPRO_E$LBXSTR,labs2010$BIOPRO_F$LBXSTR,labs2012$BIOPRO_G$LBXSTR, labs2014$BIOPRO_H$LBXSTR, labs2016$BIOPRO_I$LBXSTR, labs2018$BIOPRO_J$LBXSTR))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")

# AST ####

df2 = data.frame(SEQN=c(labs2000$LAB18$SEQN, labs2002$L40_B$SEQN, labs2004$L40_C$SEQN, labs2006$BIOPRO_D$SEQN, 
                        labs2008$BIOPRO_E$SEQN,labs2010$BIOPRO_F$SEQN,labs2012$BIOPRO_G$SEQN, labs2014$BIOPRO_H$SEQN, labs2016$BIOPRO_I$SEQN, labs2018$BIOPRO_J$SEQN),
                 AST=c(labs2000$LAB18$LBXSASSI, labs2002$L40_B$LBXSASSI, labs2004$L40_C$LBXSASSI, labs2006$BIOPRO_D$LBXSASSI, 
                       labs2008$BIOPRO_E$LBXSASSI,labs2010$BIOPRO_F$LBXSASSI,labs2012$BIOPRO_G$LBXSASSI, labs2014$BIOPRO_H$LBXSASSI, labs2016$BIOPRO_I$LBXSASSI, labs2018$BIOPRO_J$LBXSASSI))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")

# ALT ####

df2 = data.frame(SEQN=c(labs2000$LAB18$SEQN, labs2002$L40_B$SEQN, labs2004$L40_C$SEQN, labs2006$BIOPRO_D$SEQN, 
                        labs2008$BIOPRO_E$SEQN,labs2010$BIOPRO_F$SEQN,labs2012$BIOPRO_G$SEQN, labs2014$BIOPRO_H$SEQN, labs2016$BIOPRO_I$SEQN, labs2018$BIOPRO_J$SEQN),
                 ALT=c(labs2000$LAB18$LBXSATSI, labs2002$L40_B$LBXSATSI, labs2004$L40_C$LBXSATSI, labs2006$BIOPRO_D$LBXSATSI, 
                       labs2008$BIOPRO_E$LBXSATSI,labs2010$BIOPRO_F$LBXSATSI,labs2012$BIOPRO_G$LBXSATSI, labs2014$BIOPRO_H$LBXSATSI, labs2016$BIOPRO_I$LBXSATSI, labs2018$BIOPRO_J$LBXSATSI))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")

# GGT #####

df2 = data.frame(SEQN=c(labs2000$LAB18$SEQN, labs2002$L40_B$SEQN, labs2004$L40_C$SEQN, labs2006$BIOPRO_D$SEQN, 
                        labs2008$BIOPRO_E$SEQN,labs2010$BIOPRO_F$SEQN,labs2012$BIOPRO_G$SEQN, labs2014$BIOPRO_H$SEQN, labs2016$BIOPRO_I$SEQN, labs2018$BIOPRO_J$SEQN),
                 GGT=c(labs2000$LAB18$LBXSGTSI, labs2002$L40_B$LBXSGTSI, labs2004$L40_C$LBXSGTSI, labs2006$BIOPRO_D$LBXSGTSI, 
                       labs2008$BIOPRO_E$LBXSGTSI,labs2010$BIOPRO_F$LBXSGTSI,labs2012$BIOPRO_G$LBXSGTSI, labs2014$BIOPRO_H$LBXSGTSI, labs2016$BIOPRO_I$LBXSGTSI, labs2018$BIOPRO_J$LBXSGTSI))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Hemoglobina glucosilada #### 

df2 = data.frame(SEQN=c(labs2000$LAB10$SEQN, labs2002$L10_B$SEQN, labs2004$L10_C$SEQN, labs2006$GHB_D$SEQN, 
                        labs2008$GHB_E$SEQN,labs2010$GHB_F$SEQN,labs2012$GHB_G$SEQN, labs2014$GHB_H$SEQN, labs2016$GHB_I$SEQN, labs2018$GHB_J$SEQN),
                 HbA1c=c(labs2000$LAB10$LBXGH, labs2002$L10_B$LBXGH, labs2004$L10_C$LBXGH, labs2006$GHB_D$LBXGH, 
                         labs2008$GHB_E$LBXGH,labs2010$GHB_F$LBXGH,labs2012$GHB_G$LBXGH, labs2014$GHB_H$LBXGH, labs2016$GHB_I$LBXGH, labs2018$GHB_J$LBXGH))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Peptido C ####

df2 = data.frame(SEQN=c(labs2000$LAB10AM$SEQN, labs2002$L10AM_B$SEQN, labs2004$L10AM_C$SEQN),
                 cpeptide=c(labs2000$LAB10AM$LBXCPSI, labs2002$L10AM_B$LBXCPSI, labs2004$L10AM_C$LBXCPSI))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")
# Tel?meros ####

df2 = data.frame(SEQN=c(labs2000$TELO_A$SEQN, labs2002$TELO_B$SEQN),
                 Telomere_mean=c(labs2000$TELO_A$TELOMEAN, labs2002$TELO_B$TELOMEAN))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")

df2 = data.frame(SEQN=c(labs2000$TELO_A$SEQN, labs2002$TELO_B$SEQN),
                 Telomere_sd=c(labs2000$TELO_A$TELOSTD, labs2002$TELO_B$TELOSTD))

df2$SEQN<-as.character(df2$SEQN)
NHANES<- NHANES %>% left_join(df2, by="SEQN")

NHANES_LABS <- NHANES

# C?lculo de PhenoAge ####

NHANES_LABS$Glucose <- NHANES_LABS$Glucose/18
NHANES_LABS$Cr <- NHANES_LABS$Cr*	88.42

NHANES_LABS <- NHANES_LABS %>%
  mutate(PhenoAge = 141.5 + ((log(-0.00553*log(1-(1-exp((-1.51714*exp(-19.907-0.0336*NHANES_LABS$Albumin+0.0095*NHANES_LABS$Cr+0.1953*NHANES_LABS$Glucose+
                                                                        0.0954*log(NHANES_LABS$CRP)-0.0120*NHANES_LABS$LymP+0.0268*NHANES_LABS$MCV+
                                                                        0.3306*NHANES_LABS$RDW+0.00188*NHANES_LABS$ALP+0.0554*NHANES_LABS$WBC+
                                                                        0.0804*NHANES_LABS$Age))/(0.0076927)))))))/(0.09165))

write.csv(NHANES_LABS, "nhanes_labs.csv")




# UNION DE TODAS LAS BASES ####

#Anthropometric
antro<-read.csv("antro_nhanes.csv")
#Laboratory
labs<-read.csv("nhanes_labs.csv")
#DXA
dxa<-read.csv("DXA_Antro.csv")
dxa2<-as.data.frame(apply(dxa,2,as.numeric)) #TODOS LOS SUJETOS ESTAN REPETIDOS 5 VECES POR EL METODO DE IMPUTACION

###Comorb
comorb<-read.csv("comorb_short_nhanes.csv")


#UNIR BASES
b.int<-merge(antro,labs,by="SEQN")
b.int2<-merge(b.int, comorb, by="SEQN")
NHANES<-(merge(b.int2,dxa2,by="SEQN"))
NHANES2<-(NHANES%>%group_by(SEQN)%>%summarise_all(mean,na.rm=T))

write.csv(NHANES2,"nhanes4.csv")
