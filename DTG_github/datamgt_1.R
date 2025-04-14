## ----------------------------- ##
## Data Management for DTG study ##
## ----------------------------- ##

## ----------------------  ##
##     Loading Dataset     ##
## ----------------------  ##
#rm(list=ls())

library(dplyr)
library(tidyverse)
library(Hmisc)
library(purrr)
library(lubridate)
library(rms)
library(gtsummary)
library(zoo) 
library(ggthemes)
library(ggrepel)
library(sandwich)
library(lmtest)
library(splines)



visit <- read.csv("visit.csv")

follow <- read.csv("follow.csv")

art <- read.csv("art.csv")

basic <- read.csv("basic.csv")

ce <- read.csv("ce.csv")

ce_cancer <- read.csv("ce_cancer.csv")

ce_tb <- read.csv("ce_tb.csv")

center <- read.csv("center.csv")

lab_cd4 <- read.csv("lab_cd4.csv")

lab_rna <- read.csv("lab_rna.csv")

## ----------------------  ##
##   Creating Variables    ##
## ----------------------  ##

#peru_follow <-  follow %>% filter(site == "peru")
#peru_visit <- visit %>% filter(site == "peru")





## Filtering Mexico and Peru & Argentina creating dataset for cross table for 1st filter ##
art_available <- art %>% filter(site != "mexico" & site != "peru" & site != "argentina")
dem_site <- art_available %>% distinct(patient_id, .keep_all = TRUE)

##     Country as Factor   ##
dem_site <- dem_site %>% mutate(study_site.factor = as.factor(if_else(site == "brazil", "Brazil",
                                                              if_else(site == "chile", "Chile",
                                                              if_else(site == "haiti", "Haiti","Honduras")))))

basic_gender <- basic %>% select(patient_id, male, birth_d,recart_y,enrol_d)                                                                                    
dem_site <- left_join(dem_site, basic_gender, by = "patient_id")
                                                                                                                                                                        
##     Gender as Factor  ##
dem_site <- dem_site %>% mutate(gender.factor = as.factor(ifelse(male == 1, "Male", 
                                                          ifelse(male == 0 ,"Female", "Transgender"))))

##    save dem_site  ##
save(dem_site, file="dem_site.Rdata")
label(dem_site$study_site.factor)   <- "Study site"
label(dem_site$gender.factor)   <- "Gender"

## removing Transgender Individuals##
dem_gender <- dem_site %>% filter(male != 8)

## save dem_gender ##
save(dem_gender, file="dem_gender.Rdata")
label(dem_gender$study_site.factor)   <- "Study site"
label(dem_gender$gender.factor)   <- "Gender"

## Filter for age: age >= 16, age calculated from birth_d and art_sd ##
art_available <- art_available %>% group_by(patient_id) %>% slice_min(art_sd)
dem_age <- dem_gender %>% select(patient_id, birth_d, study_site.factor, gender.factor, recart_y,enrol_d)
dem_age <- left_join(dem_age, art_available, by = "patient_id")
dem_age <- dem_age %>% mutate(age = difftime(as.Date(art_sd),as.Date(birth_d),units = "days")/365.25)  %>%
  filter(age >=16) 

## save dem_gender ##
save(dem_age, file="dem_age.Rdata")
label(dem_age$study_site.factor)   <- "Study site"
label(dem_age$gender.factor)   <- "Gender"

## Removing IDs with art start before enrollment ##
dem_recart <- dem_age %>% filter(recart_y == 0)

#save dem_recart
save(dem_recart, file="dem_recart.Rdata")
label(dem_recart$study_site.factor)   <- "Study site"
label(dem_recart$gender.factor)   <- "Gender"

# AIM 1
demdata <- dem_recart %>% mutate(AVAILABLE_DT = (if_else(site == "brazil", "2017-02-01",
                                              if_else(site == "chile", "2019-08-01",
                                              if_else(site == "haiti", "2018-11-01","2018-12-01")))))


demdata$naive <- ifelse(demdata$art_sd >= demdata$AVAILABLE_DT,1,0)
follow0 <- follow %>% select(patient_id,l_alive_d)
demdata <- left_join(demdata,follow0, by = "patient_id")
demdata <- demdata %>% filter(l_alive_d >= AVAILABLE_DT) %>% select(-l_alive_d)
demdata_1 <- demdata %>% filter(naive == 1)

# Filter for Aim 1 or ART Naive
save(demdata_1, file="demdata_1.Rdata")
label(demdata_1$study_site.factor)   <- "Study site"
label(demdata_1$gender.factor)   <- "Gender"

### Creating datasets and variables for AIM1 ###
demdata_1$age <- as.numeric(word(demdata_1$age, 1, sep = fixed(' ')))


## Variables for Table 1 for Aim 1  ##
## Gender clean ##
demdata_1$gender.factor <- ifelse(demdata_1$gender.factor == "Transgender", NA , 
                            ifelse(demdata_1$gender.factor == "Male", "Male", "Female"))
demdata_1$gender <- with(demdata_1,ifelse(gender.factor == "Male",0,1))

# Year at Baseline defined from art_sd, year of clinical enrollment from enrol_d #
demdata_1$baseline_year <- word(demdata_1$art_sd, 1, sep = fixed('-'))
demdata_1$enrol_year <- as.numeric(word(demdata_1$enrol_d, 1, sep = fixed('-')))

# DTG warning variable 
demdata_1 <- demdata_1 %>% mutate(DTG_warning = if_else(art_sd <= as.Date("2018-05-31"), "Pre-warning",
                                            if_else(as.Date("2018-06-01") <= art_sd & art_sd <= as.Date("2019-07-31"), "During warning", "After warning")))

# DTG yes or No
art1 <- art %>% select(patient_id,art_sd,art_id)
art1 <- art1 %>% group_by(patient_id) %>% slice_min(art_sd) 
art1$DTG <- ifelse(grepl("DTG",art1$art_id), 1, 0)
art1 <- art1 %>%  mutate(DTG.factor = ifelse(DTG == 1, "Yes","No")) %>% select(patient_id, DTG.factor,DTG)
demdata_1 <- left_join(demdata_1, art1, by = "patient_id")


## Last Observation from clinic visit ##
dem <- demdata_1 %>% select(patient_id,art_sd)
#create clean data from followup table
dem_follow <- follow %>% select(patient_id, l_alive_d)
dem_follow <- left_join(dem,dem_follow, by = 'patient_id')
dem_follow <- dem_follow %>% group_by(patient_id) %>%  filter(l_alive_d >= art_sd) %>% slice_max(l_alive_d)
dem_follow$clinic_year <- as.numeric(word(dem_follow$l_alive_d, 1, sep = fixed('-')))
dem_follow <- dem_follow %>% select(patient_id, l_alive_d, clinic_year)



## Create CD4 values ## 
dem_cd4 <- left_join(dem, lab_cd4, by = "patient_id")
dem_cd4$upper_limit <- as.Date(dem_cd4$art_sd) + 30
dem_cd4$lower_limit <- as.Date(dem_cd4$art_sd) %m-% months(6)
dem_cd4 <- dem_cd4 %>% filter(cd4_d   <= upper_limit & cd4_d >= lower_limit)
dem_cd4 <- dem_cd4 %>% mutate(difference = as.numeric(difftime(as.Date(art_sd),as.Date(cd4_d),units ="days")))
dem_cd4 <- dem_cd4 %>% group_by(patient_id) %>% slice_min(abs(difference))
dem_cd4 <- dem_cd4 %>% select(patient_id,cd4_v)

# create HIV RNA values ##
dem_rna <- left_join(dem, lab_rna, by = "patient_id")
dem_rna$upper_limit <- as.Date(dem_rna$art_sd) + 7
dem_rna$lower_limit <- as.Date(dem_rna$art_sd) %m-% months(6)
dem_rna <- dem_rna %>% filter(rna_d >= lower_limit)  
dem_rna <- dem_rna %>% filter(rna_d <= upper_limit)
dem_rna <- dem_rna %>% mutate(difference = as.numeric(difftime(as.Date(art_sd),as.Date(rna_d),units ="days")))
dem_rna_min <- dem_rna %>% filter(difference >= 0) %>%  group_by(patient_id) %>%  slice_min(difference) %>% select(patient_id,rna_v,difference,rna_l)
dem_rna_max <- dem_rna %>% filter(difference < 0) %>% group_by(patient_id) %>% slice_max(difference) %>% select(patient_id,rna_v,difference,rna_l) 
dem_rna <- full_join(dem_rna_max,dem_rna_min, by = c("patient_id", "difference", "rna_v"))
dem_rna <- dem_rna %>% group_by(patient_id) %>% slice_max(difference)
dem_rna$log_rna <- log10(abs(dem_rna$rna_v))
dem_rna$viral_supression <- ifelse((dem_rna$rna_v) >= 50, "Detectable", "Undetectable")
dem_rna <- dem_rna %>% distinct(patient_id, .keep_all = TRUE) %>% filter(viral_supression != "Undetectable") %>% select(patient_id,rna_v, log_rna, viral_supression,rna_l.x)


## History of Other AIDS Defining Illness ##
dem_ce <- left_join(dem, ce, by = "patient_id")
dem_ce$upper_limit <- as.Date(dem_ce$art_sd) + 30
dem_ce$lower_limit <- as.Date(dem_ce$art_sd) %m-% months(6)
dem_ce <- dem_ce %>% filter(ce_d  <= upper_limit & ce_d>= lower_limit)
dem_ce <- dem_ce %>% mutate(difference = as.numeric(difftime(as.Date(art_sd),as.Date(ce_d),units ="days")))
dem_ce <- dem_ce %>% mutate(ade_type= if_else(grepl('^ade', ce_id) & ce_id != "ade_tuberculosis" , 'Yes', 'No'))
dem_ce <- dem_ce %>% filter(ade_type == "Yes")
dem_ce <- dem_ce %>% group_by(patient_id) %>% slice_min(abs(difference))
dem_ce <- dem_ce %>% select(patient_id,ce_id,ade_type)

## History of TB ##
dem_tb <- left_join(dem, ce_tb, by = "patient_id")
dem_tb$upper_limit <- as.Date(dem_tb$art_sd) + 30
dem_tb$lower_limit <- as.Date(dem_tb$art_sd) %m-% months(6)
dem_tb <- dem_tb %>% filter(tbdiagnosis_d  <= upper_limit & tbdiagnosis_d >= lower_limit)
dem_tb <- dem_tb %>% mutate(difference = as.numeric(difftime(as.Date(art_sd),as.Date(tbdiagnosis_d),units ="days")))
dem_tb <- dem_tb %>% group_by(patient_id) %>% slice_min(abs(difference))
dem_tb <- dem_tb %>% mutate(tb = "Yes")
dem_tb <- dem_tb %>% select(patient_id,tb)

demdata_1 <- list(demdata_1, dem_rna, dem_cd4, dem_ce, dem_follow, dem_tb) %>% 
           reduce(left_join, by = "patient_id")
demdata_1 <- demdata_1 %>% distinct(patient_id, .keep_all = TRUE)

demdata_1 <- demdata_1 %>% mutate(tb = (if_else(is.na(tb), "No", "Yes")))
demdata_1 <- demdata_1 %>% mutate(ade_type = (if_else(is.na(ade_type), "No", "Yes")))

## TB yes or No ##
demdata_1 <- demdata_1 %>% mutate(tb = if_else(grepl('^ade_tuberculosis', ce_id) | tb == "Yes" , "Yes", 'No'))



##Creating ART Class category ##
demdata_1 <- demdata_1 %>% mutate(base = case_when(grepl('ABC', art_id) | grepl('ddI', art_id) | grepl('FTC', art_id) | grepl('3TC', art_id)|
                                                   grepl('d4T', art_id) | grepl('ddC', art_id) | grepl('AZT', art_id) | grepl('TDF', art_id) | 
                                                   grepl('TAF', art_id)  ~ 'NRTI', 
                                                   grepl('ETR', art_id) | grepl('DLV', art_id) | grepl('EFV', art_id) | grepl('NVP', art_id)|
                                                   grepl('RPV', art_id)  ~ 'NNRTI', 
                                                   grepl('RAL', art_id) | grepl('EVG', art_id) | grepl('DTG', art_id) | grepl('CAB', art_id)  ~ 'InSTI', 
                                                   grepl('APV', art_id) | grepl('ATZ', art_id) | grepl('DRV', art_id) | grepl('FPV', art_id)|
                                                   grepl('IDV', art_id) | grepl('LPV', art_id) | grepl('NFV', art_id) | grepl('RTV', art_id) | 
                                                   grepl('SQV', art_id) | grepl('TPV', art_id)  ~ 'PI',
                                                   grepl('ENF', art_id)   ~ 'FI' , 
                                                   grepl('VCV', art_id) | grepl('MVC', art_id)  ~ 'CA'))


# demdata_1 <- demdata_1 %>% mutate(ART_class = (if_else(DTG.factor == "Yes", "DTG",
#                                                if_else(nnrti !=0 | nnrti1 != 0 | nnrti2 !=0 & nrti != 0, "NNRTI",
#                                                if_else(nrti != 0 & pi != 0, "PI",
#                                                if_else(ii1 != 0 | ii2 !=0 & DTG.factor == "No" & nrti != 0 , "II",
#                                                if_else(nnrti ==0 & nnrti1 == 0 & nnrti2 ==0 & nrti == 0 & pi ==0 & ii1 ==0 & ii2 == 0, "Other", "None")))))))


demdata_1 <- demdata_1 %>% mutate(ART_class = ifelse(DTG.factor == "Yes", "DTG",
                                              ifelse(ii1 > 0 | ii2 > 0, "IINSTI-based",
                                              ifelse(nnrti > 0 & pi==0, "NNRTI-based",
                                              ifelse(pi > 0, "PI-based", "Other")))))


#Time till art Start 

demdata_1$duration_art <- as.numeric(difftime(as.Date(demdata_1$art_sd),as.Date(demdata_1$enrol_d),units="days"))

# CD4 category 
demdata_1$cd4.factor <- as.factor(ifelse(demdata_1$cd4_v > 350, ">350 cells/ mm3", 
                                  ifelse(demdata_1$cd4_v <= 350,"<=350 cells/ mm3", "Missing")))

demdata_1$sq_cd4 <- sqrt(demdata_1$cd4_v)

# variables Graph
demdata_1$DTG_year <- word(demdata_1$art_sd, 1, sep = fixed('-'))
demdata_1$DTG_month <- word(demdata_1$art_sd, 2, sep = fixed('-'))
demdata_1$month<- ifelse(demdata_1$DTG_month == "01" | demdata_1$DTG_month == "02" | demdata_1$DTG_month == "03"|
                              demdata_1$DTG_month == "04" | demdata_1$DTG_month == "05" | demdata_1$DTG_month == "06", "06", "12")

demdata_1 <- demdata_1 %>% mutate(year_semi = ifelse( month == "06",paste(DTG_year, month,30, sep = "-"),
                                              ifelse( month == "12",paste(DTG_year, month,31, sep = "-"), NA)))


demdata_1$year_semi <- as.Date(demdata_1$year_semi,format = "%Y-%m-%d")
                   

## -------------------------------------------  ##
##          un-adjusted RR for AIM 1 DTG        ##
## -------------------------------------------  ##




## creating levels and variables for interaction term(dichotomous age,DTG,gender)
#demdata_1 <- demdata_1 %>% select(age , gender.factor , study_site.factor ,  DTG_warning ,  sq_cd4 ,  log_rna  , ade_type , tb, DTG.factor)
#Dichotomous age, less then 50 vs Greater then 50 with "Greater then 50 as reference"
demdata_1 <- demdata_1 %>% select(-naive,-recart_y,-art_rs3,-art_rs4,-rna_l.x)
demdata_1$age_gt50 <- ifelse(demdata_1$age >= 50, 1,0)
demdata_1$age_gt50.factor <- factor(demdata_1$age_gt50,
                                    levels=c(1,0),
                                    labels=c("Greater then or equal to 50","Less then 50"))

#Dichotomous DTG, PRE and During DTG combined vs Post DTG, post DTG as reference
demdata_1$DTG_pre <- ifelse(demdata_1$DTG_warning == "Pre-warning" | demdata_1$DTG_warning == "During warning", 1, 0)
demdata_1$DTG_pre.factor <- factor(demdata_1$DTG_pre,
                                   levels=c(1,0),
                                   labels=c("Pre and during DTG","Post DTG"))

#Female vs Male, Male as reference
demdata_1$female <- ifelse(demdata_1$gender.factor == "Female", 1,0)
demdata_1$female.factor <- factor(demdata_1$female,
                                  levels=c(1,0),
                                  labels=c("Female","Male"))


#TB vs No TB, No TB as reference
demdata_1$tb_yes <- ifelse(demdata_1$tb == "Yes", 1,0)
demdata_1$tb_yes.factor <- factor(demdata_1$tb_yes,
                                  levels=c(1,0),
                                  labels=c("Yes","No"))

#ADE vs No ADE, No ADE as reference
demdata_1$ade_yes <- ifelse(demdata_1$ade_type == "Yes", 1,0)
demdata_1$ade_yes.factor <- factor(demdata_1$ade_yes,
                                   levels=c(1,0),
                                   labels=c("Yes","No"))

#CD4 Dichotomous greater then sqrt(350) vs less then sqrt(350), less then sqrt(350) as reference
demdata_1$cd4_gt350 <- ifelse(demdata_1$sq_cd4 >=  sqrt(350), 1,0)
demdata_1$cd4_gt350.factor <- factor(demdata_1$cd4_gt350,
                                     levels=c(1,0),
                                     labels=c("Greater then 350","Less then 350"))
#log RNA dichitimous, greater then 6500 vs less then 6500, with less then log(6500) as refrence 
demdata_1$rna_gt6500 <- ifelse(demdata_1$log_rna >=  log(6500,10), 1,0)
demdata_1$rna_gt6500.factor <- factor(demdata_1$rna_gt6500,
                                      levels=c(1,0),
                                      labels=c("Greater then 6500","Less then 6500"))

# creating new interaction term variable
demdata_1$age_gender_period <- with(demdata_1,ifelse(age_gt50.factor == "Greater then or equal to 50" & DTG_pre.factor == "Pre and during DTG" & gender.factor == "Male", "older male pre",
                                                     ifelse(age_gt50.factor == "Less then 50"                & DTG_pre.factor == "Pre and during DTG" & gender.factor == "Male", "younger male pre",
                                                            ifelse(age_gt50.factor == "Greater then or equal to 50" & DTG_pre.factor == "Post DTG"           & gender.factor == "Male", "older male post",
                                                                   ifelse(age_gt50.factor == "Greater then or equal to 50" & DTG_pre.factor == "Pre and during DTG" & gender.factor == "Female",  "older female pre",
                                                                          ifelse(age_gt50.factor == "Less then 50"                & DTG_pre.factor == "Post DTG"           & gender.factor == "Male","younger male post",
                                                                                 ifelse(age_gt50.factor == "Less then 50"                & DTG_pre.factor == "Pre and during DTG" & gender.factor == "Female", "younger female pre",
                                                                                        ifelse(age_gt50.factor == "Greater then or equal to 50" & DTG_pre.factor ==  "Post DTG"          & gender.factor == "Female", "older female post",
                                                                                               ifelse(age_gt50.factor == "Less then 50"                & DTG_pre.factor == "Post DTG"           & gender.factor == "Female","younger female post", "99")))))))))      


## Setting reference "young male post"
demdata_1$DTG <- ifelse(demdata_1$DTG.factor == "No", 0, 1) 
dd <- datadist(demdata_1)
options(datadist="dd")
dd$limits$DTG.factor[2] <- "No"
# Setting reference levels for regular model
demdata_1 <- within(demdata_1, DTG.factor <- relevel(factor(DTG.factor), ref = "No"))
demdata_1 <- within(demdata_1, gender.factor <- relevel(factor(gender.factor), ref = "Male"))
demdata_1 <- within(demdata_1, study_site.factor <- relevel(study_site.factor, ref = "Haiti"))
demdata_1 <- within(demdata_1, age_gender_period <- relevel(factor(age_gender_period), ref = "younger male post"))
demdata_1 <- within(demdata_1, age_gt50.factor <- relevel(factor(age_gt50.factor), ref = "Less then 50"))
demdata_1 <- within(demdata_1, DTG_pre.factor <- relevel(factor(DTG_pre.factor), ref = "Post DTG"))
demdata_1 <- within(demdata_1, female.factor <- relevel(factor(female.factor), ref = "Male"))
demdata_1 <- within(demdata_1, tb_yes.factor <- relevel(factor(tb_yes.factor), ref = "No"))
demdata_1 <- within(demdata_1, tb <- relevel(factor(tb), ref = "No"))
demdata_1 <- within(demdata_1, ade_yes <- relevel(factor(ade_yes.factor), ref = "No"))
demdata_1 <- within(demdata_1,cd4_gt350 <- relevel(factor(cd4_gt350.factor), ref = "Less then 350"))
demdata_1 <- within(demdata_1,rna_gt6500 <- relevel(factor(rna_gt6500.factor), ref = "Less then 6500"))

#setting reference levels for Glm model
dd$limits$DTG.factor[2] <- "No"
dd$limits$gender.factor[2] <- "Male"
dd$limits$female.factor[2] <- "Male"
dd$limits$female[2] <- 0
dd$limits$age_gt50.factor[2] <- "Less then 50"
dd$limits$tb_yes[2] <- "No"
dd$limits$tb_yes.factor[2] <- "No"
dd$limits$ade_yes[2] <- "No"
dd$limits$cd4_gt350[2] <- "Less then 350"
dd$limits$rna_gt6500[2] <- "Less then 6500"
dd$limits$study_site.factor[2] <- "Haiti"
dd$limits$age_gender_period[2] <- "younger male post"
dd$limits$DTG_pre.factor[2] <- "Post DTG"
dd$limits$DTG[2] <- 0


## ----------------------  ##
##         Model 1         ##
## ----------------------  ##

## Unadjusted RR for individual variables 
#female
m1_un_female <- glm(DTG ~   female.factor, data = demdata_1, family=poisson(link="log"))
m1_un_Female <- Glm(DTG ~   female.factor, data = demdata_1, family=poisson(link="log"))
summ_un_female <- as.data.frame(my_summary.rms(object=m1_un_Female, object2 = m1_un_female))
exposures <- " "
anova_un_female <- as.data.frame(anova_w_sandwich(object = m1_un_Female, objectglm = m1_un_female,coVars_woInt = c("female")))

#age_gt50
m1_un_age_gt50 <- glm(DTG ~   age_gt50.factor, data = demdata_1, family=poisson(link="log"))
m1_un_Age_gt50 <- Glm(DTG ~   age_gt50.factor, data = demdata_1, family=poisson(link="log"))
summ_un_age_gt50 <- as.data.frame(my_summary.rms(object=m1_un_Age_gt50, object2 = m1_un_age_gt50))
exposures <- " "
anova_un_age_gt50 <- as.data.frame(anova_w_sandwich(object = m1_un_Age_gt50, objectglm = m1_un_age_gt50,coVars_woInt = c("age_gt50")))

#TB
m1_un_tb <- glm(DTG ~   tb_yes.factor, data = demdata_1, family=poisson(link="log"))
m1_un_Tb <- Glm(DTG ~   tb_yes.factor, data = demdata_1, family=poisson(link="log"))
summ_un_tb <- as.data.frame(my_summary.rms(object=m1_un_Tb, object2 = m1_un_tb))
exposures <- " "
anova_un_tb <- as.data.frame(anova_w_sandwich(object = m1_un_Tb, objectglm = m1_un_tb,coVars_woInt = c("tb_yes")))

#ADE
# m1_un_ade <- glm(DTG ~   ade_yes.factor, data = demdata_1, family=poisson(link="log"))
# m1_un_Ade <- Glm(DTG ~   ade_yes.factor, data = demdata_1, family=poisson(link="log"))
# summ_un_ade <- as.data.frame(my_summary.rms(object=m1_un_Ade, object2 = m1_un_ade))
# exposures <- " "
# anova_un_ade <- as.data.frame(anova_w_sandwich(object = m1_un_Ade, objectglm = m1_un_ade,coVars_woInt = c("ade_yes")))

#cd4
m1_un_cd4 <- glm(DTG ~   cd4_gt350.factor, data = demdata_1, family=poisson(link="log"))
m1_un_Cd4 <- Glm(DTG ~   cd4_gt350.factor, data = demdata_1, family=poisson(link="log"))
summ_un_cd4 <- as.data.frame(my_summary.rms(object=m1_un_Cd4, object2 = m1_un_cd4))
exposures <- " "
anova_un_cd4 <- as.data.frame(anova_w_sandwich(object = m1_un_Cd4, objectglm = m1_un_cd4,coVars_woInt = c("cd4_gt350")))


#rna
m1_un_rna <- glm(DTG ~   rna_gt6500.factor, data = demdata_1, family=poisson(link="log"))
m1_un_Rna <- Glm(DTG ~   rna_gt6500.factor, data = demdata_1, family=poisson(link="log"))
summ_un_rna <- as.data.frame(my_summary.rms(object=m1_un_Rna, object2 = m1_un_rna))
exposures <- " "
anova_un_rna <- as.data.frame(anova_w_sandwich(object = m1_un_Rna, objectglm = m1_un_rna,coVars_woInt = c("rna_gt6500")))


#site
m1_un_site <- glm(DTG ~   study_site.factor, data = demdata_1, family=poisson(link="log"))
m1_un_Site <- Glm(DTG ~   study_site.factor, data = demdata_1, family=poisson(link="log"))
summ_un_site <- as.data.frame(my_summary.rms(object=m1_un_Site, object2 = m1_un_site))
exposures <- " "
anova_un_site <- as.data.frame(anova_w_sandwich(object = m1_un_Site, objectglm = m1_un_site,coVars_woInt = c("study_site.factor")))


#DTG
m1_un_dtg <- glm(DTG ~   DTG_pre.factor, data = demdata_1, family=poisson(link="log"))
m1_un_Dtg <- Glm(DTG ~   DTG_pre.factor, data = demdata_1, family=poisson(link="log"))
summ_un_dtg <- as.data.frame(my_summary.rms(object=m1_un_Dtg, object2 = m1_un_dtg))
exposures <- " "
anova_un_dtg <- as.data.frame(anova_w_sandwich(object = m1_un_Dtg, objectglm = m1_un_dtg,coVars_woInt = c("DTG_pre")))

##Unadjusted Model
variable <- c("Female", "Greater then or equal to 50 ref(Less then 50)", "Tb (ref = No)", "Site Brazil (ref= Haiti)", "Site Chile (ref= Haiti)", "Site Honduras (ref= Haiti)",
              "CD4 Greater then 350( ref = Less then 350)", "Greater then 6500 (ref =Less then 6500", "DTG (ref = Post DTG)")

rr <- c(round(exp(summ_un_female$Effect[1]), digits =3),round(exp(summ_un_age_gt50$Effect[1]), digits =3),round(exp(summ_un_tb$Effect[1]), digits =3),
        round(exp(summ_un_site$Effect[1]), digits =3),round(exp(summ_un_site$Effect[2]), digits =3),round(exp(summ_un_site$Effect[3]), digits =3),
        round(exp(summ_un_cd4$Effect[1]), digits =3),round(exp(summ_un_rna$Effect[1]), digits =3),round(exp(summ_un_dtg$Effect[1]), digits =3))                                                

ll <- c(round(exp(summ_un_female$`Lower 0.95`[1]), digits =3),round(exp(summ_un_age_gt50$`Lower 0.95`[1]), digits =3),round(exp(summ_un_tb$`Lower 0.95`[1]), digits =3),
        round(exp(summ_un_site$`Lower 0.95`[1]), digits =3),round(exp(summ_un_site$`Lower 0.95`[2]), digits =3),round(exp(summ_un_site$`Lower 0.95`[3]), digits =3),
        round(exp(summ_un_cd4$`Lower 0.95`[1]), digits =3),round(exp(summ_un_rna$`Lower 0.95`[1]), digits =3),round(exp(summ_un_dtg$`Lower 0.95`[1]), digits =3))                                                


ul <- c(round(exp(summ_un_female$`Upper 0.95`[1]), digits =3),round(exp(summ_un_age_gt50$`Upper 0.95`[1]), digits =3),round(exp(summ_un_tb$`Upper 0.95`[1]), digits =3),
        round(exp(summ_un_site$`Upper 0.95`[1]), digits =3),round(exp(summ_un_site$`Upper 0.95`[2]), digits =3),round(exp(summ_un_site$`Upper 0.95`[3]), digits =3),
        round(exp(summ_un_cd4$`Upper 0.95`[1]), digits =3),round(exp(summ_un_rna$`Upper 0.95`[1]), digits =3),round(exp(summ_un_dtg$`Upper 0.95`[1]), digits =3)) 

p <- c(anova_un_female$stats[1], anova_un_age_gt50$stats[1], anova_un_tb$stats[1], anova_un_site$stats[1], "","",anova_un_cd4$stats[1], anova_un_rna$stats[1],anova_un_dtg$stats[1])

unadjusted <- data.frame(variable,rr,ll,ul,p)


## ----------------------  ##
##         Model 2         ##
## ----------------------  ##
# dummy coding
m1_dummy_interaction <- glm(DTG ~   age_gender_period, data = demdata_1, family=poisson(link="log"))
m1_dummy_Interaction <- Glm(DTG ~   age_gender_period, data = demdata_1, family=poisson(link="log"))
summ_dummy_int <- as.data.frame(my_summary.rms(object=m1_dummy_Interaction, object2 = m1_dummy_interaction))
object2 = m1_dummy_interaction
var_dummy_int <- sandwich(object2)[-1, -1] # Removes Intercept row/column
exposures <- " "
anova_dummy_int <- as.data.frame(anova_w_sandwich(object = m1_dummy_Interaction, objectglm = m1_dummy_interaction,coVars_woInt = c("age_gender_period")))


variables <- c("older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(round(exp(summ_dummy_int$Effect[1]), digits =3),round(exp(summ_dummy_int$Effect[2]), digits =3),round(exp(summ_dummy_int$Effect[3]), digits =3)
        ,round(exp(summ_dummy_int$Effect[4]), digits =3),round(exp(summ_dummy_int$Effect[5]), digits =3),round(exp(summ_dummy_int$Effect[6]), digits =3)
        ,round(exp(summ_dummy_int$Effect[7]), digits =3))
ll <- c(round(exp(summ_dummy_int$`Lower 0.95`[1]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[2]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[3]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[4]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[5]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[6]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[7]), digits =3))

ul <- c(round(exp(summ_dummy_int$`Upper 0.95`[1]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[2]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[3]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[4]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[5]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[6]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[7]), digits =3))

z1 <- summ_dummy_int$Effect[1]/summ_dummy_int$S.E.[1]
z2 <- summ_dummy_int$Effect[2]/summ_dummy_int$S.E.[2]
z3 <- summ_dummy_int$Effect[3]/summ_dummy_int$S.E.[3]
z4 <- summ_dummy_int$Effect[4]/summ_dummy_int$S.E.[4]
z5 <- summ_dummy_int$Effect[5]/summ_dummy_int$S.E.[5]
z6 <- summ_dummy_int$Effect[6]/summ_dummy_int$S.E.[6]
z7 <- summ_dummy_int$Effect[7]/summ_dummy_int$S.E.[7]


p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3))



dummy_interaction <-  data.frame(variables,rr,ll,ul,p)


## ----------------------  ##
##         Model 3         ##
## ----------------------  ##

## Interaction terms model
m1_interaction <- glm(DTG ~    age_gt50.factor*DTG_pre.factor*gender.factor , data= demdata_1, family=poisson(link="log"))
summ_int <- summary(m1_interaction)
var <- sandwich(m1_interaction)



#OMP
ind_test = 2
test_beta = summ_int$coefficients[ind_test,1]
OMP_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
OMP_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(var[ind_test,ind_test])



#OFP
ind_test = c(2,4,6)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(2,3,4,5,6,7,8)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(2,3,5)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 4
test_beta = summ_int$coefficients[ind_test]
YFP_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
YFP_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(var[ind_test,ind_test])

#YFPr
ind_test = c(3,4,7)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = summ_int$coefficients[ind_test]
YMPr_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(var[ind_test,ind_test])


## young female pre vs young female post
ind_test = c(3,7)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


interaction1 <-  data.frame(variable, rr,ll,ul,p)



## -------------------------------------------  ##
##           Adjusted RR for AIM 2 DTG          ##
## -------------------------------------------  ##




#Dummy interaction term
mod_dummy_a <- glm(DTG ~    study_site.factor + tb_yes.factor + age_gender_period , data= demdata_1, family=poisson(link="log"))
mod_Dummy_a <- Glm(DTG ~  study_site.factor  + tb_yes.factor  + age_gender_period , data= demdata_1, family=poisson(link="log"))
summ_dummy_a <- as.data.frame(my_summary.rms(object = mod_Dummy_a, object2 = mod_dummy_a))
exposures <- " "
anova_cat_a <- as.data.frame(anova_w_sandwich(object = mod_Dummy_a, objectglm = mod_dummy_a,coVars_woInt = c(" tb_yes.factor "," study_site.factor  ","age_gender_period")))

# Table for adjusted model with dummy variables
variable <- c("Brazil","Chile","Honduras","Tb (ref =no)","older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")

rr <- c(round(exp(summ_dummy_a$Effect[1]), digits =3),round(exp(summ_dummy_a$Effect[2]), digits =3),round(exp(summ_dummy_a$Effect[3]), digits =3)
        ,round(exp(summ_dummy_a$Effect[4]), digits =3),round(exp(summ_dummy_a$Effect[5]), digits =3),round(exp(summ_dummy_a$Effect[6]), digits =3)
        ,round(exp(summ_dummy_a$Effect[7]), digits =3),round(exp(summ_dummy_a$Effect[8]), digits =3),round(exp(summ_dummy_a$Effect[9]), digits =3)
        ,round(exp(summ_dummy_a$Effect[10]), digits =3),round(exp(summ_dummy_a$Effect[11]), digits =3))

ll <- c(round(exp(summ_dummy_a$`Lower 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[8]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[9]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[10]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[11]), digits =3))

ul <- c(round(exp(summ_dummy_a$`Upper 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[8]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[9]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[10]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[11]), digits =3))                                                                                                                                                                                   

z1 <- summ_dummy_a$Effect[1]/summ_dummy_a$S.E.[1]
z2 <- summ_dummy_a$Effect[2]/summ_dummy_a$S.E.[2]
z3 <- summ_dummy_a$Effect[3]/summ_dummy_a$S.E.[3]
z4 <- summ_dummy_a$Effect[4]/summ_dummy_a$S.E.[4]
z5 <- summ_dummy_a$Effect[5]/summ_dummy_a$S.E.[5]
z6 <- summ_dummy_a$Effect[6]/summ_dummy_a$S.E.[6]
z7 <- summ_dummy_a$Effect[7]/summ_dummy_a$S.E.[7]
z8 <- summ_dummy_a$Effect[8]/summ_dummy_a$S.E.[8]
z9 <- summ_dummy_a$Effect[9]/summ_dummy_a$S.E.[9]
z10 <- summ_dummy_a$Effect[10]/summ_dummy_a$S.E.[10]
z11<- summ_dummy_a$Effect[11]/summ_dummy_a$S.E.[11]



p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits = 3),round((1 - pnorm(abs(z9))) * 2,digits =3),round((1 - pnorm(abs(z10))) * 2,digits=3),
       round((1 - pnorm(abs(z11))) * 2,digits =3))

interaction_dummy_a <- data.frame(variable,rr,ll,ul,p)


#Categorical interaction term
mod_cat <- glm(DTG ~   gender.factor + age_gt50.factor + tb_yes.factor + study_site.factor   + DTG_pre.factor + age_gt50.factor*DTG_pre.factor*gender.factor , data= demdata_1, family=poisson(link="log"))
summ_overall_a <- summary(mod_cat)
var_overall_a <- sandwich(mod_cat)


## Reporting overall Model with categorical interaction term


#YMP
ind = 3
beta = summ_overall_a$coefficients[ind,1]
YMP_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
YMP_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
YMP_point = exp(beta*(-1))
z1 <- beta/sqrt(var_overall_a[ind,ind])



#OFP
ind = c(2,3,10)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(2,3,8,9,10,11,12)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(3,8,9)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = summ_overall_a$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
YFP_point = exp(beta)
z5 <- beta/sqrt(var_overall_a[ind,ind])

#YFPr
ind = c(2,8,11)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 8
beta = summ_overall_a$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
YMPr_point = exp(beta)
z7 <- beta/sqrt(var_overall_a[ind,ind])


## young female pre vs young female post
ind = c(8,11)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(3,10)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta =  sum((-1)*summ_overall_a$coefficients[ind])
YFP_OFP_UB = exp(beta + 1.96*sqrt(se))
YFP_OFP_LB = exp(beta - 1.96*sqrt(se))
YFP_OFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,10)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,11)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,10,11,12)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)




#YFPr OFPr
ind = c(3,9,10,12)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
OFPr_YFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_YFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_YFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = summ_overall_a$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(var_overall_a[ind,ind])

ind = 5
beta = summ_overall_a$coefficients[ind]
bra_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
bra_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
bra_point = exp(beta)
bra <- beta/sqrt(var_overall_a[ind,ind])

ind = 6
beta = summ_overall_a$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(var_overall_a[ind,ind])

ind = 7
beta = summ_overall_a$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(var_overall_a[ind,ind])


#YFPr OFPr
ind = c(3,9)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
YMPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
YMPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
YMPr_OMPr_point = exp(beta)
z14 <- beta/sqrt(se)

ind = c(8,9,11,12)
var = as.numeric(var_overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_OFP_UB = exp(beta + 1.96*sqrt(se))
OFPr_OFP_LB = exp(beta - 1.96*sqrt(se))
OFPr_OFP_point = exp(beta)
z15 <- beta/sqrt(se)

ind = c(8,9)
var = as.numeric(var_overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OMPr_OMP_UB = exp(beta + 1.96*sqrt(se))
OMPr_OMP_LB = exp(beta - 1.96*sqrt(se))
OMPr_OMP_point = exp(beta)
z16 <- beta/sqrt(se)



var_overall <- c("Brazil","Chile","Honduras",
                 "TB (ref=no)",
                 "Younger male post DTG warning (ref = Older male post)", "older female postDTG warning (ref = younger male post)",
                 "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
                 "young female pre : young female post DTG warning", "younger Female post DTG warning : older female post DTG warning",
                 "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                 "older female pre  & during DTG warning: older male pre  & during DTG warning",
                 "younger female pre & during DTG warning : older female pre & during DTG warning","Younger Male Pre & during DTG warning:
                 Older Male pre & during warning", "Older female Pre: Older female Post", "Older Male Pre : Older Male Post")
rr <- c(bra_point,chile_point,hon_point,tb_point,YMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,YFP_OFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,OFPr_YFPr_point,YMPr_OMPr_point,OFPr_OFP_point,OMPr_OMP_point)

ll <- c(bra_LB,chile_LB,hon_LB,tb_LB,YMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,YFP_OFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, OFPr_YFPr_LB,YMPr_OMPr_LB,OFPr_OFP_LB,OMPr_OMP_LB)
ul <- c(bra_UB,chile_UB,hon_UB,tb_UB,YMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,YFP_OFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, OFPr_YFPr_UB,YMPr_OMPr_UB,OFPr_OFP_UB,OMPr_OMP_UB)

p <- c(round((1 - pnorm(abs(bra))) * 2,digits = 3),
       round((1 - pnorm(abs(chile))) * 2,digits = 3),
       round((1 - pnorm(abs(hon))) * 2,digits = 3),
       round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(z14))) * 2,digits=3),
       round((1 - pnorm(abs(z15))) * 2,digits=3),
       round((1 - pnorm(abs(z16))) * 2,digits=3))


interaction_a <-  data.frame(var_overall, rr,ll,ul,p)



# Inclusion for AIM 3 

demdata_1$l_alive_d <- as.Date(demdata_1$l_alive_d)
demdata_1$l_alive_d_n <- as.numeric(demdata_1$l_alive_d)
quantile(demdata_1$l_alive_d_n[demdata_1$study_site.factor == "Brazil"],c(0.9,0.95, 0.97))
demdata_1$cd <- NULL
demdata_1$cd[demdata_1$study_site.factor == "Brazil"] <- as.Date(18988)


quantile(demdata_1$l_alive_d_n[demdata_1$study_site.factor == "Honduras"],c(0.9,0.95, 0.97))
demdata_1$cd[demdata_1$study_site.factor == "Honduras"] <- as.Date(19006.4)

quantile(demdata_1$l_alive_d_n[demdata_1$study_site.factor == "Chile"],c(0.9,0.95, 0.97))
demdata_1$cd[demdata_1$study_site.factor == "Chile"] <- as.Date(19266.5)

quantile(demdata_1$l_alive_d_n[demdata_1$study_site.factor == "Haiti"],c(0.9,0.95, 0.97))
demdata_1$cd[demdata_1$study_site.factor == "Haiti"] <- as.Date(19069)


demdata_1$include3 <- ifelse(as.Date(demdata_1$art_sd) + 365 < demdata_1$cd, 1,0)



demdata_3 <- demdata_1 %>% filter(include3 == 1) %>% select(-rna_v)

test <- left_join(demdata_3, lab_rna , by = "patient_id")
test$undetectable0 <- ifelse(test$rna_v < 50, 1,0)
test$lower <- as.Date(test$art_sd) + 30
test$upper <- as.Date(test$art_sd) + 365
test <- test %>% filter(as.Date(rna_d) > lower & as.Date(rna_d) < upper)
test$undetectable1 <- ifelse(test$undetectable0 == 1,1,0)
test <- test %>% group_by(patient_id) %>% mutate(count = sum(undetectable1))
test$undetectable2 <- ifelse(test$count >= 1,1,0)
test <- test %>% distinct(patient_id,.keep_all = TRUE) %>% select(patient_id,undetectable2)

demdata_3 <- merge(demdata_3, test, by = "patient_id",all.x=TRUE)
demdata_3$undetectable3 <- with(demdata_3,ifelse(is.na(undetectable2), "Missing",
                                                 ifelse(undetectable2 == 1, "Undetectable","Detectable")))


demdata_3$undetectable4 <- with(demdata_3,ifelse(undetectable3 == "Missing" | undetectable3 == "Detectable" , 0,1))

## -------------------------------------------  ##
##           Un-adjusted RR for AIM 3            ##
## -------------------------------------------  ##

## Setting reference "young male post"
demdata_3$DTG <- ifelse(demdata_3$DTG.factor == "No", 0, 1) 
dd <- datadist(demdata_3)
options(datadist="dd")
dd$limits$DTG.factor[2] <- "No"
# Setting reference levels for regular model
demdata_3 <- within(demdata_3, DTG.factor <- relevel(factor(DTG.factor), ref = "No"))
demdata_3 <- within(demdata_3, gender.factor <- relevel(factor(gender.factor), ref = "Male"))
demdata_3 <- within(demdata_3, age_gender_period <- relevel(factor(age_gender_period), ref = "younger male post"))
demdata_3 <- within(demdata_3, age_gt50.factor <- relevel(factor(age_gt50.factor), ref = "Less then 50"))
demdata_3 <- within(demdata_3, DTG_pre.factor <- relevel(factor(DTG_pre.factor), ref = "Post DTG"))
demdata_3 <- within(demdata_3, female.factor <- relevel(factor(female.factor), ref = "Male"))
demdata_3 <- within(demdata_3, tb_yes.factor <- relevel(factor(tb_yes.factor), ref = "No"))
demdata_3 <- within(demdata_3, tb <- relevel(factor(tb), ref = "No"))
demdata_3 <- within(demdata_3, ade_yes <- relevel(factor(ade_yes.factor), ref = "No"))
demdata_3 <- within(demdata_3,cd4_gt350 <- relevel(factor(cd4_gt350.factor), ref = "Less then 350"))
demdata_3 <- within(demdata_3,rna_gt6500 <- relevel(factor(rna_gt6500.factor), ref = "Less then 6500"))

#setting reference levels for Glm model
dd$limits$DTG.factor[2] <- "No"
dd$limits$gender.factor[2] <- "Male"
dd$limits$female.factor[2] <- "Male"
dd$limits$female[2] <- 0
dd$limits$age_gt50.factor[2] <- "Less then 50"
dd$limits$tb_yes[2] <- "No"
dd$limits$tb_yes.factor[2] <- "No"
dd$limits$ade_yes[2] <- "No"
dd$limits$cd4_gt350[2] <- "Less then 350"
dd$limits$rna_gt6500[2] <- "Less then 6500"
dd$limits$age_gender_period[2] <- "younger male post"
dd$limits$DTG_pre.factor[2] <- "Post DTG"
dd$limits$DTG[2] <- 0
aim3_uni <- data.frame(Variable= c(
  c("Age (ref = Less then 50) ", "Greater then 50"),
  c("Age (ref = 35) ", 20,25,30,40,45,50,55,60),
  c("Site (ref= Haiti)", "Brazil", "Chile", "Honduras"),
  c("TB (ref = No)", "TB: Yes"),
  c("Gender (ref = Male)", "Females"),
  c("DTG (ref = no)", " DTG : Yes"),
  c("Timing (ref = Post DTG warning)", "Pre & during DTG warning")),
  RR= NA, Lower=NA, Upper=NA,  p=NA)



#age as categorical
m1 <- glm(undetectable4 ~ age_gt50.factor, demdata_3, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
age50_beta = s1$coefficients[ind_atest]
age50_UB = exp(age50_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
age50_LB = exp(age50_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
age50_pointest = exp(age50_beta)
z1 <- age50_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni$RR[aim3_uni$Variable %in% "Greater then 50"] <- age50_pointest
aim3_uni$Lower[aim3_uni$Variable %in% "Greater then 50"] <- age50_LB
aim3_uni$Upper[aim3_uni$Variable %in% "Greater then 50"] <- age50_UB
aim3_uni$p[aim3_uni$Variable %in% "Greater then 50"] <- p1

#age as spline
ages <- round(min(demdata_3$age)):round(max(demdata_3$age))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")



ms1 <- glm(undetectable4 ~ ns(age,df =4), demdata_3, family= poisson(link="log"))
summ_spline  <- summary(ms1)
v_spline <- sandwich(ms1)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr <- function(myage,model){
  exp(model$Coefs[model$Variables %in% "ns(age, df = 4)1"]*
        (sp[rownames(sp) %in% myage, 1]-sp[rownames(sp) %in% 35, 1]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)2"]*
        (sp[rownames(sp) %in% myage, 2]- sp[rownames(sp) %in% 35, 2]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)3"]*
        (sp[rownames(sp) %in% myage, 3]- sp[rownames(sp) %in% 35, 3]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)4"]*
        (sp[rownames(sp) %in% myage, 4]- sp[rownames(sp) %in% 35, 4]))
  
}

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)


# my function to calculate variance and SE
spline_se_fun <- function(myage){
  
  tmp <- v_spline
  a <- sp$X1[rownames(sp) %in% myage]-sp$X1[rownames(sp) %in% 35]
  b <- sp$X2[rownames(sp) %in% myage]-sp$X2[rownames(sp) %in% 35]
  c <- sp$X3[rownames(sp) %in% myage]-sp$X3[rownames(sp) %in% 35]
  d <- sp$X4[rownames(sp) %in% myage]-sp$X4[rownames(sp) %in% 35]
  
  Var_x <- tmp["ns(age, df = 4)1", "ns(age, df = 4)1"]
  Var_y <- tmp["ns(age, df = 4)2", "ns(age, df = 4)2"]
  Var_z <- tmp["ns(age, df = 4)3", "ns(age, df = 4)3"]
  Var_w <- tmp["ns(age, df = 4)4", "ns(age, df = 4)4"]
  
  Cov_xy <- tmp["ns(age, df = 4)1", "ns(age, df = 4)2"]
  Cov_xz <- tmp["ns(age, df = 4)1", "ns(age, df = 4)3"]
  Cov_xw <- tmp["ns(age, df = 4)1", "ns(age, df = 4)4"]
  Cov_yz <- tmp["ns(age, df = 4)2", "ns(age, df = 4)3"]
  Cov_yw <- tmp["ns(age, df = 4)2", "ns(age, df = 4)4"]
  Cov_zw <- tmp["ns(age, df = 4)3", "ns(age, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI <- function(rr,SE) 
{ exp(log(rr) + c(-1,1)*1.96*SE)}


CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)


aim3_uni$RR[aim3_uni$Variable %in% 20] <- rr_age_20
aim3_uni$RR[aim3_uni$Variable %in% 25] <- rr_age_25
aim3_uni$RR[aim3_uni$Variable %in% 30] <- rr_age_30
aim3_uni$RR[aim3_uni$Variable %in% 40] <- rr_age_40
aim3_uni$RR[aim3_uni$Variable %in% 45] <- rr_age_45
aim3_uni$RR[aim3_uni$Variable %in% 50] <- rr_age_50
aim3_uni$RR[aim3_uni$Variable %in% 55] <- rr_age_55
aim3_uni$RR[aim3_uni$Variable %in% 60] <- rr_age_60

aim3_uni$Lower[aim3_uni$Variable %in% 20] <- CI_20[1]
aim3_uni$Lower[aim3_uni$Variable %in% 25] <- CI_25[1]
aim3_uni$Lower[aim3_uni$Variable %in% 30] <- CI_30[1]
aim3_uni$Lower[aim3_uni$Variable %in% 40] <- CI_40[1]
aim3_uni$Lower[aim3_uni$Variable %in% 45] <- CI_45[1]
aim3_uni$Lower[aim3_uni$Variable %in% 50] <- CI_50[1]
aim3_uni$Lower[aim3_uni$Variable %in% 55] <- CI_55[1]
aim3_uni$Lower[aim3_uni$Variable %in% 60] <- CI_60[1]

aim3_uni$Upper[aim3_uni$Variable %in% 20] <- CI_20[2]
aim3_uni$Upper[aim3_uni$Variable %in% 25] <- CI_25[2]
aim3_uni$Upper[aim3_uni$Variable %in% 30] <- CI_30[2]
aim3_uni$Upper[aim3_uni$Variable %in% 40] <- CI_40[2]
aim3_uni$Upper[aim3_uni$Variable %in% 45] <- CI_45[2]
aim3_uni$Upper[aim3_uni$Variable %in% 50] <- CI_50[2]
aim3_uni$Upper[aim3_uni$Variable %in% 55] <- CI_55[2]
aim3_uni$Upper[aim3_uni$Variable %in% 60] <- CI_60[2]

p <- function(rr,SE)
{
  round((1 - pnorm(abs(log(rr)/SE))) * 2,digits = 3)
}


aim3_uni$p[aim3_uni$Variable %in% 20] <- p(rr_age_20,SE_20)
aim3_uni$p[aim3_uni$Variable %in% 25] <- p(rr_age_25,SE_25)
aim3_uni$p[aim3_uni$Variable %in% 30] <- p(rr_age_30,SE_30)
aim3_uni$p[aim3_uni$Variable %in% 40] <- p(rr_age_40,SE_40)
aim3_uni$p[aim3_uni$Variable %in% 45] <- p(rr_age_45,SE_45)
aim3_uni$p[aim3_uni$Variable %in% 50] <- p(rr_age_50,SE_55)
aim3_uni$p[aim3_uni$Variable %in% 55] <- p(rr_age_55,SE_55)
aim3_uni$p[aim3_uni$Variable %in% 60] <- p(rr_age_60,SE_60)

# Study Site 
m_site <- glm(undetectable4 ~   study_site.factor, data = demdata_3, family=poisson(link="log"))
summ_site  <- summary(m_site)
vsite <- sandwich(m_site)




ind_atest = 2
brazil_beta = m_site$coefficients[ind_atest]
brazil_UB = exp(brazil_beta  + 1.96*sqrt(vsite[ind_atest,ind_atest]))
brazil_LB = exp(brazil_beta  - 1.96*sqrt(vsite[ind_atest,ind_atest]))
brazil_pointest = exp(brazil_beta )
z3 <- brazil_beta /sqrt(vsite[ind_atest,ind_atest])
p3 <- round((1 - pnorm(abs(z3))) * 2,digits = 3)

aim3_uni$RR[aim3_uni$Variable %in% "Brazil"] <- brazil_pointest
aim3_uni$Lower[aim3_uni$Variable %in% "Brazil" ]<- brazil_LB
aim3_uni$Upper[aim3_uni$Variable %in% "Brazil"] <- brazil_UB
aim3_uni$p[aim3_uni$Variable %in% "Brazil"] <- p3

ind_atest = 3
chile_beta = m_site$coefficients[ind_atest]
chile_UB = exp(chile_beta  + 1.96*sqrt(vsite[ind_atest,ind_atest]))
chile_LB = exp(chile_beta  - 1.96*sqrt(vsite[ind_atest,ind_atest]))
chile_pointest = exp(chile_beta )
z4 <- chile_beta /sqrt(vsite[ind_atest,ind_atest])
p4 <- round((1 - pnorm(abs(z4))) * 2,digits = 3)

aim3_uni$RR[aim3_uni$Variable %in% "Chile"] <- chile_pointest
aim3_uni$Lower[aim3_uni$Variable %in% "Chile" ]<- chile_LB
aim3_uni$Upper[aim3_uni$Variable %in% "Chile"] <- chile_UB
aim3_uni$p[aim3_uni$Variable %in% "Chile"] <- p4

ind_atest = 4
honduras_beta = m_site$coefficients[ind_atest]
honduras_UB = exp(honduras_beta  + 1.96*sqrt(vsite[ind_atest,ind_atest]))
honduras_LB = exp(honduras_beta  - 1.96*sqrt(vsite[ind_atest,ind_atest]))
honduras_pointest = exp(honduras_beta )
z5 <- honduras_beta /sqrt(vsite[ind_atest,ind_atest])
p5 <- round((1 - pnorm(abs(z5))) * 2,digits = 3)

aim3_uni$RR[aim3_uni$Variable %in% "Honduras"] <- honduras_pointest
aim3_uni$Lower[aim3_uni$Variable %in% "Honduras" ]<- honduras_LB
aim3_uni$Upper[aim3_uni$Variable %in% "Honduras"] <- honduras_UB
aim3_uni$p[aim3_uni$Variable %in% "Honduras"] <- p5

#Other
# demdata_3$ade_yes <- ifelse(demdata_3$ade_type == "Yes", 1,0)
# demdata_3$ade_yes <- factor(demdata_3$ade_yes,
#                             levels=c(1,0),
#                             labels=c("Yes","No"))
# m1 <- glm(undetectable4 ~ ade_yes, demdata_3, family= poisson(link="log"))
# s1  <- summary(m1)
# v1 <- sandwich(m1)
# 
# 
# ind_atest = 2
# ade_beta = s1$coefficients[ind_atest]
# ade_UB = exp(ade_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
# ade_LB = exp(ade_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
# ade_pointest = exp(ade_beta)
# z1 <- ade_beta/sqrt(v1[ind_atest,ind_atest])
# p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)
# 
# aim3_uni$RR[aim3_uni$Variable %in% "Other AIDS defining illness: Yes"] <- ade_pointest
# aim3_uni$Lower[aim3_uni$Variable %in% "Other AIDS defining illness: Yes"] <- ade_LB
# aim3_uni$Upper[aim3_uni$Variable %in% "Other AIDS defining illness: Yes"] <- ade_UB
# aim3_uni$p[aim3_uni$Variable %in% "Other AIDS defining illness: Yes"] <- p1
# 

#TB
demdata_3$tb_yes <- ifelse(demdata_3$tb == "Yes", 1,0)
demdata_3$tb_yes <- factor(demdata_3$tb_yes,
                           levels=c(1,0),
                           labels=c("Yes","No"))
m1 <- glm(undetectable4 ~ tb_yes, demdata_3, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
tb_beta = s1$coefficients[ind_atest]
tb_UB = exp(tb_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
tb_LB = exp(tb_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
tb_pointest = exp(tb_beta)
z1 <- tb_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni$RR[aim3_uni$Variable %in% "TB: Yes"] <- tb_pointest
aim3_uni$Lower[aim3_uni$Variable %in% "TB: Yes"] <- tb_LB
aim3_uni$Upper[aim3_uni$Variable %in% "TB: Yes"] <- tb_UB
aim3_uni$p[aim3_uni$Variable %in% "TB: Yes"] <- p1

m1 <- glm(undetectable4 ~ gender.factor, demdata_3, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
gender_beta = s1$coefficients[ind_atest]
gender_UB = exp(gender_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
gender_LB = exp(gender_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
gender_pointest = exp(gender_beta)
z1 <- gender_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni$RR[aim3_uni$Variable %in% "Females"] <- gender_pointest
aim3_uni$Lower[aim3_uni$Variable %in% "Females"] <- gender_LB
aim3_uni$Upper[aim3_uni$Variable %in% "Females"] <- gender_UB
aim3_uni$p[aim3_uni$Variable %in% "Females"] <- p1

demdata_3 <- within(demdata_3, DTG_pre.factor <- relevel(factor(DTG_pre.factor), ref = "Post DTG"))
m1 <- glm(undetectable4 ~ DTG, demdata_3, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
DTG_beta = s1$coefficients[ind_atest]
DTG_UB = exp(DTG_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
DTG_LB = exp(DTG_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
DTG_pointest = exp(DTG_beta)
z1 <- DTG_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni$RR[aim3_uni$Variable %in% " DTG : Yes"] <- DTG_pointest
aim3_uni$Lower[aim3_uni$Variable %in% " DTG : Yes"] <- DTG_LB
aim3_uni$Upper[aim3_uni$Variable %in% " DTG : Yes"] <- DTG_UB
aim3_uni$p[aim3_uni$Variable %in% " DTG : Yes"] <- p1

## Timing
demdata_3$DTG_pre <- ifelse(demdata_3$DTG_warning == "Pre-warning" | demdata_3$DTG_warning == "During warning", 1, 0)
demdata_3$DTG_pre <- factor(demdata_3$DTG_pre,
                            levels=c(1,0),
                            labels=c("Pre and during DTG", "Post DTG" ))

m1 <- glm(undetectable4 ~ DTG_pre, demdata_3, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
timing_beta = s1$coefficients[ind_atest]
timing_UB = exp(timing_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
timing_LB = exp(timing_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
timing_pointest = exp(timing_beta)
z1 <- timing_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni$RR[aim3_uni$Variable %in% "Pre & during DTG warning"] <- timing_pointest
aim3_uni$Lower[aim3_uni$Variable %in% "Pre & during DTG warning"] <- timing_LB
aim3_uni$Upper[aim3_uni$Variable %in% "Pre & during DTG warning"] <- timing_UB
aim3_uni$p[aim3_uni$Variable %in% "Pre & during DTG warning"] <- p1

demdata_3$age_gender_period <- with(demdata_3,ifelse(age_gt50 == 0 & DTG_pre == "Pre and during DTG" & gender.factor == "Male", "older male pre",
                                                     ifelse(age_gt50 == 1                & DTG_pre == "Pre and during DTG" & gender.factor == "Male", "younger male pre",
                                                            ifelse(age_gt50 == 0 & DTG_pre == "Post DTG"           & gender.factor == "Male", "older male post",
                                                                   ifelse(age_gt50 == 0 & DTG_pre == "Pre and during DTG" & gender.factor == "Female",  "older female pre",
                                                                          ifelse(age_gt50 == 1                & DTG_pre == "Post DTG" & gender.factor == "Male","younger male post",
                                                                                 ifelse(age_gt50 == 1               & DTG_pre == "Pre and during DTG" & gender.factor == "Female", "younger female pre",
                                                                                        ifelse(age_gt50 == 0 & DTG_pre ==  "Post DTG" &  gender.factor == "Female", "older female post",
                                                                                               ifelse(age_gt50 == 1                & DTG_pre == "Post DTG" & gender.factor == "Female","younger female post", "Missing")))))))))      



m1 <- glm(undetectable4 ~  age_gender_period , demdata_3, family= poisson(link="log"))
s1 <- summary(m1)
v1 <- sandwich(m1)

ind_atest = 2
ofp_beta = s1$coefficients[ind_atest]
ofp_UB = exp(ofp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ofp_LB = exp(ofp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ofp_pointest = exp(ofp_beta)
z1 <- ofp_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 3
ofpr_beta = s1$coefficients[ind_atest]
ofpr_UB = exp(ofpr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ofpr_LB = exp(ofpr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ofpr_pointest = exp(ofpr_beta)
z1 <- ofpr_beta/sqrt(v1[ind_atest,ind_atest])
p2 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 4
omp_beta = s1$coefficients[ind_atest]
omp_UB = exp(omp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
omp_LB = exp(omp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
omp_pointest = exp(omp_beta)
z1 <- omp_beta/sqrt(v1[ind_atest,ind_atest])
p3 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 5
ompr_beta = s1$coefficients[ind_atest]
ompr_UB = exp(ompr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ompr_LB = exp(ompr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ompr_pointest = exp(ompr_beta)
z1 <- ompr_beta/sqrt(v1[ind_atest,ind_atest])
p4 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 6
yfp_beta = s1$coefficients[ind_atest]
yfp_UB = exp(yfp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
yfp_LB = exp(yfp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
yfp_pointest = exp(yfp_beta)
z1 <- yfp_beta/sqrt(v1[ind_atest,ind_atest])
p5 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)


ind_atest = 7
yfpr_beta = s1$coefficients[ind_atest]
yfpr_UB = exp(yfpr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
yfpr_LB = exp(yfpr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
yfpr_pointest = exp(yfpr_beta)
z1 <- yfpr_beta/sqrt(v1[ind_atest,ind_atest])
p6 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 8
ympr_beta = s1$coefficients[ind_atest]
ympr_UB = exp(ympr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ympr_LB = exp(ympr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ympr_pointest = exp(ympr_beta)
z1 <- ympr_beta/sqrt(v1[ind_atest,ind_atest])
p7 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

variables <- c("older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(round(ofp_pointest, digits =3),round(ofpr_pointest, digits =3),round(omp_pointest, digits =3)
        ,round(ompr_pointest, digits =3),round(yfp_pointest, digits =3),round(yfpr_pointest, digits =3)
        ,round(ympr_pointest, digits =3))
ll <- c(round(ofp_LB, digits =3),
        round(ofpr_LB, digits =3),
        round(omp_LB, digits =3),
        round(ompr_LB, digits =3),
        round(yfp_LB, digits =3),
        round(yfpr_LB, digits =3),
        round(ympr_LB, digits =3))

ul <- c(round(ofp_UB, digits =3),
        round(ofpr_UB, digits =3),
        round(omp_UB, digits =3),
        round(ompr_UB, digits =3),
        round(yfp_UB, digits =3),
        round(yfpr_UB, digits =3),
        round(ympr_UB, digits =3))




p <- c(p1,p2,p3,p4,p5,p6,p7)

aim3_int <-  data.frame(variables,rr,ll,ul,p)



i3 <- glm(undetectable4 ~   study_site.factor + tb_yes + age_gender_period , demdata_3, family= poisson(link="log"))
is3 <- summary(i3)
iv3 <- sandwich(i3)



ind = 2
beta = is3$coefficients[ind]
bra_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
bra_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
bra_point = exp(beta)
bra <- beta/sqrt(iv3[ind,ind])

ind = 3
beta = is3$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(iv3[ind,ind])

ind = 4
beta = is3$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(iv3[ind,ind])

# ind = 5
# beta = is3$coefficients[ind]
# ade_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
# ade_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
# ade_point = exp(beta)
# ade <- beta/sqrt(iv3[ind,ind])

ind = 5
beta = is3$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv3[ind,ind])

ind = 6
beta = is3$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv3[ind,ind])

ind = 7
beta = is3$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv3[ind,ind])

ind = 8
beta = is3$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv3[ind,ind])

ind = 9
beta = is3$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv3[ind,ind])

ind = 10
beta = is3$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv3[ind,ind])

ind = 11
beta = is3$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv3[ind,ind])

ind = 12
beta = is3$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv3[ind,ind])

variables <- c("Brazil","Chile","Honduras","TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(bra_point,chile_point,hon_point,tb_point,
        OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point)
ll <- c(bra_LB,chile_LB,hon_LB,tb_LB,
        OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB)
ul <- c(bra_UB,chile_UB,hon_UB,tb_UB,
        OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB)
p <- c(round((1 - pnorm(abs(bra))) * 2,digits = 3),round((1 - pnorm(abs(chile))) * 2,digits = 3),round((1 - pnorm(abs(hon))) * 2,digits = 3),round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3))



aim3_dummy <-  data.frame(variables, rr,ll,ul,p)

# Overall model with DTG
i4 <- glm(undetectable4 ~   study_site.factor + tb_yes + age_gender_period  + DTG, demdata_3, family= poisson(link="log"))
is4 <- summary(i4)
iv4<- sandwich(i4)


ind = 2
beta = is4$coefficients[ind]
bra_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
bra_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
bra_point = exp(beta)
bra <- beta/sqrt(iv4[ind,ind])

ind = 3
beta = is4$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(iv4[ind,ind])

ind = 4
beta = is4$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(iv4[ind,ind])

# ind = 5
# beta = is4$coefficients[ind]
# ade_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
# ade_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
# ade_point = exp(beta)
# ade <- beta/sqrt(iv4[ind,ind])

ind = 5
beta = is4$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv4[ind,ind])

ind = 6
beta = is4$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv4[ind,ind])

ind = 7
beta = is4$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv4[ind,ind])

ind = 8
beta = is4$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv4[ind,ind])

ind = 9
beta = is4$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv4[ind,ind])

ind = 10
beta = is4$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv4[ind,ind])

ind = 11
beta = is4$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv4[ind,ind])

ind = 12
beta = is4$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv4[ind,ind])

ind = 13
beta = is4$coefficients[ind]
DTG_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
DTG_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
DTG_point = exp(beta)
z_DTG <- beta/sqrt(iv4[ind,ind])


variables <- c("Brazil","Chile","Honduras","TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)", "DTG : Yes/No")
rr <- c(bra_point,chile_point,hon_point,tb_point,
        OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point, DTG_point)
ll <- c(bra_LB,chile_LB,hon_LB,tb_LB,
        OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB, DTG_LB)
ul <- c(bra_UB,chile_UB,hon_UB,tb_UB,
        OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB, DTG_UB)
p <- c(round((1 - pnorm(abs(bra))) * 2,digits = 3),round((1 - pnorm(abs(chile))) * 2,digits = 3),round((1 - pnorm(abs(hon))) * 2,digits = 3),round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3),round((1 - pnorm(abs(z_DTG))) * 2,digits=3))


aim3_dummy_DTG <-  data.frame(variables, rr,ll,ul,p)

# UNdetectable2
i3 <- glm(undetectable2 ~   study_site.factor + tb_yes + age_gender_period , demdata_3, family= poisson(link="log"))
is3 <- summary(i3)
iv3 <- sandwich(i3)



ind = 2
beta = is3$coefficients[ind]
bra_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
bra_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
bra_point = exp(beta)
bra <- beta/sqrt(iv3[ind,ind])

ind = 3
beta = is3$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(iv3[ind,ind])

ind = 4
beta = is3$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(iv3[ind,ind])

# ind = 5
# beta = is3$coefficients[ind]
# ade_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
# ade_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
# ade_point = exp(beta)
# ade <- beta/sqrt(iv3[ind,ind])

ind = 5
beta = is3$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv3[ind,ind])

ind = 6
beta = is3$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv3[ind,ind])

ind = 7
beta = is3$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv3[ind,ind])

ind = 8
beta = is3$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv3[ind,ind])

ind = 9
beta = is3$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv3[ind,ind])

ind = 10
beta = is3$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv3[ind,ind])

ind = 11
beta = is3$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv3[ind,ind])

ind = 12
beta = is3$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv3[ind,ind])

variables <- c("Brazil","Chile","Honduras","TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(bra_point,chile_point,hon_point,tb_point,
        OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point)
ll <- c(bra_LB,chile_LB,hon_LB,tb_LB,
        OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB)
ul <- c(bra_UB,chile_UB,hon_UB,tb_UB,
        OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB)
p <- c(round((1 - pnorm(abs(bra))) * 2,digits = 3),round((1 - pnorm(abs(chile))) * 2,digits = 3),round((1 - pnorm(abs(hon))) * 2,digits = 3),round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3))



undetectable2_overall <-  data.frame(variables, rr,ll,ul,p)

# Overall model with DTG
i4 <- glm(undetectable2 ~   study_site.factor + tb_yes + age_gender_period  + DTG, demdata_3, family= poisson(link="log"))
is4 <- summary(i4)
iv4<- sandwich(i4)


ind = 2
beta = is4$coefficients[ind]
bra_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
bra_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
bra_point = exp(beta)
bra <- beta/sqrt(iv4[ind,ind])

ind = 3
beta = is4$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(iv4[ind,ind])

ind = 4
beta = is4$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(iv4[ind,ind])

# ind = 5
# beta = is4$coefficients[ind]
# ade_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
# ade_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
# ade_point = exp(beta)
# ade <- beta/sqrt(iv4[ind,ind])

ind = 5
beta = is4$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv4[ind,ind])

ind = 6
beta = is4$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv4[ind,ind])

ind = 7
beta = is4$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv4[ind,ind])

ind = 8
beta = is4$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv4[ind,ind])

ind = 9
beta = is4$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv4[ind,ind])

ind = 10
beta = is4$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv4[ind,ind])

ind = 11
beta = is4$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv4[ind,ind])

ind = 12
beta = is4$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv4[ind,ind])

ind = 13
beta = is4$coefficients[ind]
DTG_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
DTG_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
DTG_point = exp(beta)
z_DTG <- beta/sqrt(iv4[ind,ind])


variables <- c("Brazil","Chile","Honduras","TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)", "DTG : Yes/No")
rr <- c(bra_point,chile_point,hon_point,tb_point,
        OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point, DTG_point)
ll <- c(bra_LB,chile_LB,hon_LB,tb_LB,
        OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB, DTG_LB)
ul <- c(bra_UB,chile_UB,hon_UB,tb_UB,
        OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB, DTG_UB)
p <- c(round((1 - pnorm(abs(bra))) * 2,digits = 3),round((1 - pnorm(abs(chile))) * 2,digits = 3),round((1 - pnorm(abs(hon))) * 2,digits = 3),round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3),round((1 - pnorm(abs(z_DTG))) * 2,digits=3))


undetectable2_DTG <-  data.frame(variables, rr,ll,ul,p)

#levels for undetectable 2
demdata_3$undetectable2.factor <- factor(demdata_3$undetectable2,
                                         levels=c(1,0),
                                         labels=c("Undetectable","Detectable"))



## ----------------------  ##
##          Labels         ##
## ----------------------  ##
label(demdata_1$age)   <- "Age at baseline"
label(demdata_1$age_gt50.factor) <- "Age(category) at baseline"
label(demdata_1$study_site.factor)   <- "Study site"
label(demdata_1$baseline_year) <- "Year of baseline"
label(demdata_1$enrol_year)   <- "Year of clinic enrollment"
label(demdata_1$DTG_warning) <- "DTG warning calendar era at baseline"
label(demdata_1$clinic_year) <- "Year of last clinic visit"
label(demdata_1$cd4_v)  <- "CD4 at baseline "
label(demdata_1$viral_supression) <-  "HIV-1 RNA at baseline categorized as detectable vs undetectable"
label(demdata_1$log_rna) <- "HIV-1 RNA at baseline (log10-transformed) "
label(demdata_1$ade_type) <- "History of other AIDS-Defining illness (not TB)"
label(demdata_1$DTG.factor) <- "DTG Given Yes/No"
label(demdata_1$tb) <- "TB Yes/No"
label(demdata_1$ART_class) <- "ART Class"
label(demdata_1$gender.factor) <- "Gender  "
label(demdata_1$duration_art)  <- "Time from site enrollment to ART initiation in days"
label(demdata_1$cd4.factor) <- "CD4 at baseline "
label(demdata_3$undetectable3) <- "Undetectable viral load"
label(demdata_3$undetectable2.factor) <- "Undetectable viral load(without missing viral load)"
label(demdata_3$age)   <- "Age at baseline"
label(demdata_3$study_site.factor)   <- "Study site"
label(demdata_3$baseline_year) <- "Year of baseline"
label(demdata_3$enrol_year)   <- "Year of clinic enrollment"
label(demdata_3$DTG_warning) <- "DTG warning calendar era at baseline"
label(demdata_3$clinic_year) <- "Year of last clinic visit"
label(demdata_3$cd4_v)  <- "CD4 at baseline "
label(demdata_3$viral_supression) <-  "HIV-1 RNA at baseline categorized as detectable vs undetectable"
label(demdata_3$log_rna) <- "HIV-1 RNA at baseline (log10-transformed) "
label(demdata_3$ade_type) <- "History of other AIDS-Defining illness (not TB)"
label(demdata_3$DTG.factor) <- "DTG Given Yes/No"
label(demdata_3$tb) <- "TB Yes/No"
label(demdata_3$ART_class) <- "ART Class"
label(demdata_3$gender.factor) <- "Gender  "
label(demdata_3$duration_art)  <- "Time from site enrollment to ART initiation in years"
label(demdata_3$cd4.factor) <- "CD4 at baseline "
label(demdata_3$age_gt50.factor) <- "Age(category) at baseline"



save(demdata_1, file="demdata_1.Rdata")
save(demdata_3, file="demdata_3.Rdata")
save(unadjusted,file ="unadjusted.Rdata")
save(interaction1,file ="interaction1.Rdata")
save(dummy_interaction,file ="dummy_interaction.Rdata")
#save(overall,file ="overall.Rdata")
save(interaction_a,file ="interaction_a.Rdata")
save(interaction_dummy_a,file ="interaction_dummy_a.Rdata")
save(aim3_uni,file ="aim3_uni.Rdata")
save(aim3_int,file ="aim3_int.Rdata")
save(aim3_dummy,file ="aim3_dummy.Rdata")
save(aim3_dummy_DTG,file ="aim3_dummy_DTG.Rdata")
save(undetectable2_overall,file ="undetectable2_overall.Rdata")
save(undetectable2_DTG,file ="undetectable2_DTG.Rdata")





## Different analysis for Haiti
demdata_h <- demdata_1 %>% filter(study_site.factor == "Haiti")
demdata_h <- demdata_h %>% select(-t20,-ccr5,-ce_id,-cd)
## Setting reference "young male post"
demdata_h$DTG <- ifelse(demdata_h$DTG.factor == "No", 0, 1) 
dd <- datadist(demdata_h)
options(datadist="dd")
dd$limits$DTG.factor[2] <- "No"
# Setting reference levels for regular model
demdata_h <- within(demdata_h, DTG.factor <- relevel(factor(DTG.factor), ref = "No"))
demdata_h <- within(demdata_h, gender.factor <- relevel(factor(gender.factor), ref = "Male"))
demdata_h <- within(demdata_h, age_gender_period <- relevel(factor(age_gender_period), ref = "younger male post"))
demdata_h <- within(demdata_h, age_gt50.factor <- relevel(factor(age_gt50.factor), ref = "Less then 50"))
demdata_h <- within(demdata_h, DTG_pre.factor <- relevel(factor(DTG_pre.factor), ref = "Post DTG"))
demdata_h <- within(demdata_h, female.factor <- relevel(factor(female.factor), ref = "Male"))
demdata_h <- within(demdata_h, tb_yes.factor <- relevel(factor(tb_yes.factor), ref = "No"))
demdata_h <- within(demdata_h, tb <- relevel(factor(tb), ref = "No"))
demdata_h <- within(demdata_h, ade_yes <- relevel(factor(ade_yes.factor), ref = "No"))
demdata_h <- within(demdata_h,cd4_gt350 <- relevel(factor(cd4_gt350.factor), ref = "Less then 350"))
demdata_h <- within(demdata_h,rna_gt6500 <- relevel(factor(rna_gt6500.factor), ref = "Less then 6500"))

#setting reference levels for Glm model
dd$limits$DTG.factor[2] <- "No"
dd$limits$gender.factor[2] <- "Male"
dd$limits$female.factor[2] <- "Male"
dd$limits$female[2] <- 0
dd$limits$age_gt50.factor[2] <- "Less then 50"
dd$limits$tb_yes[2] <- "No"
dd$limits$tb_yes.factor[2] <- "No"
dd$limits$ade_yes[2] <- "No"
dd$limits$cd4_gt350[2] <- "Less then 350"
dd$limits$rna_gt6500[2] <- "Less then 6500"
dd$limits$age_gender_period[2] <- "younger male post"
dd$limits$DTG_pre.factor[2] <- "Post DTG"
dd$limits$DTG[2] <- 0


## Unadjusted RR for individual variables 
#female
m1_un_female <- glm(DTG ~   female.factor, data = demdata_h, family=poisson(link="log"))
m1_un_Female <- Glm(DTG ~   female.factor, data = demdata_h, family=poisson(link="log"))
summ_un_female <- as.data.frame(my_summary.rms(object=m1_un_Female, object2 = m1_un_female))
exposures <- " "
anova_un_female <- as.data.frame(anova_w_sandwich(object = m1_un_Female, objectglm = m1_un_female,coVars_woInt = c("female")))

#age_gt50
m1_un_age_gt50 <- glm(DTG ~   age_gt50.factor, data = demdata_h, family=poisson(link="log"))
m1_un_Age_gt50 <- Glm(DTG ~   age_gt50.factor, data = demdata_h, family=poisson(link="log"))
summ_un_age_gt50 <- as.data.frame(my_summary.rms(object=m1_un_Age_gt50, object2 = m1_un_age_gt50))
exposures <- " "
anova_un_age_gt50 <- as.data.frame(anova_w_sandwich(object = m1_un_Age_gt50, objectglm = m1_un_age_gt50,coVars_woInt = c("age_gt50")))

#TB
m1_un_tb <- glm(DTG ~   tb_yes.factor, data = demdata_h, family=poisson(link="log"))
m1_un_Tb <- Glm(DTG ~   tb_yes.factor, data = demdata_h, family=poisson(link="log"))
summ_un_tb <- as.data.frame(my_summary.rms(object=m1_un_Tb, object2 = m1_un_tb))
exposures <- " "
anova_un_tb <- as.data.frame(anova_w_sandwich(object = m1_un_Tb, objectglm = m1_un_tb,coVars_woInt = c("tb_yes")))


#cd4
m1_un_cd4 <- glm(DTG ~   cd4_gt350.factor, data = demdata_h, family=poisson(link="log"))
m1_un_Cd4 <- Glm(DTG ~   cd4_gt350.factor, data = demdata_h, family=poisson(link="log"))
summ_un_cd4 <- as.data.frame(my_summary.rms(object=m1_un_Cd4, object2 = m1_un_cd4))
exposures <- " "
anova_un_cd4 <- as.data.frame(anova_w_sandwich(object = m1_un_Cd4, objectglm = m1_un_cd4,coVars_woInt = c("cd4_gt350")))


#rna
m1_un_rna <- glm(DTG ~   rna_gt6500.factor, data = demdata_h, family=poisson(link="log"))
m1_un_Rna <- Glm(DTG ~   rna_gt6500.factor, data = demdata_h, family=poisson(link="log"))
summ_un_rna <- as.data.frame(my_summary.rms(object=m1_un_Rna, object2 = m1_un_rna))
exposures <- " "
anova_un_rna <- as.data.frame(anova_w_sandwich(object = m1_un_Rna, objectglm = m1_un_rna,coVars_woInt = c("rna_gt6500")))



#DTG
m1_un_dtg <- glm(DTG ~   DTG_pre.factor, data = demdata_h, family=poisson(link="log"))
m1_un_Dtg <- Glm(DTG ~   DTG_pre.factor, data = demdata_h, family=poisson(link="log"))
summ_un_dtg <- as.data.frame(my_summary.rms(object=m1_un_Dtg, object2 = m1_un_dtg))
exposures <- " "
anova_un_dtg <- as.data.frame(anova_w_sandwich(object = m1_un_Dtg, objectglm = m1_un_dtg,coVars_woInt = c("DTG_pre")))

##Unadjusted Model
variable <- c("Female", "Greater then or equal to 50 ref(Less then 50)", "Tb (ref = No)",  
              "CD4 Greater then 350( ref = Less then 350)", "Greater then 6500 (ref =Less then 6500", "DTG (ref = Post DTG)")

rr <- c(round(exp(summ_un_female$Effect[1]), digits =3),round(exp(summ_un_age_gt50$Effect[1]), digits =3),round(exp(summ_un_tb$Effect[1]), digits =3),
        round(exp(summ_un_cd4$Effect[1]), digits =3),round(exp(summ_un_rna$Effect[1]), digits =3),round(exp(summ_un_dtg$Effect[1]), digits =3))                                                

ll <- c(round(exp(summ_un_female$`Lower 0.95`[1]), digits =3),round(exp(summ_un_age_gt50$`Lower 0.95`[1]), digits =3),round(exp(summ_un_tb$`Lower 0.95`[1]), digits =3),
        round(exp(summ_un_cd4$`Lower 0.95`[1]), digits =3),round(exp(summ_un_rna$`Lower 0.95`[1]), digits =3),round(exp(summ_un_dtg$`Lower 0.95`[1]), digits =3))                                                


ul <- c(round(exp(summ_un_female$`Upper 0.95`[1]), digits =3),round(exp(summ_un_age_gt50$`Upper 0.95`[1]), digits =3),round(exp(summ_un_tb$`Upper 0.95`[1]), digits =3),
        round(exp(summ_un_cd4$`Upper 0.95`[1]), digits =3),round(exp(summ_un_rna$`Upper 0.95`[1]), digits =3),round(exp(summ_un_dtg$`Upper 0.95`[1]), digits =3)) 

p <- c(anova_un_female$stats[1], anova_un_age_gt50$stats[1], anova_un_tb$stats[1],anova_un_cd4$stats[1], anova_un_rna$stats[1],anova_un_dtg$stats[1])

unadjusted_h <- data.frame(variable,rr,ll,ul,p)


## ----------------------  ##
##         Model 2         ##
## ----------------------  ##
# dummy coding
m1_dummy_interaction <- glm(DTG ~   age_gender_period, data = demdata_h, family=poisson(link="log"))
m1_dummy_Interaction <- Glm(DTG ~   age_gender_period, data = demdata_h, family=poisson(link="log"))
summ_dummy_int <- as.data.frame(my_summary.rms(object=m1_dummy_Interaction, object2 = m1_dummy_interaction))
object2 = m1_dummy_interaction
var_dummy_int <- sandwich(object2)[-1, -1] # Removes Intercept row/column
exposures <- " "
anova_dummy_int <- as.data.frame(anova_w_sandwich(object = m1_dummy_Interaction, objectglm = m1_dummy_interaction,coVars_woInt = c("age_gender_period")))


variables <- c("older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(round(exp(summ_dummy_int$Effect[1]), digits =3),round(exp(summ_dummy_int$Effect[2]), digits =3),round(exp(summ_dummy_int$Effect[3]), digits =3)
        ,round(exp(summ_dummy_int$Effect[4]), digits =3),round(exp(summ_dummy_int$Effect[5]), digits =3),round(exp(summ_dummy_int$Effect[6]), digits =3)
        ,round(exp(summ_dummy_int$Effect[7]), digits =3))
ll <- c(round(exp(summ_dummy_int$`Lower 0.95`[1]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[2]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[3]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[4]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[5]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[6]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[7]), digits =3))

ul <- c(round(exp(summ_dummy_int$`Upper 0.95`[1]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[2]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[3]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[4]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[5]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[6]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[7]), digits =3))

z1 <- summ_dummy_int$Effect[1]/summ_dummy_int$S.E.[1]
z2 <- summ_dummy_int$Effect[2]/summ_dummy_int$S.E.[2]
z3 <- summ_dummy_int$Effect[3]/summ_dummy_int$S.E.[3]
z4 <- summ_dummy_int$Effect[4]/summ_dummy_int$S.E.[4]
z5 <- summ_dummy_int$Effect[5]/summ_dummy_int$S.E.[5]
z6 <- summ_dummy_int$Effect[6]/summ_dummy_int$S.E.[6]
z7 <- summ_dummy_int$Effect[7]/summ_dummy_int$S.E.[7]


p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3))



dummy_interaction_h <-  data.frame(variables,rr,ll,ul,p)


## ----------------------  ##
##         Model 3         ##
## ----------------------  ##

## Interaction terms model
m1_interaction <- glm(DTG ~    age_gt50.factor*DTG_pre.factor*gender.factor , data= demdata_h, family=poisson(link="log"))
summ_int <- summary(m1_interaction)
var <- sandwich(m1_interaction)



#OMP
ind_test = 2
test_beta = summ_int$coefficients[ind_test,1]
OMP_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
OMP_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(var[ind_test,ind_test])



#OFP
ind_test = c(2,4,6)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(2,3,4,5,6,7,8)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(2,3,5)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 4
test_beta = summ_int$coefficients[ind_test]
YFP_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
YFP_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(var[ind_test,ind_test])

#YFPr
ind_test = c(3,4,7)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = summ_int$coefficients[ind_test]
YMPr_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(var[ind_test,ind_test])


## young female pre vs young female post
ind_test = c(3,7)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


interaction_h <-  data.frame(variable, rr,ll,ul,p)



## -------------------------------------------  ##
##           Adjusted RR for AIM 2 DTG          ##
## -------------------------------------------  ##

#Dummy interaction term
mod_dummy_a <- glm(DTG ~   tb_yes.factor + age_gender_period , data= demdata_h, family=poisson(link="log"))
mod_Dummy_a <- Glm(DTG ~   tb_yes.factor  + age_gender_period , data= demdata_h, family=poisson(link="log"))
summ_dummy_a <- as.data.frame(my_summary.rms(object = mod_Dummy_a, object2 = mod_dummy_a))
exposures <- " "
anova_cat_a <- as.data.frame(anova_w_sandwich(object = mod_Dummy_a, objectglm = mod_dummy_a,coVars_woInt = c(" tb_yes.factor ","age_gender_period")))

# Table for adjusted model with dummy variables
variable <- c("Tb (ref =no)","older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")

rr <- c(round(exp(summ_dummy_a$Effect[1]), digits =3),round(exp(summ_dummy_a$Effect[2]), digits =3),round(exp(summ_dummy_a$Effect[3]), digits =3)
        ,round(exp(summ_dummy_a$Effect[4]), digits =3),round(exp(summ_dummy_a$Effect[5]), digits =3),round(exp(summ_dummy_a$Effect[6]), digits =3)
        ,round(exp(summ_dummy_a$Effect[7]), digits =3),round(exp(summ_dummy_a$Effect[8]), digits =3))

ll <- c(round(exp(summ_dummy_a$`Lower 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[8]), digits =3))

ul <- c(round(exp(summ_dummy_a$`Upper 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[8]), digits =3))                                                                                                                                                                                   

z1 <- summ_dummy_a$Effect[1]/summ_dummy_a$S.E.[1]
z2 <- summ_dummy_a$Effect[2]/summ_dummy_a$S.E.[2]
z3 <- summ_dummy_a$Effect[3]/summ_dummy_a$S.E.[3]
z4 <- summ_dummy_a$Effect[4]/summ_dummy_a$S.E.[4]
z5 <- summ_dummy_a$Effect[5]/summ_dummy_a$S.E.[5]
z6 <- summ_dummy_a$Effect[6]/summ_dummy_a$S.E.[6]
z7 <- summ_dummy_a$Effect[7]/summ_dummy_a$S.E.[7]
z8 <- summ_dummy_a$Effect[8]/summ_dummy_a$S.E.[8]


p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits = 3))

interaction_dummy_a_h <- data.frame(variable,rr,ll,ul,p)


#Categorical interaction term
mod_cat <- glm(DTG ~   gender.factor + age_gt50.factor + tb_yes.factor + DTG_pre.factor + age_gt50.factor*DTG_pre.factor*gender.factor , data= demdata_h, family=poisson(link="log"))
summ_overall_a <- summary(mod_cat)
var_overall_a <- sandwich(mod_cat)


## Reporting overall Model with categorical interaction term

#OMP
ind = 3
beta = summ_overall_a$coefficients[ind,1]
OMP_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
OMP_point = exp(beta)
z1 <- beta/sqrt(var_overall_a[ind,ind])



#OFP
ind = c(2,3,7)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(2,3,5,6,7,8,9)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(3,5,6)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = summ_overall_a$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
YFP_point = exp(beta)
z5 <- beta/sqrt(var_overall_a[ind,ind])

#YFPr
ind = c(2,5,8)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 5
beta = summ_overall_a$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
YMPr_point = exp(beta)
z7 <- beta/sqrt(var_overall_a[ind,ind])


## young female pre vs young female post
ind = c(5,8)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(3,7)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,8)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,8,9)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)




#OFPr YFPr
ind = c(3,6,7,9)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
OFPr_YFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_YFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_YFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = summ_overall_a$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(var_overall_a[ind,ind])

#YFPr OFPr
ind = c(3,6)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
YMPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
YMPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
YMPr_OMPr_point = exp(beta)
z14 <- beta/sqrt(se)

ind = c(5,8,9,6)
var = as.numeric(var_overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_OFP_UB = exp(beta + 1.96*sqrt(se))
OFPr_OFP_LB = exp(beta - 1.96*sqrt(se))
OFPr_OFP_point = exp(beta)
z15 <- beta/sqrt(se)

ind = c(5,6)
var = as.numeric(var_overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OMPr_OMP_UB = exp(beta + 1.96*sqrt(se))
OMPr_OMP_LB = exp(beta - 1.96*sqrt(se))
OMPr_OMP_point = exp(beta)
z16 <- beta/sqrt(se)







var_overall <- c("TB (ref=no)",
                 "older male post DTG warning (ref = younger male post)", "older female postDTG warning (ref = younger male post)",
                 "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
                 "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
                 "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                 "older female pre  & during DTG warning: older male pre  & during DTG warning",
                 "younger female pre & during DTG warning : older female pre & during DTG warning","Younger Male Pre & during DTG warning:
                 Older Male pre & during warning", "Older female Pre: Older female Post", "Older Male Pre : Older Male Post")
rr <- c(tb_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,OFPr_YFPr_point,YMPr_OMPr_point,OFPr_OFP_point,OMPr_OMP_point)

ll <- c(tb_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, OFPr_YFPr_LB,YMPr_OMPr_LB,OFPr_OFP_LB,OMPr_OMP_LB)
ul <- c(tb_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, OFPr_YFPr_UB,YMPr_OMPr_UB,OFPr_OFP_UB,OMPr_OMP_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(z14))) * 2,digits=3),
       round((1 - pnorm(abs(z15))) * 2,digits=3),
       round((1 - pnorm(abs(z16))) * 2,digits=3))


haiti_a <-  data.frame(var_overall, rr,ll,ul,p)

## ANALYSIS for other sites except Haiti for AIM 1
demdata_1$study_site <- with(demdata_1, ifelse(study_site.factor == "Haiti",1,
                                   ifelse(study_site.factor == "Brazil",2,
                                          ifelse(study_site.factor == "Chile",3,4))))

demdata_o <- demdata_1 %>% filter(study_site != 1)
demdata_o$study_site.factor <- factor(demdata_o$study_site,
                                      levels=c(2,3,4),
                                      labels=c("Brazil", "Chile","Honduras" ))
                                    

## Setting reference "young male post"
demdata_o$DTG <- ifelse(demdata_o$DTG.factor == "No", 0, 1) 
dd <- datadist(demdata_o)
options(datadist="dd")
dd$limits$DTG.factor[2] <- "No"
# Setting reference levels for regular model
demdata_o <- within(demdata_o, DTG.factor <- relevel(factor(DTG.factor), ref = "No"))
demdata_o <- within(demdata_o, gender.factor <- relevel(factor(gender.factor), ref = "Male"))
demdata_o <- within(demdata_o, study_site.factor <- relevel(study_site.factor, ref = "Brazil"))
demdata_o <- within(demdata_o, age_gender_period <- relevel(factor(age_gender_period), ref = "younger male post"))
demdata_o <- within(demdata_o, age_gt50.factor <- relevel(factor(age_gt50.factor), ref = "Less then 50"))
demdata_o <- within(demdata_o, DTG_pre.factor <- relevel(factor(DTG_pre.factor), ref = "Post DTG"))
demdata_o <- within(demdata_o, female.factor <- relevel(factor(female.factor), ref = "Male"))
demdata_o <- within(demdata_o, tb_yes.factor <- relevel(factor(tb_yes.factor), ref = "No"))
demdata_o <- within(demdata_o, tb <- relevel(factor(tb), ref = "No"))
demdata_o <- within(demdata_o, ade_yes <- relevel(factor(ade_yes.factor), ref = "No"))
demdata_o <- within(demdata_o,cd4_gt350 <- relevel(factor(cd4_gt350.factor), ref = "Less then 350"))
demdata_o <- within(demdata_o,rna_gt6500 <- relevel(factor(rna_gt6500.factor), ref = "Less then 6500"))

#setting reference levels for Glm model
dd$limits$DTG.factor[2] <- "No"
dd$limits$gender.factor[2] <- "Male"
dd$limits$female.factor[2] <- "Male"
dd$limits$female[2] <- 0
dd$limits$age_gt50.factor[2] <- "Less then 50"
dd$limits$tb_yes[2] <- "No"
dd$limits$tb_yes.factor[2] <- "No"
dd$limits$ade_yes[2] <- "No"
dd$limits$cd4_gt350[2] <- "Less then 350"
dd$limits$rna_gt6500[2] <- "Less then 6500"
dd$limits$study_site.factor[2] <- "Brazil"
dd$limits$age_gender_period[2] <- "younger male post"
dd$limits$DTG_pre.factor[2] <- "Post DTG"
dd$limits$DTG[2] <- 0


## Unadjusted RR for individual variables 
#female
m1_un_female <- glm(DTG ~   female.factor, data = demdata_o, family=poisson(link="log"))
m1_un_Female <- Glm(DTG ~   female.factor, data = demdata_o, family=poisson(link="log"))
summ_un_female <- as.data.frame(my_summary.rms(object=m1_un_Female, object2 = m1_un_female))
exposures <- " "
anova_un_female <- as.data.frame(anova_w_sandwich(object = m1_un_Female, objectglm = m1_un_female,coVars_woInt = c("female")))

#age_gt50
m1_un_age_gt50 <- glm(DTG ~   age_gt50.factor, data = demdata_o, family=poisson(link="log"))
m1_un_Age_gt50 <- Glm(DTG ~   age_gt50.factor, data = demdata_o, family=poisson(link="log"))
summ_un_age_gt50 <- as.data.frame(my_summary.rms(object=m1_un_Age_gt50, object2 = m1_un_age_gt50))
exposures <- " "
anova_un_age_gt50 <- as.data.frame(anova_w_sandwich(object = m1_un_Age_gt50, objectglm = m1_un_age_gt50,coVars_woInt = c("age_gt50")))

#TB
m1_un_tb <- glm(DTG ~   tb_yes.factor, data = demdata_o, family=poisson(link="log"))
m1_un_Tb <- Glm(DTG ~   tb_yes.factor, data = demdata_o, family=poisson(link="log"))
summ_un_tb <- as.data.frame(my_summary.rms(object=m1_un_Tb, object2 = m1_un_tb))
exposures <- " "
anova_un_tb <- as.data.frame(anova_w_sandwich(object = m1_un_Tb, objectglm = m1_un_tb,coVars_woInt = c("tb_yes")))


#cd4
m1_un_cd4 <- glm(DTG ~   cd4_gt350.factor, data = demdata_o, family=poisson(link="log"))
m1_un_Cd4 <- Glm(DTG ~   cd4_gt350.factor, data = demdata_o, family=poisson(link="log"))
summ_un_cd4 <- as.data.frame(my_summary.rms(object=m1_un_Cd4, object2 = m1_un_cd4))
exposures <- " "
anova_un_cd4 <- as.data.frame(anova_w_sandwich(object = m1_un_Cd4, objectglm = m1_un_cd4,coVars_woInt = c("cd4_gt350")))


#rna
m1_un_rna <- glm(DTG ~   rna_gt6500.factor, data = demdata_o, family=poisson(link="log"))
m1_un_Rna <- Glm(DTG ~   rna_gt6500.factor, data = demdata_o, family=poisson(link="log"))
summ_un_rna <- as.data.frame(my_summary.rms(object=m1_un_Rna, object2 = m1_un_rna))
exposures <- " "
anova_un_rna <- as.data.frame(anova_w_sandwich(object = m1_un_Rna, objectglm = m1_un_rna,coVars_woInt = c("rna_gt6500")))


#site
m1_un_site <- glm(DTG ~   study_site.factor, data = demdata_o, family=poisson(link="log"))
m1_un_Site <- Glm(DTG ~   study_site.factor, data = demdata_o, family=poisson(link="log"))
summ_un_site <- as.data.frame(my_summary.rms(object=m1_un_Site, object2 = m1_un_site))
exposures <- " "
anova_un_site <- as.data.frame(anova_w_sandwich(object = m1_un_Site, objectglm = m1_un_site,coVars_woInt = c("study_site.factor")))


#DTG
m1_un_dtg <- glm(DTG ~   DTG_pre.factor, data = demdata_o, family=poisson(link="log"))
m1_un_Dtg <- Glm(DTG ~   DTG_pre.factor, data = demdata_o, family=poisson(link="log"))
summ_un_dtg <- as.data.frame(my_summary.rms(object=m1_un_Dtg, object2 = m1_un_dtg))
exposures <- " "
anova_un_dtg <- as.data.frame(anova_w_sandwich(object = m1_un_Dtg, objectglm = m1_un_dtg,coVars_woInt = c("DTG_pre")))

##Unadjusted Model
variable <- c("Female", "Greater then or equal to 50 ref(Less then 50)", "Tb (ref = No)",  "Site Chile (ref = Brazil)", "Site Honduras (ref = Brazil)",
              "CD4 Greater then 350( ref = Less then 350)", "Greater then 6500 (ref =Less then 6500", "DTG (ref = Post DTG)")

rr <- c(round(exp(summ_un_female$Effect[1]), digits =3),round(exp(summ_un_age_gt50$Effect[1]), digits =3),round(exp(summ_un_tb$Effect[1]), digits =3),
        round(exp(summ_un_site$Effect[1]), digits =3),round(exp(summ_un_site$Effect[2]), digits =3),
        round(exp(summ_un_cd4$Effect[1]), digits =3),round(exp(summ_un_rna$Effect[1]), digits =3),round(exp(summ_un_dtg$Effect[1]), digits =3))                                                

ll <- c(round(exp(summ_un_female$`Lower 0.95`[1]), digits =3),round(exp(summ_un_age_gt50$`Lower 0.95`[1]), digits =3),round(exp(summ_un_tb$`Lower 0.95`[1]), digits =3),
        round(exp(summ_un_site$`Lower 0.95`[1]), digits =3),round(exp(summ_un_site$`Lower 0.95`[2]), digits =3),
        round(exp(summ_un_cd4$`Lower 0.95`[1]), digits =3),round(exp(summ_un_rna$`Lower 0.95`[1]), digits =3),round(exp(summ_un_dtg$`Lower 0.95`[1]), digits =3))                                                


ul <- c(round(exp(summ_un_female$`Upper 0.95`[1]), digits =3),round(exp(summ_un_age_gt50$`Upper 0.95`[1]), digits =3),round(exp(summ_un_tb$`Upper 0.95`[1]), digits =3),
        round(exp(summ_un_site$`Upper 0.95`[1]), digits =3),round(exp(summ_un_site$`Upper 0.95`[2]), digits =3),
        round(exp(summ_un_cd4$`Upper 0.95`[1]), digits =3),round(exp(summ_un_rna$`Upper 0.95`[1]), digits =3),round(exp(summ_un_dtg$`Upper 0.95`[1]), digits =3)) 

p <- c(anova_un_female$stats[1], anova_un_age_gt50$stats[1], anova_un_tb$stats[1], anova_un_site$stats[1], "",anova_un_cd4$stats[1], anova_un_rna$stats[1],anova_un_dtg$stats[1])

unadjusted_o <- data.frame(variable,rr,ll,ul,p)

# dummy coding
m1_dummy_interaction <- glm(DTG ~   age_gender_period, data = demdata_o, family=poisson(link="log"))
m1_dummy_Interaction <- Glm(DTG ~   age_gender_period, data = demdata_o, family=poisson(link="log"))
summ_dummy_int <- as.data.frame(my_summary.rms(object=m1_dummy_Interaction, object2 = m1_dummy_interaction))
object2 = m1_dummy_interaction
var_dummy_int <- sandwich(object2)[-1, -1] # Removes Intercept row/column
exposures <- " "
anova_dummy_int <- as.data.frame(anova_w_sandwich(object = m1_dummy_Interaction, objectglm = m1_dummy_interaction,coVars_woInt = c("age_gender_period")))


variables <- c("older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(round(exp(summ_dummy_int$Effect[1]), digits =3),round(exp(summ_dummy_int$Effect[2]), digits =3),round(exp(summ_dummy_int$Effect[3]), digits =3)
        ,round(exp(summ_dummy_int$Effect[4]), digits =3),round(exp(summ_dummy_int$Effect[5]), digits =3),round(exp(summ_dummy_int$Effect[6]), digits =3)
        ,round(exp(summ_dummy_int$Effect[7]), digits =3))
ll <- c(round(exp(summ_dummy_int$`Lower 0.95`[1]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[2]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[3]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[4]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[5]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[6]), digits =3),
        round(exp(summ_dummy_int$`Lower 0.95`[7]), digits =3))

ul <- c(round(exp(summ_dummy_int$`Upper 0.95`[1]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[2]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[3]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[4]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[5]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[6]), digits =3),
        round(exp(summ_dummy_int$`Upper 0.95`[7]), digits =3))

z1 <- summ_dummy_int$Effect[1]/summ_dummy_int$S.E.[1]
z2 <- summ_dummy_int$Effect[2]/summ_dummy_int$S.E.[2]
z3 <- summ_dummy_int$Effect[3]/summ_dummy_int$S.E.[3]
z4 <- summ_dummy_int$Effect[4]/summ_dummy_int$S.E.[4]
z5 <- summ_dummy_int$Effect[5]/summ_dummy_int$S.E.[5]
z6 <- summ_dummy_int$Effect[6]/summ_dummy_int$S.E.[6]
z7 <- summ_dummy_int$Effect[7]/summ_dummy_int$S.E.[7]


p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3))



dummy_interaction_o <-  data.frame(variables,rr,ll,ul,p)


## ----------------------  ##
##         Model 3         ##
## ----------------------  ##

## Interaction terms model
m1_interaction <- glm(DTG ~    age_gt50.factor*DTG_pre.factor*gender.factor , data= demdata_o, family=poisson(link="log"))
summ_int <- summary(m1_interaction)
var <- sandwich(m1_interaction)



#OMP
ind_test = 2
test_beta = summ_int$coefficients[ind_test,1]
OMP_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
OMP_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(var[ind_test,ind_test])



#OFP
ind_test = c(2,4,6)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(2,3,4,5,6,7,8)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(2,3,5)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 4
test_beta = summ_int$coefficients[ind_test]
YFP_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
YFP_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(var[ind_test,ind_test])

#YFPr
ind_test = c(3,4,7)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = summ_int$coefficients[ind_test]
YMPr_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(var[ind_test,ind_test])


## young female pre vs young female post
ind_test = c(3,7)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


interaction_o <-  data.frame(variable, rr,ll,ul,p)


#Dummy interaction term
mod_dummy_a <- glm(DTG ~    study_site.factor + tb_yes.factor + age_gender_period , data= demdata_o, family=poisson(link="log"))
mod_Dummy_a <- Glm(DTG ~  study_site.factor  + tb_yes.factor  + age_gender_period , data= demdata_o, family=poisson(link="log"))
summ_dummy_a <- as.data.frame(my_summary.rms(object = mod_Dummy_a, object2 = mod_dummy_a))
exposures <- " "
anova_cat_a <- as.data.frame(anova_w_sandwich(object = mod_Dummy_a, objectglm = mod_dummy_a,coVars_woInt = c(" tb_yes.factor "," study_site.factor  ","age_gender_period")))

# Table for adjusted model with dummy variables
variable <- c("Chile","Honduras","Tb (ref =no)","older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")

rr <- c(round(exp(summ_dummy_a$Effect[1]), digits =3),round(exp(summ_dummy_a$Effect[2]), digits =3),round(exp(summ_dummy_a$Effect[3]), digits =3)
        ,round(exp(summ_dummy_a$Effect[4]), digits =3),round(exp(summ_dummy_a$Effect[5]), digits =3),round(exp(summ_dummy_a$Effect[6]), digits =3)
        ,round(exp(summ_dummy_a$Effect[7]), digits =3),round(exp(summ_dummy_a$Effect[8]), digits =3),round(exp(summ_dummy_a$Effect[9]), digits =3)
        ,round(exp(summ_dummy_a$Effect[10]), digits =3))

ll <- c(round(exp(summ_dummy_a$`Lower 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[8]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[9]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[10]), digits =3))

ul <- c(round(exp(summ_dummy_a$`Upper 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[8]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[9]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[10]), digits =3))                                                                                                                                                                                   

z1 <- summ_dummy_a$Effect[1]/summ_dummy_a$S.E.[1]
z2 <- summ_dummy_a$Effect[2]/summ_dummy_a$S.E.[2]
z3 <- summ_dummy_a$Effect[3]/summ_dummy_a$S.E.[3]
z4 <- summ_dummy_a$Effect[4]/summ_dummy_a$S.E.[4]
z5 <- summ_dummy_a$Effect[5]/summ_dummy_a$S.E.[5]
z6 <- summ_dummy_a$Effect[6]/summ_dummy_a$S.E.[6]
z7 <- summ_dummy_a$Effect[7]/summ_dummy_a$S.E.[7]
z8 <- summ_dummy_a$Effect[8]/summ_dummy_a$S.E.[8]
z9 <- summ_dummy_a$Effect[9]/summ_dummy_a$S.E.[9]
z10 <- summ_dummy_a$Effect[10]/summ_dummy_a$S.E.[10]


p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits = 3),round((1 - pnorm(abs(z9))) * 2,digits =3),round((1 - pnorm(abs(z10))) * 2,digits=3))

interaction_dummy_a_o <- data.frame(variable,rr,ll,ul,p)

#Categorical interaction term
mod_cat <- glm(DTG ~   gender.factor + age_gt50.factor + tb_yes.factor + study_site.factor   + DTG_pre.factor + age_gt50.factor*DTG_pre.factor*gender.factor , data= demdata_o, family=poisson(link="log"))
summ_overall_a <- summary(mod_cat)
var_overall_a <- sandwich(mod_cat)


## Reporting overall Model with categorical interaction term


#OMP
ind = 3
beta = summ_overall_a$coefficients[ind,1]
OMP_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
OMP_point = exp(beta)
z1 <- beta/sqrt(var_overall_a[ind,ind])



#OFP
ind = c(2,3,9)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(2,3,7,8,9,10,11)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(3,7,8)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = summ_overall_a$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
YFP_point = exp(beta)
z5 <- beta/sqrt(var_overall_a[ind,ind])

#YFPr
ind = c(2,7,10)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 7
beta = summ_overall_a$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
YMPr_point = exp(beta)
z7 <- beta/sqrt(var_overall_a[ind,ind])


## young female pre vs young female post
ind = c(7,10)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(3,9)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,9)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,10)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,9,10,11)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)




#OFPr YFPr
ind = c(3,8,9,11)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
OFPr_YFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_YFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_YFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = summ_overall_a$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(var_overall_a[ind,ind])

# ind = 5
# beta = summ_overall_a$coefficients[ind]
# ade_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
# ade_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
# ade_point = exp(beta)
# ade <- beta/sqrt(var_overall_a[ind,ind])



# ind = 5
# beta = summ_overall_a$coefficients[ind]
# bra_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
# bra_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
# bra_point = exp(beta)
# bra <- beta/sqrt(var_overall_a[ind,ind])

ind = 5
beta = summ_overall_a$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(var_overall_a[ind,ind])

ind = 6
beta = summ_overall_a$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(var_overall_a[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(var_overall_a[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(var_overall_a[ind,ind])

#YFPr OFPr
ind = c(3,8)
var = var_overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*summ_overall_a$coefficients[ind])
YMPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
YMPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
YMPr_OMPr_point = exp(beta)
z14 <- beta/sqrt(se)

ind = c(7,8,10,11)
var = as.numeric(var_overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OFPr_OFP_UB = exp(beta + 1.96*sqrt(se))
OFPr_OFP_LB = exp(beta - 1.96*sqrt(se))
OFPr_OFP_point = exp(beta)
z15 <- beta/sqrt(se)

ind = c(7,8)
var = as.numeric(var_overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(summ_overall_a$coefficients[ind])
OMPr_OMP_UB = exp(beta + 1.96*sqrt(se))
OMPr_OMP_LB = exp(beta - 1.96*sqrt(se))
OMPr_OMP_point = exp(beta)
z16 <- beta/sqrt(se)


var_overall <- c("Chile","Honduras",
                 "TB (ref=no)",
                 "older male post DTG warning (ref = younger male post)", "older female postDTG warning (ref = younger male post)",
                 "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
                 "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
                 "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                 "older female pre  & during DTG warning: older male pre  & during DTG warning",
                 "younger female pre & during DTG warning : older female pre & during DTG warning", "Younger Male Pre & during DTG warning:
                 Older Male pre & during warning", "Older female Pre: Older female Post", "Older Male Pre : Older Male Post")
rr <- c(chile_point,hon_point,tb_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,OFPr_YFPr_point,YMPr_OMPr_point,OFPr_OFP_point,OMPr_OMP_point)

ll <- c(chile_LB,hon_LB,tb_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, OFPr_YFPr_LB,YMPr_OMPr_LB,OFPr_OFP_LB,OMPr_OMP_LB)
ul <- c(chile_UB,hon_UB,tb_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, OFPr_YFPr_UB,YMPr_OMPr_UB,OFPr_OFP_UB,OMPr_OMP_UB)

p <- c(round((1 - pnorm(abs(chile))) * 2,digits = 3),
       round((1 - pnorm(abs(hon))) * 2,digits = 3),
       round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(z14))) * 2,digits=3),
       round((1 - pnorm(abs(z15))) * 2,digits=3),
       round((1 - pnorm(abs(z16))) * 2,digits=3))


interaction_a_o <-  data.frame(var_overall, rr,ll,ul,p)

save(unadjusted_h,file ="unadjusted_h.Rdata")
save(dummy_interaction_h,file ="dummy_interaction_h.Rdata")
save(interaction_h,file ="interaction_h.Rdata")
save(interaction_dummy_a_h,file ="interaction_dummy_a_h.Rdata")
save(haiti_a,file ="haiti_a.Rdata")
save(unadjusted_o,file ="unadjusted_o.Rdata")
save(dummy_interaction_o,file ="dummy_interaction_o.Rdata")
save(interaction_o,file ="interaction_o.Rdata")
save(interaction_dummy_a_o,file ="interaction_dummy_a_o.Rdata")
save(interaction_a_o,file ="interaction_a_o.Rdata")

# AIM 3 Haiti only
demdata3_h <- demdata_3 %>% filter(study_site.factor == "Haiti")

aim3_uni_h <- data.frame(Variable= c(
  c("Age (ref = Less then 50) ", "Greater then 50"),
  c("Age (ref = 35) ", 20,25,30,40,45,50,55,60),
  c("TB (ref = No)", "TB: Yes"),
  c("Gender (ref = Male)", "Females"),
  c("DTG (ref = no)", " DTG : Yes"),
  c("Timing (ref = Post DTG warning)", "Pre & during DTG warning")),
  RR= NA, Lower=NA, Upper=NA,  p=NA)



#age as categorical
m1 <- glm(undetectable4 ~ age_gt50.factor, demdata3_h, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
age50_beta = s1$coefficients[ind_atest]
age50_UB = exp(age50_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
age50_LB = exp(age50_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
age50_pointest = exp(age50_beta)
z1 <- age50_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_h$RR[aim3_uni_h$Variable %in% "Greater then 50"] <- age50_pointest
aim3_uni_h$Lower[aim3_uni_h$Variable %in% "Greater then 50"] <- age50_LB
aim3_uni_h$Upper[aim3_uni_h$Variable %in% "Greater then 50"] <- age50_UB
aim3_uni_h$p[aim3_uni_h$Variable %in% "Greater then 50"] <- p1

#age as spline
ages <- round(min(demdata3_h$age)):round(max(demdata3_h$age))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")



ms1 <- glm(undetectable4 ~ ns(age,df =4), demdata3_h, family= poisson(link="log"))
summ_spline  <- summary(ms1)
v_spline <- sandwich(ms1)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr <- function(myage,model){
  exp(model$Coefs[model$Variables %in% "ns(age, df = 4)1"]*
        (sp[rownames(sp) %in% myage, 1]-sp[rownames(sp) %in% 35, 1]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)2"]*
        (sp[rownames(sp) %in% myage, 2]- sp[rownames(sp) %in% 35, 2]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)3"]*
        (sp[rownames(sp) %in% myage, 3]- sp[rownames(sp) %in% 35, 3]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)4"]*
        (sp[rownames(sp) %in% myage, 4]- sp[rownames(sp) %in% 35, 4]))
  
}

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)


# my function to calculate variance and SE
spline_se_fun <- function(myage){
  
  tmp <- v_spline
  a <- sp$X1[rownames(sp) %in% myage]-sp$X1[rownames(sp) %in% 35]
  b <- sp$X2[rownames(sp) %in% myage]-sp$X2[rownames(sp) %in% 35]
  c <- sp$X3[rownames(sp) %in% myage]-sp$X3[rownames(sp) %in% 35]
  d <- sp$X4[rownames(sp) %in% myage]-sp$X4[rownames(sp) %in% 35]
  
  Var_x <- tmp["ns(age, df = 4)1", "ns(age, df = 4)1"]
  Var_y <- tmp["ns(age, df = 4)2", "ns(age, df = 4)2"]
  Var_z <- tmp["ns(age, df = 4)3", "ns(age, df = 4)3"]
  Var_w <- tmp["ns(age, df = 4)4", "ns(age, df = 4)4"]
  
  Cov_xy <- tmp["ns(age, df = 4)1", "ns(age, df = 4)2"]
  Cov_xz <- tmp["ns(age, df = 4)1", "ns(age, df = 4)3"]
  Cov_xw <- tmp["ns(age, df = 4)1", "ns(age, df = 4)4"]
  Cov_yz <- tmp["ns(age, df = 4)2", "ns(age, df = 4)3"]
  Cov_yw <- tmp["ns(age, df = 4)2", "ns(age, df = 4)4"]
  Cov_zw <- tmp["ns(age, df = 4)3", "ns(age, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI <- function(rr,SE) 
{ exp(log(rr) + c(-1,1)*1.96*SE)}


CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)


aim3_uni_h$RR[aim3_uni_h$Variable %in% 20] <- rr_age_20
aim3_uni_h$RR[aim3_uni_h$Variable %in% 25] <- rr_age_25
aim3_uni_h$RR[aim3_uni_h$Variable %in% 30] <- rr_age_30
aim3_uni_h$RR[aim3_uni_h$Variable %in% 40] <- rr_age_40
aim3_uni_h$RR[aim3_uni_h$Variable %in% 45] <- rr_age_45
aim3_uni_h$RR[aim3_uni_h$Variable %in% 50] <- rr_age_50
aim3_uni_h$RR[aim3_uni_h$Variable %in% 55] <- rr_age_55
aim3_uni_h$RR[aim3_uni_h$Variable %in% 60] <- rr_age_60

aim3_uni_h$Lower[aim3_uni_h$Variable %in% 20] <- CI_20[1]
aim3_uni_h$Lower[aim3_uni_h$Variable %in% 25] <- CI_25[1]
aim3_uni_h$Lower[aim3_uni_h$Variable %in% 30] <- CI_30[1]
aim3_uni_h$Lower[aim3_uni_h$Variable %in% 40] <- CI_40[1]
aim3_uni_h$Lower[aim3_uni_h$Variable %in% 45] <- CI_45[1]
aim3_uni_h$Lower[aim3_uni_h$Variable %in% 50] <- CI_50[1]
aim3_uni_h$Lower[aim3_uni_h$Variable %in% 55] <- CI_55[1]
aim3_uni_h$Lower[aim3_uni_h$Variable %in% 60] <- CI_60[1]

aim3_uni_h$Upper[aim3_uni_h$Variable %in% 20] <- CI_20[2]
aim3_uni_h$Upper[aim3_uni_h$Variable %in% 25] <- CI_25[2]
aim3_uni_h$Upper[aim3_uni_h$Variable %in% 30] <- CI_30[2]
aim3_uni_h$Upper[aim3_uni_h$Variable %in% 40] <- CI_40[2]
aim3_uni_h$Upper[aim3_uni_h$Variable %in% 45] <- CI_45[2]
aim3_uni_h$Upper[aim3_uni_h$Variable %in% 50] <- CI_50[2]
aim3_uni_h$Upper[aim3_uni_h$Variable %in% 55] <- CI_55[2]
aim3_uni_h$Upper[aim3_uni_h$Variable %in% 60] <- CI_60[2]

p <- function(rr,SE)
{
  round((1 - pnorm(abs(log(rr)/SE))) * 2,digits = 3)
}


aim3_uni_h$p[aim3_uni_h$Variable %in% 20] <- p(rr_age_20,SE_20)
aim3_uni_h$p[aim3_uni_h$Variable %in% 25] <- p(rr_age_25,SE_25)
aim3_uni_h$p[aim3_uni_h$Variable %in% 30] <- p(rr_age_30,SE_30)
aim3_uni_h$p[aim3_uni_h$Variable %in% 40] <- p(rr_age_40,SE_40)
aim3_uni_h$p[aim3_uni_h$Variable %in% 45] <- p(rr_age_45,SE_45)
aim3_uni_h$p[aim3_uni_h$Variable %in% 50] <- p(rr_age_50,SE_55)
aim3_uni_h$p[aim3_uni_h$Variable %in% 55] <- p(rr_age_55,SE_55)
aim3_uni_h$p[aim3_uni_h$Variable %in% 60] <- p(rr_age_60,SE_60)


#TB
demdata3_h$tb_yes <- ifelse(demdata3_h$tb == "Yes", 1,0)
demdata3_h$tb_yes <- factor(demdata3_h$tb_yes,
                            levels=c(1,0),
                            labels=c("Yes","No"))
m1 <- glm(undetectable4 ~ tb_yes, demdata3_h, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
tb_beta = s1$coefficients[ind_atest]
tb_UB = exp(tb_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
tb_LB = exp(tb_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
tb_pointest = exp(tb_beta)
z1 <- tb_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_h$RR[aim3_uni_h$Variable %in% "TB: Yes"] <- tb_pointest
aim3_uni_h$Lower[aim3_uni_h$Variable %in% "TB: Yes"] <- tb_LB
aim3_uni_h$Upper[aim3_uni_h$Variable %in% "TB: Yes"] <- tb_UB
aim3_uni_h$p[aim3_uni_h$Variable %in% "TB: Yes"] <- p1

m1 <- glm(undetectable4 ~ gender.factor, demdata3_h, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
gender_beta = s1$coefficients[ind_atest]
gender_UB = exp(gender_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
gender_LB = exp(gender_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
gender_pointest = exp(gender_beta)
z1 <- gender_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_h$RR[aim3_uni_h$Variable %in% "Females"] <- gender_pointest
aim3_uni_h$Lower[aim3_uni_h$Variable %in% "Females"] <- gender_LB
aim3_uni_h$Upper[aim3_uni_h$Variable %in% "Females"] <- gender_UB
aim3_uni_h$p[aim3_uni_h$Variable %in% "Females"] <- p1

demdata3_h <- within(demdata3_h, DTG_pre <- relevel(factor(DTG_pre), ref = "Post DTG"))
m1 <- glm(undetectable4 ~ DTG, demdata3_h, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
DTG_beta = s1$coefficients[ind_atest]
DTG_UB = exp(DTG_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
DTG_LB = exp(DTG_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
DTG_pointest = exp(DTG_beta)
z1 <- DTG_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_h$RR[aim3_uni_h$Variable %in% " DTG : Yes"] <- DTG_pointest
aim3_uni_h$Lower[aim3_uni_h$Variable %in% " DTG : Yes"] <- DTG_LB
aim3_uni_h$Upper[aim3_uni_h$Variable %in% " DTG : Yes"] <- DTG_UB
aim3_uni_h$p[aim3_uni_h$Variable %in% " DTG : Yes"] <- p1

## Timing
demdata3_h$DTG_pre <- ifelse(demdata3_h$DTG_warning == "Pre-warning" | demdata3_h$DTG_warning == "During warning", 1, 0)
demdata3_h$DTG_pre <- factor(demdata3_h$DTG_pre,
                             levels=c(1,0),
                             labels=c("Pre and during DTG", "Post DTG" ))

m1 <- glm(undetectable4 ~ DTG_pre, demdata3_h, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
timing_beta = s1$coefficients[ind_atest]
timing_UB = exp(timing_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
timing_LB = exp(timing_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
timing_pointest = exp(timing_beta)
z1 <- timing_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_h$RR[aim3_uni_h$Variable %in% "Pre & during DTG warning"] <- timing_pointest
aim3_uni_h$Lower[aim3_uni_h$Variable %in% "Pre & during DTG warning"] <- timing_LB
aim3_uni_h$Upper[aim3_uni_h$Variable %in% "Pre & during DTG warning"] <- timing_UB
aim3_uni_h$p[aim3_uni_h$Variable %in% "Pre & during DTG warning"] <- p1

demdata3_h$age_gender_period <- with(demdata3_h,ifelse(age_gt50 == 0 & DTG_pre == "Pre and during DTG" & gender == 0, "older male pre",
                                                       ifelse(age_gt50 == 1                & DTG_pre == "Pre and during DTG" & gender == 0, "younger male pre",
                                                              ifelse(age_gt50 == 0 & DTG_pre == "Post DTG"           & gender == 0, "older male post",
                                                                     ifelse(age_gt50 == 0 & DTG_pre == "Pre and during DTG" & gender == 1,  "older female pre",
                                                                            ifelse(age_gt50 == 1                & DTG_pre == "Post DTG"           & gender == 0,"younger male post",
                                                                                   ifelse(age_gt50 == 1               & DTG_pre == "Pre and during DTG" & gender == 1, "younger female pre",
                                                                                          ifelse(age_gt50 == 0 & DTG_pre ==  "Post DTG" & gender == 1, "older female post",
                                                                                                 ifelse(age_gt50 == 1                & DTG_pre == "Post DTG"           & gender == 1,"younger female post", "Missing")))))))))      



m1 <- glm(undetectable4 ~  age_gender_period , demdata3_h, family= poisson(link="log"))
s1 <- summary(m1)
v1 <- sandwich(m1)

ind_atest = 2
ofp_beta = s1$coefficients[ind_atest]
ofp_UB = exp(ofp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ofp_LB = exp(ofp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ofp_pointest = exp(ofp_beta)
z1 <- ofp_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 3
ofpr_beta = s1$coefficients[ind_atest]
ofpr_UB = exp(ofpr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ofpr_LB = exp(ofpr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ofpr_pointest = exp(ofpr_beta)
z1 <- ofpr_beta/sqrt(v1[ind_atest,ind_atest])
p2 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 4
omp_beta = s1$coefficients[ind_atest]
omp_UB = exp(omp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
omp_LB = exp(omp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
omp_pointest = exp(omp_beta)
z1 <- omp_beta/sqrt(v1[ind_atest,ind_atest])
p3 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 5
ompr_beta = s1$coefficients[ind_atest]
ompr_UB = exp(ompr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ompr_LB = exp(ompr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ompr_pointest = exp(ompr_beta)
z1 <- ompr_beta/sqrt(v1[ind_atest,ind_atest])
p4 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 6
yfp_beta = s1$coefficients[ind_atest]
yfp_UB = exp(yfp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
yfp_LB = exp(yfp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
yfp_pointest = exp(yfp_beta)
z1 <- yfp_beta/sqrt(v1[ind_atest,ind_atest])
p5 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)


ind_atest = 7
yfpr_beta = s1$coefficients[ind_atest]
yfpr_UB = exp(yfpr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
yfpr_LB = exp(yfpr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
yfpr_pointest = exp(yfpr_beta)
z1 <- yfpr_beta/sqrt(v1[ind_atest,ind_atest])
p6 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 8
ympr_beta = s1$coefficients[ind_atest]
ympr_UB = exp(ympr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ympr_LB = exp(ympr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ympr_pointest = exp(ympr_beta)
z1 <- ympr_beta/sqrt(v1[ind_atest,ind_atest])
p7 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

variables <- c("older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(round(ofp_pointest, digits =3),round(ofpr_pointest, digits =3),round(omp_pointest, digits =3)
        ,round(ompr_pointest, digits =3),round(yfp_pointest, digits =3),round(yfpr_pointest, digits =3)
        ,round(ympr_pointest, digits =3))
ll <- c(round(ofp_LB, digits =3),
        round(ofpr_LB, digits =3),
        round(omp_LB, digits =3),
        round(ompr_LB, digits =3),
        round(yfp_LB, digits =3),
        round(yfpr_LB, digits =3),
        round(ympr_LB, digits =3))

ul <- c(round(ofp_UB, digits =3),
        round(ofpr_UB, digits =3),
        round(omp_UB, digits =3),
        round(ompr_UB, digits =3),
        round(yfp_UB, digits =3),
        round(yfpr_UB, digits =3),
        round(ympr_UB, digits =3))




p <- c(p1,p2,p3,p4,p5,p6,p7)

aim3_int_h <-  data.frame(variables,rr,ll,ul,p)



i3 <- glm(undetectable4 ~    tb_yes + age_gender_period , demdata3_h, family= poisson(link="log"))
is3 <- summary(i3)
iv3 <- sandwich(i3)


ind = 2
beta = is3$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv3[ind,ind])

ind = 3
beta = is3$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv3[ind,ind])

ind = 4
beta = is3$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv3[ind,ind])

ind = 5
beta = is3$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv3[ind,ind])

ind = 6
beta = is3$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv3[ind,ind])

ind = 7
beta = is3$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv3[ind,ind])

ind = 8
beta = is3$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv3[ind,ind])

ind = 9
beta = is3$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv3[ind,ind])

variables <- c("TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(tb_point,OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point)
ll <- c(tb_LB,OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB)
ul <- c(tb_UB,OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB)
p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3))



aim3_dummy_h <-  data.frame(variables, rr,ll,ul,p)

# Overall model with DTG
i4 <- glm(undetectable4 ~  tb_yes + age_gender_period  + DTG, demdata3_h, family= poisson(link="log"))
is4 <- summary(i4)
iv4<- sandwich(i4)



ind = 2
beta = is4$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv4[ind,ind])

ind = 3
beta = is4$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv4[ind,ind])

ind = 4
beta = is4$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv4[ind,ind])

ind = 5
beta = is4$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv4[ind,ind])

ind = 6
beta = is4$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv4[ind,ind])

ind = 7
beta = is4$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv4[ind,ind])

ind = 8
beta = is4$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv4[ind,ind])

ind = 9
beta = is4$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv4[ind,ind])

ind = 10
beta = is4$coefficients[ind]
DTG_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
DTG_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
DTG_point = exp(beta)
z_DTG <- beta/sqrt(iv4[ind,ind])


variables <- c("TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)", "DTG : Yes/No")
rr <- c(tb_point,OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point, DTG_point)
ll <- c(tb_LB,OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB, DTG_LB)
ul <- c(tb_UB,OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB, DTG_UB)
p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3),round((1 - pnorm(abs(z_DTG))) * 2,digits=3))


aim3_dummy_DTG_h <-  data.frame(variables, rr,ll,ul,p)

# UNdetectable2
i3 <- glm(undetectable2 ~    tb_yes + age_gender_period , demdata3_h, family= poisson(link="log"))
is3 <- summary(i3)
iv3 <- sandwich(i3)



ind = 2
beta = is3$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv3[ind,ind])

ind = 3
beta = is3$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv3[ind,ind])

ind = 4
beta = is3$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv3[ind,ind])

ind = 5
beta = is3$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv3[ind,ind])

ind = 6
beta = is3$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv3[ind,ind])

ind = 7
beta = is3$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv3[ind,ind])

ind = 8
beta = is3$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv3[ind,ind])

ind = 9
beta = is3$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv3[ind,ind])

variables <- c("TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(tb_point,OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point)
ll <- c(tb_LB,OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB)
ul <- c(tb_UB,OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB)
p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3))



undetectable2_overall_h <-  data.frame(variables, rr,ll,ul,p)


# Overall model with DTG
i4 <- glm(undetectable2 ~  tb_yes + age_gender_period  + DTG, demdata3_h, family= poisson(link="log"))
is4 <- summary(i4)
iv4<- sandwich(i4)



ind = 2
beta = is4$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv4[ind,ind])

ind = 3
beta = is4$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv4[ind,ind])

ind = 4
beta = is4$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv4[ind,ind])

ind = 5
beta = is4$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv4[ind,ind])

ind = 6
beta = is4$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv4[ind,ind])

ind = 7
beta = is4$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv4[ind,ind])

ind = 8
beta = is4$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv4[ind,ind])

ind = 9
beta = is4$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv4[ind,ind])

ind = 10
beta = is4$coefficients[ind]
DTG_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
DTG_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
DTG_point = exp(beta)
z_DTG <- beta/sqrt(iv4[ind,ind])


variables <- c("TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)", "DTG : Yes/No")
rr <- c(tb_point,OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point, DTG_point)
ll <- c(tb_LB,OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB, DTG_LB)
ul <- c(tb_UB,OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB, DTG_UB)
p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3),round((1 - pnorm(abs(z_DTG))) * 2,digits=3))

undetectable2_DTG_h <-  data.frame(variables, rr,ll,ul,p)




## AIM 3 for all sites except Haiti
demdata_3$study_site <- with(demdata_3, ifelse(study_site.factor == "Haiti",1,
                                                 ifelse(study_site.factor == "Brazil",2,
                                                        ifelse(study_site.factor == "Chile",3,4))))

demdata3_o <- demdata_3 %>% filter(study_site != 1)
demdata3_o$study_site.factor <- factor(demdata3_o$study_site,
                                       levels=c(2,3,4),
                                       labels=c("Brazil", "Chile","Honduras" ))

aim3_uni_o <- data.frame(Variable= c(
  c("Age (ref = Less then 50) ", "Greater then 50"),
  c("Age (ref = 35) ", 20,25,30,40,45,50,55,60),
  c("Site (ref= Brazil)", "Chile", "Honduras"),
  c("TB (ref = No)", "TB: Yes"),
  c("Gender (ref = Male)", "Females"),
  c("DTG (ref = no)", " DTG : Yes"),
  c("Timing (ref = Post DTG warning)", "Pre & during DTG warning")),
  RR= NA, Lower=NA, Upper=NA,  p=NA)



#age as categorical
m1 <- glm(undetectable4 ~ age_gt50.factor, demdata3_o, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
age50_beta = s1$coefficients[ind_atest]
age50_UB = exp(age50_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
age50_LB = exp(age50_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
age50_pointest = exp(age50_beta)
z1 <- age50_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_o$RR[aim3_uni_o$Variable %in% "Greater then 50"] <- age50_pointest
aim3_uni_o$Lower[aim3_uni_o$Variable %in% "Greater then 50"] <- age50_LB
aim3_uni_o$Upper[aim3_uni_o$Variable %in% "Greater then 50"] <- age50_UB
aim3_uni_o$p[aim3_uni_o$Variable %in% "Greater then 50"] <- p1

#age as spline
ages <- round(min(demdata3_o$age)):round(max(demdata3_o$age))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")



ms1 <- glm(undetectable4 ~ ns(age,df =4), demdata3_o, family= poisson(link="log"))
summ_spline  <- summary(ms1)
v_spline <- sandwich(ms1)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr <- function(myage,model){
  exp(model$Coefs[model$Variables %in% "ns(age, df = 4)1"]*
        (sp[rownames(sp) %in% myage, 1]-sp[rownames(sp) %in% 35, 1]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)2"]*
        (sp[rownames(sp) %in% myage, 2]- sp[rownames(sp) %in% 35, 2]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)3"]*
        (sp[rownames(sp) %in% myage, 3]- sp[rownames(sp) %in% 35, 3]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)4"]*
        (sp[rownames(sp) %in% myage, 4]- sp[rownames(sp) %in% 35, 4]))
  
}

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)


# my function to calculate variance and SE
spline_se_fun <- function(myage){
  
  tmp <- v_spline
  a <- sp$X1[rownames(sp) %in% myage]-sp$X1[rownames(sp) %in% 35]
  b <- sp$X2[rownames(sp) %in% myage]-sp$X2[rownames(sp) %in% 35]
  c <- sp$X3[rownames(sp) %in% myage]-sp$X3[rownames(sp) %in% 35]
  d <- sp$X4[rownames(sp) %in% myage]-sp$X4[rownames(sp) %in% 35]
  
  Var_x <- tmp["ns(age, df = 4)1", "ns(age, df = 4)1"]
  Var_y <- tmp["ns(age, df = 4)2", "ns(age, df = 4)2"]
  Var_z <- tmp["ns(age, df = 4)3", "ns(age, df = 4)3"]
  Var_w <- tmp["ns(age, df = 4)4", "ns(age, df = 4)4"]
  
  Cov_xy <- tmp["ns(age, df = 4)1", "ns(age, df = 4)2"]
  Cov_xz <- tmp["ns(age, df = 4)1", "ns(age, df = 4)3"]
  Cov_xw <- tmp["ns(age, df = 4)1", "ns(age, df = 4)4"]
  Cov_yz <- tmp["ns(age, df = 4)2", "ns(age, df = 4)3"]
  Cov_yw <- tmp["ns(age, df = 4)2", "ns(age, df = 4)4"]
  Cov_zw <- tmp["ns(age, df = 4)3", "ns(age, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI <- function(rr,SE) 
{ exp(log(rr) + c(-1,1)*1.96*SE)}


CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)


aim3_uni_o$RR[aim3_uni_o$Variable %in% 20] <- rr_age_20
aim3_uni_o$RR[aim3_uni_o$Variable %in% 25] <- rr_age_25
aim3_uni_o$RR[aim3_uni_o$Variable %in% 30] <- rr_age_30
aim3_uni_o$RR[aim3_uni_o$Variable %in% 40] <- rr_age_40
aim3_uni_o$RR[aim3_uni_o$Variable %in% 45] <- rr_age_45
aim3_uni_o$RR[aim3_uni_o$Variable %in% 50] <- rr_age_50
aim3_uni_o$RR[aim3_uni_o$Variable %in% 55] <- rr_age_55
aim3_uni_o$RR[aim3_uni_o$Variable %in% 60] <- rr_age_60

aim3_uni_o$Lower[aim3_uni_o$Variable %in% 20] <- CI_20[1]
aim3_uni_o$Lower[aim3_uni_o$Variable %in% 25] <- CI_25[1]
aim3_uni_o$Lower[aim3_uni_o$Variable %in% 30] <- CI_30[1]
aim3_uni_o$Lower[aim3_uni_o$Variable %in% 40] <- CI_40[1]
aim3_uni_o$Lower[aim3_uni_o$Variable %in% 45] <- CI_45[1]
aim3_uni_o$Lower[aim3_uni_o$Variable %in% 50] <- CI_50[1]
aim3_uni_o$Lower[aim3_uni_o$Variable %in% 55] <- CI_55[1]
aim3_uni_o$Lower[aim3_uni_o$Variable %in% 60] <- CI_60[1]

aim3_uni_o$Upper[aim3_uni_o$Variable %in% 20] <- CI_20[2]
aim3_uni_o$Upper[aim3_uni_o$Variable %in% 25] <- CI_25[2]
aim3_uni_o$Upper[aim3_uni_o$Variable %in% 30] <- CI_30[2]
aim3_uni_o$Upper[aim3_uni_o$Variable %in% 40] <- CI_40[2]
aim3_uni_o$Upper[aim3_uni_o$Variable %in% 45] <- CI_45[2]
aim3_uni_o$Upper[aim3_uni_o$Variable %in% 50] <- CI_50[2]
aim3_uni_o$Upper[aim3_uni_o$Variable %in% 55] <- CI_55[2]
aim3_uni_o$Upper[aim3_uni_o$Variable %in% 60] <- CI_60[2]

p <- function(rr,SE)
{
  round((1 - pnorm(abs(log(rr)/SE))) * 2,digits = 3)
}


aim3_uni_o$p[aim3_uni_o$Variable %in% 20] <- p(rr_age_20,SE_20)
aim3_uni_o$p[aim3_uni_o$Variable %in% 25] <- p(rr_age_25,SE_25)
aim3_uni_o$p[aim3_uni_o$Variable %in% 30] <- p(rr_age_30,SE_30)
aim3_uni_o$p[aim3_uni_o$Variable %in% 40] <- p(rr_age_40,SE_40)
aim3_uni_o$p[aim3_uni_o$Variable %in% 45] <- p(rr_age_45,SE_45)
aim3_uni_o$p[aim3_uni_o$Variable %in% 50] <- p(rr_age_50,SE_55)
aim3_uni_o$p[aim3_uni_o$Variable %in% 55] <- p(rr_age_55,SE_55)
aim3_uni_o$p[aim3_uni_o$Variable %in% 60] <- p(rr_age_60,SE_60)

# Study Site 
m_site <- glm(undetectable4 ~   study_site.factor, data = demdata3_o, family=poisson(link="log"))
summ_site  <- summary(m_site)
vsite <- sandwich(m_site)

ind_atest = 2
chile_beta = m_site$coefficients[ind_atest]
chile_UB = exp(chile_beta  + 1.96*sqrt(vsite[ind_atest,ind_atest]))
chile_LB = exp(chile_beta  - 1.96*sqrt(vsite[ind_atest,ind_atest]))
chile_pointest = exp(chile_beta )
z4 <- chile_beta /sqrt(vsite[ind_atest,ind_atest])
p4 <- round((1 - pnorm(abs(z4))) * 2,digits = 3)

aim3_uni_o$RR[aim3_uni_o$Variable %in% "Chile"] <- chile_pointest
aim3_uni_o$Lower[aim3_uni_o$Variable %in% "Chile" ]<- chile_LB
aim3_uni_o$Upper[aim3_uni_o$Variable %in% "Chile"] <- chile_UB
aim3_uni_o$p[aim3_uni_o$Variable %in% "Chile"] <- p4

ind_atest = 3
honduras_beta = m_site$coefficients[ind_atest]
honduras_UB = exp(honduras_beta  + 1.96*sqrt(vsite[ind_atest,ind_atest]))
honduras_LB = exp(honduras_beta  - 1.96*sqrt(vsite[ind_atest,ind_atest]))
honduras_pointest = exp(honduras_beta )
z5 <- honduras_beta /sqrt(vsite[ind_atest,ind_atest])
p5 <- round((1 - pnorm(abs(z5))) * 2,digits = 3)

aim3_uni_o$RR[aim3_uni_o$Variable %in% "Honduras"] <- honduras_pointest
aim3_uni_o$Lower[aim3_uni_o$Variable %in% "Honduras" ]<- honduras_LB
aim3_uni_o$Upper[aim3_uni_o$Variable %in% "Honduras"] <- honduras_UB
aim3_uni_o$p[aim3_uni_o$Variable %in% "Honduras"] <- p5


#TB
demdata3_o$tb_yes <- ifelse(demdata3_o$tb == "Yes", 1,0)
demdata3_o$tb_yes <- factor(demdata3_o$tb_yes,
                           levels=c(1,0),
                           labels=c("Yes","No"))
m1 <- glm(undetectable4 ~ tb_yes, demdata3_o, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
tb_beta = s1$coefficients[ind_atest]
tb_UB = exp(tb_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
tb_LB = exp(tb_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
tb_pointest = exp(tb_beta)
z1 <- tb_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_o$RR[aim3_uni_o$Variable %in% "TB: Yes"] <- tb_pointest
aim3_uni_o$Lower[aim3_uni_o$Variable %in% "TB: Yes"] <- tb_LB
aim3_uni_o$Upper[aim3_uni_o$Variable %in% "TB: Yes"] <- tb_UB
aim3_uni_o$p[aim3_uni_o$Variable %in% "TB: Yes"] <- p1

m1 <- glm(undetectable4 ~ gender.factor, demdata3_o, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
gender_beta = s1$coefficients[ind_atest]
gender_UB = exp(gender_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
gender_LB = exp(gender_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
gender_pointest = exp(gender_beta)
z1 <- gender_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_o$RR[aim3_uni_o$Variable %in% "Females"] <- gender_pointest
aim3_uni_o$Lower[aim3_uni_o$Variable %in% "Females"] <- gender_LB
aim3_uni_o$Upper[aim3_uni_o$Variable %in% "Females"] <- gender_UB
aim3_uni_o$p[aim3_uni_o$Variable %in% "Females"] <- p1

demdata3_o <- within(demdata3_o, DTG_pre <- relevel(factor(DTG_pre), ref = "Post DTG"))
m1 <- glm(undetectable4 ~ DTG, demdata3_o, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
DTG_beta = s1$coefficients[ind_atest]
DTG_UB = exp(DTG_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
DTG_LB = exp(DTG_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
DTG_pointest = exp(DTG_beta)
z1 <- DTG_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_o$RR[aim3_uni_o$Variable %in% " DTG : Yes"] <- DTG_pointest
aim3_uni_o$Lower[aim3_uni_o$Variable %in% " DTG : Yes"] <- DTG_LB
aim3_uni_o$Upper[aim3_uni_o$Variable %in% " DTG : Yes"] <- DTG_UB
aim3_uni_o$p[aim3_uni_o$Variable %in% " DTG : Yes"] <- p1

## Timing
demdata3_o$DTG_pre <- ifelse(demdata3_o$DTG_warning == "Pre-warning" | demdata3_o$DTG_warning == "During warning", 1, 0)
demdata3_o$DTG_pre <- factor(demdata3_o$DTG_pre,
                            levels=c(1,0),
                            labels=c("Pre and during DTG", "Post DTG" ))

m1 <- glm(undetectable4 ~ DTG_pre, demdata3_o, family= poisson(link="log"))
s1  <- summary(m1)
v1 <- sandwich(m1)


ind_atest = 2
timing_beta = s1$coefficients[ind_atest]
timing_UB = exp(timing_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
timing_LB = exp(timing_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
timing_pointest = exp(timing_beta)
z1 <- timing_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

aim3_uni_o$RR[aim3_uni_o$Variable %in% "Pre & during DTG warning"] <- timing_pointest
aim3_uni_o$Lower[aim3_uni_o$Variable %in% "Pre & during DTG warning"] <- timing_LB
aim3_uni_o$Upper[aim3_uni_o$Variable %in% "Pre & during DTG warning"] <- timing_UB
aim3_uni_o$p[aim3_uni_o$Variable %in% "Pre & during DTG warning"] <- p1

# demdata3_o$age_gender_period <- with(demdata3_o,ifelse(age_gt50 == 0 & DTG_pre == "Pre and during DTG" & gender == 0, "older male pre",
#                                                      ifelse(age_gt50 == 1                & DTG_pre == "Pre and during DTG" & gender == 0, "younger male pre",
#                                                             ifelse(age_gt50 == 0 & DTG_pre == "Post DTG"           & gender == 0, "older male post",
#                                                                    ifelse(age_gt50 == 0 & DTG_pre == "Pre and during DTG" & gender == 1,  "older female pre",
#                                                                           ifelse(age_gt50 == 1                & DTG_pre == "Post DTG"           & gender == 0,"younger male post",
#                                                                                  ifelse(age_gt50 == 1               & DTG_pre == "Pre and during DTG" & gender == 1, "younger female pre",
#                                                                                         ifelse(age_gt50 == 0 & DTG_pre ==  "Post DTG" & gender == 1, "older female post",
#                                                                                                ifelse(age_gt50 == 1                & DTG_pre == "Post DTG"           & gender == 1,"younger female post", "Missing")))))))))      
# 


m1 <- glm(undetectable4 ~  age_gender_period , demdata3_o, family= poisson(link="log"))
s1 <- summary(m1)
v1 <- sandwich(m1)

ind_atest = 2
ofp_beta = s1$coefficients[ind_atest]
ofp_UB = exp(ofp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ofp_LB = exp(ofp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ofp_pointest = exp(ofp_beta)
z1 <- ofp_beta/sqrt(v1[ind_atest,ind_atest])
p1 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 3
ofpr_beta = s1$coefficients[ind_atest]
ofpr_UB = exp(ofpr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ofpr_LB = exp(ofpr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ofpr_pointest = exp(ofpr_beta)
z1 <- ofpr_beta/sqrt(v1[ind_atest,ind_atest])
p2 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 4
omp_beta = s1$coefficients[ind_atest]
omp_UB = exp(omp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
omp_LB = exp(omp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
omp_pointest = exp(omp_beta)
z1 <- omp_beta/sqrt(v1[ind_atest,ind_atest])
p3 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 5
ompr_beta = s1$coefficients[ind_atest]
ompr_UB = exp(ompr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ompr_LB = exp(ompr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ompr_pointest = exp(ompr_beta)
z1 <- ompr_beta/sqrt(v1[ind_atest,ind_atest])
p4 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 6
yfp_beta = s1$coefficients[ind_atest]
yfp_UB = exp(yfp_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
yfp_LB = exp(yfp_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
yfp_pointest = exp(yfp_beta)
z1 <- yfp_beta/sqrt(v1[ind_atest,ind_atest])
p5 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)


ind_atest = 7
yfpr_beta = s1$coefficients[ind_atest]
yfpr_UB = exp(yfpr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
yfpr_LB = exp(yfpr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
yfpr_pointest = exp(yfpr_beta)
z1 <- yfpr_beta/sqrt(v1[ind_atest,ind_atest])
p6 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

ind_atest = 8
ympr_beta = s1$coefficients[ind_atest]
ympr_UB = exp(ympr_beta + 1.96*sqrt(v1[ind_atest,ind_atest]))
ympr_LB = exp(ympr_beta - 1.96*sqrt(v1[ind_atest,ind_atest]))
ympr_pointest = exp(ympr_beta)
z1 <- ympr_beta/sqrt(v1[ind_atest,ind_atest])
p7 <- round((1 - pnorm(abs(z1))) * 2,digits = 3)

variables <- c("older female post (ref = younger male post)", "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(round(ofp_pointest, digits =3),round(ofpr_pointest, digits =3),round(omp_pointest, digits =3)
        ,round(ompr_pointest, digits =3),round(yfp_pointest, digits =3),round(yfpr_pointest, digits =3)
        ,round(ympr_pointest, digits =3))
ll <- c(round(ofp_LB, digits =3),
        round(ofpr_LB, digits =3),
        round(omp_LB, digits =3),
        round(ompr_LB, digits =3),
        round(yfp_LB, digits =3),
        round(yfpr_LB, digits =3),
        round(ympr_LB, digits =3))

ul <- c(round(ofp_UB, digits =3),
        round(ofpr_UB, digits =3),
        round(omp_UB, digits =3),
        round(ompr_UB, digits =3),
        round(yfp_UB, digits =3),
        round(yfpr_UB, digits =3),
        round(ympr_UB, digits =3))




p <- c(p1,p2,p3,p4,p5,p6,p7)

aim3_int_o <-  data.frame(variables,rr,ll,ul,p)



i3 <- glm(undetectable4 ~   study_site.factor + tb_yes + age_gender_period , demdata3_o, family= poisson(link="log"))
is3 <- summary(i3)
iv3 <- sandwich(i3)


ind = 2
beta = is3$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(iv3[ind,ind])

ind = 3
beta = is3$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(iv3[ind,ind])


ind = 4
beta = is3$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv3[ind,ind])

ind = 5
beta = is3$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv3[ind,ind])

ind = 6
beta = is3$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv3[ind,ind])

ind = 7
beta = is3$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv3[ind,ind])

ind = 8
beta = is3$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv3[ind,ind])

ind = 9
beta = is3$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv3[ind,ind])

ind = 10
beta = is3$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv3[ind,ind])

ind = 11
beta = is3$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv3[ind,ind])

variables <- c("Chile","Honduras","TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(chile_point,hon_point,tb_point,
        OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point)
ll <- c(chile_LB,hon_LB,tb_LB,
        OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB)
ul <- c(chile_UB,hon_UB,tb_UB,
        OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB)
p <- c(round((1 - pnorm(abs(chile))) * 2,digits = 3),round((1 - pnorm(abs(hon))) * 2,digits = 3),round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3))

aim3_dummy_o <-  data.frame(variables, rr,ll,ul,p)

# Overall model with DTG
i4 <- glm(undetectable4 ~   study_site.factor + tb_yes + age_gender_period  + DTG, demdata3_o, family= poisson(link="log"))
is4 <- summary(i4)
iv4<- sandwich(i4)


ind = 2
beta = is4$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(iv4[ind,ind])

ind = 3
beta = is4$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(iv4[ind,ind])

ind = 4
beta = is4$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv4[ind,ind])

ind = 5
beta = is4$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv4[ind,ind])

ind = 6
beta = is4$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv4[ind,ind])

ind = 7
beta = is4$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv4[ind,ind])

ind = 8
beta = is4$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv4[ind,ind])

ind = 9
beta = is4$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv4[ind,ind])

ind = 10
beta = is4$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv4[ind,ind])

ind = 11
beta = is4$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv4[ind,ind])

ind = 12
beta = is4$coefficients[ind]
DTG_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
DTG_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
DTG_point = exp(beta)
z_DTG <- beta/sqrt(iv4[ind,ind])


variables <- c("Chile","Honduras","TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)", "DTG : Yes/No")
rr <- c(chile_point,hon_point,tb_point,
        OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point, DTG_point)
ll <- c(chile_LB,hon_LB,tb_LB,
        OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB, DTG_LB)
ul <- c(chile_UB,hon_UB,tb_UB,
        OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB, DTG_UB)
p <- c(round((1 - pnorm(abs(chile))) * 2,digits = 3),round((1 - pnorm(abs(hon))) * 2,digits = 3),round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3),round((1 - pnorm(abs(z_DTG))) * 2,digits=3))


aim3_dummy_DTG_o <-  data.frame(variables, rr,ll,ul,p)

# UNdetectable2
i3 <- glm(undetectable2 ~   study_site.factor + tb_yes + age_gender_period , demdata3_o, family= poisson(link="log"))
is3 <- summary(i3)
iv3 <- sandwich(i3)


ind = 2
beta = is3$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(iv3[ind,ind])

ind = 3
beta = is3$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(iv3[ind,ind])

ind = 4
beta = is3$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv3[ind,ind])

ind = 5
beta = is3$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv3[ind,ind])

ind = 6
beta = is3$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv3[ind,ind])

ind = 7
beta = is3$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv3[ind,ind])

ind = 8
beta = is3$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv3[ind,ind])

ind = 9
beta = is3$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv3[ind,ind])

ind = 10
beta = is3$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv3[ind,ind])

ind = 11
beta = is3$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv3[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv3[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv3[ind,ind])

variables <- c("Chile","Honduras","TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(chile_point,hon_point,tb_point,
        OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point)
ll <- c(chile_LB,hon_LB,tb_LB,
        OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB)
ul <- c(chile_UB,hon_UB,tb_UB,
        OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB)
p <- c(round((1 - pnorm(abs(chile))) * 2,digits = 3),round((1 - pnorm(abs(hon))) * 2,digits = 3),round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3))



undetectable2_overall_o <-  data.frame(variables, rr,ll,ul,p)

# Overall model with DTG
i4 <- glm(undetectable2 ~   study_site.factor + tb_yes + age_gender_period  + DTG, demdata3_o, family= poisson(link="log"))
is4 <- summary(i4)
iv4<- sandwich(i4)

ind = 2
beta = is4$coefficients[ind]
chile_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
chile_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
chile_point = exp(beta)
chile <- beta/sqrt(iv4[ind,ind])

ind = 3
beta = is4$coefficients[ind]
hon_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
hon_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
hon_point = exp(beta)
hon <- beta/sqrt(iv4[ind,ind])

ind = 4
beta = is4$coefficients[ind]
tb_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
tb_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
tb_point = exp(beta)
tb <- beta/sqrt(iv4[ind,ind])

ind = 5
beta = is4$coefficients[ind]
OFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFP_point = exp(beta)
OFP <- beta/sqrt(iv4[ind,ind])

ind = 6
beta = is4$coefficients[ind]
OFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OFPr_point = exp(beta)
OFPr <- beta/sqrt(iv4[ind,ind])

ind = 7
beta = is4$coefficients[ind]
OMP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMP_point = exp(beta)
OMP <- beta/sqrt(iv4[ind,ind])

ind = 8
beta = is4$coefficients[ind]
OMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
OMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
OMPr_point = exp(beta)
OMPr <- beta/sqrt(iv4[ind,ind])

ind = 9
beta = is4$coefficients[ind]
YFP_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFP_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFP_point = exp(beta)
YFP <- beta/sqrt(iv4[ind,ind])

ind = 10
beta = is4$coefficients[ind]
YFPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YFPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YFPr_point = exp(beta)
YFPr <- beta/sqrt(iv4[ind,ind])

ind = 11
beta = is4$coefficients[ind]
YMPr_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
YMPr_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
YMPr_point = exp(beta)
YMPr <- beta/sqrt(iv4[ind,ind])

ind = 12
beta = is4$coefficients[ind]
DTG_UB = exp(beta + 1.96*sqrt(iv4[ind,ind]))
DTG_LB = exp(beta - 1.96*sqrt(iv4[ind,ind]))
DTG_point = exp(beta)
z_DTG <- beta/sqrt(iv4[ind,ind])


variables <- c("Chile","Honduras","TB (ref=no)","older female post (ref = younger male post)", 
               "older female pre (ref = younger male post)","older male post (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)", "DTG : Yes/No")
rr <- c(chile_point,hon_point,tb_point,
        OFP_point,OFPr_point,OMP_point,OMPr_point,YFP_point,YFPr_point,YMPr_point, DTG_point)
ll <- c(chile_LB,hon_LB,tb_LB,
        OFP_LB,OFPr_LB,OMP_LB,OMPr_LB,YFP_LB,YFPr_LB,YMPr_LB, DTG_LB)
ul <- c(chile_UB,hon_UB,tb_UB,
        OFP_UB,OFPr_UB,OMP_UB,OMPr_UB,YFP_UB,YFPr_UB,YMPr_UB, DTG_UB)
p <- c(round((1 - pnorm(abs(chile))) * 2,digits = 3),round((1 - pnorm(abs(hon))) * 2,digits = 3),round((1 - pnorm(abs(tb))) * 2,digits = 3),round((1 - pnorm(abs(OFP))) * 2,digits = 3),round((1 - pnorm(abs(OFPr))) * 2,digits =3),round((1 - pnorm(abs(OMP))) * 2,digits=3),
       round((1 - pnorm(abs(OMPr))) * 2,digits =3),round((1 - pnorm(abs(YFP))) * 2,digits=3),round((1 - pnorm(abs(YFPr))) * 2,digits =3),round((1 - pnorm(abs(YMPr))) * 2,digits=3),round((1 - pnorm(abs(z_DTG))) * 2,digits=3))


undetectable2_DTG_o <-  data.frame(variables, rr,ll,ul,p)

# Save files
save(aim3_uni_h,file ="aim3_uni_h.Rdata")
save(aim3_int_h,file ="aim3_int_h.Rdata")
save(aim3_dummy_h,file ="aim3_dummy_h.Rdata")
save(aim3_dummy_DTG_h,file ="aim3_dummy_DTG_h.Rdata")
save(undetectable2_overall_h,file ="undetectable2_overall_h.Rdata")
save(undetectable2_DTG_h,file ="undetectable2_DTG_h.Rdata")

save(aim3_uni_o,file ="aim3_uni_o.Rdata")
save(aim3_int_o,file ="aim3_int_o.Rdata")
save(aim3_dummy_o,file ="aim3_dummy_o.Rdata")
save(aim3_dummy_DTG_o,file ="aim3_dummy_DTG_o.Rdata")
save(undetectable2_overall_o,file ="undetectable2_overall_o.Rdata")
save(undetectable2_DTG_o,file ="undetectable2_DTG_o.Rdata")

## ------------------------------------ ##
##     Data Management for DTG  study   ##
## ------------------------------------ ##

## ----------------------  ##
##     Loading Dataset     ##
## ----------------------  ##
library(dplyr)
library(tidyverse)
library(Hmisc)
library(purrr)
library(lubridate)
library(sandwich)
library(cmprsk)
library(survminer)

### Creating datasets and variables for AIM2 ###

#Filtering for ART Experienced with visit date after art available and baseline(min of visit_d,lab results,art switch date)


demdata_2 <- demdata  %>% filter(naive == 0) %>% select(-art_sd,-art_id,-site) 
dem2 <- demdata_2 %>% select(patient_id,AVAILABLE_DT)
b <- left_join(dem2,visit,by = "patient_id") %>% filter(visit_d >= AVAILABLE_DT) %>% group_by(patient_id) %>% slice_min(visit_d)
c <- left_join(dem2,lab_cd4, by = "patient_id") %>% filter(cd4_d >= AVAILABLE_DT) %>% group_by(patient_id) %>% slice_min(cd4_d)
d <- left_join(dem2,lab_rna,by = "patient_id") %>% filter(rna_d >= AVAILABLE_DT) %>% group_by(patient_id) %>% slice_min(rna_d)
e <- left_join(dem2,art,by = "patient_id") %>% filter(art_sd >= AVAILABLE_DT) %>% group_by(patient_id) %>% slice_min(art_sd)
b <- list(b,c,d,e) %>% reduce(full_join, by = "patient_id")
b <- b %>% group_by(patient_id) %>% mutate(baseline = min(visit_d, cd4_d,rna_d,art_sd, na.rm = T)) %>% select(patient_id,baseline)
demdata_2 <- left_join(b , demdata_2)
dem2 <- demdata_2 %>% select(patient_id,AVAILABLE_DT)
visit2 <- left_join(dem2, visit)
visit2$aim2 <-  ifelse(visit2$visit_d >= visit2$AVAILABLE_DT, 1,0)
visit2 <- visit2  %>% filter(visit_d >= AVAILABLE_DT ) 
visit2 <- visit2 %>% group_by(patient_id) %>% mutate(gap.time0 = as.numeric(difftime(as.Date(visit_d),lag(as.Date(visit_d)), units="weeks")/52.25))
visit2 <- visit2 %>% group_by(patient_id) %>% mutate(gap.time = max(gap.time0,na.rm=T))
visit2 <- visit2 %>% group_by(patient_id)  %>% slice_min(visit_d) %>%  select(patient_id,gap.time)
demdata_2 <- left_join(demdata_2,visit2, by = "patient_id") 
dem2 <- demdata_2 %>% select(patient_id)

## computing DTG switch and removing those with DTG as their first art regiment and started DTG before date of DTG availability 
art2 <- art %>% select(patient_id,art_sd,art_id)
art2 <- left_join(dem2,art2, by = "patient_id") 
art2$DTG0 <- ifelse(grepl("DTG",art2$art_id), 1, 0)
art2 <- art2 %>% group_by(patient_id) %>% mutate(DTG = sum(DTG0)) %>% mutate(switch = ifelse(DTG >= 1,1,0))  %>%
  mutate(DTG_start = ifelse(DTG0 == 1,art_sd, NA)) 
art2 <- art2 %>% group_by(patient_id) %>% mutate(DTG_d =  DTG_start[which(!is.na(DTG_start))[1]]) 
art2 <- art2 %>% 
  group_by(patient_id) %>% 
  slice_min(art_sd) %>% 
  mutate(art1 = ifelse(grepl("DTG",art_id), 1, 0)) %>% 
  filter(art1 == 0)
demdata_2 <- left_join(art2,demdata_2, by = "patient_id")
demdata_2$check <- ifelse(is.na(demdata_2$DTG_d) | demdata_2$DTG_d > demdata_2$AVAILABLE_DT, 1,0)
demdata_2 <- demdata_2 %>% filter(check ==1)



## Variables for Aim 2  ##


# age at baseline, since baseline is defined at first clinical visit
demdata_2$ageb <- round(as.numeric(difftime(as.Date(demdata_2$baseline),as.Date(demdata_2$birth_d), units="weeks")/52.25),digits = 0)
demdata_2$dtg_time <- abs(round(as.numeric(difftime(as.Date(demdata_2$baseline),as.Date(demdata_2$DTG_d), units="weeks")/52.25),digits = 0))
demdata_2$duration_art <- as.numeric(round(difftime(as.Date(demdata_2$baseline),as.Date(demdata_2$art_sd), units="weeks")/52.25),digits = 0)


# Year at Baseline defined from visit after DTG available, year of clinical enrollment from enrol_d #
demdata_2$baseline_year <- as.numeric(word(demdata_2$baseline, 1, sep = fixed('-')))
demdata_2$enrol_year <- as.numeric(word(demdata_2$enrol_d, 1, sep = fixed('-')))


# DTG warning variable 
demdata_2 <- demdata_2 %>% mutate(DTG_warning = if_else(baseline <= as.Date("2018-05-31"), "Pre-warning",
                                                        if_else(as.Date("2018-06-01") <= baseline & baseline <= as.Date("2019-07-31"), "During warning", "After warning")))




## Last Observation from clinic visit ##
#demdata_2$clinic_year <- as.numeric(word(demdata_2$l_alive_d, 1, sep = fixed('-')))
dem2 <- demdata_2 %>% select(patient_id,baseline)
## Create CD4 values ## 
dem2_cd4 <- left_join(dem2, lab_cd4, by = "patient_id")
dem2_cd4$upper_limit <- as.Date(dem2_cd4$baseline) + 30
dem2_cd4$lower_limit <- as.Date(dem2_cd4$baseline) %m-% months(6)
dem2_cd4 <- dem2_cd4 %>% filter(cd4_d   <= upper_limit & cd4_d >= lower_limit)
dem2_cd4 <- dem2_cd4 %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(cd4_d),units ="days")))
dem2_cd4 <- dem2_cd4 %>% group_by(patient_id) %>% slice_min(abs(difference))
dem2_cd4 <- dem2_cd4 %>% select(patient_id,cd4_v)

# create HIV RNA values ##
dem2_rna <- left_join(dem2, lab_rna, by = "patient_id")
dem2_rna$upper_limit <- as.Date(dem2_rna$baseline) + 30
dem2_rna$lower_limit <- as.Date(dem2_rna$baseline) %m-% months(6)
dem2_rna <- dem2_rna %>% filter(rna_d >= lower_limit)  
dem2_rna <- dem2_rna %>% filter(rna_d <= upper_limit)
dem2_rna <- dem2_rna %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(rna_d),units ="days")))
dem2_rna_min <- dem2_rna %>% filter(difference >= 0) %>%  group_by(patient_id) %>%  slice_min(difference) %>% select(patient_id,rna_v,difference)
dem2_rna_max <- dem2_rna %>% filter(difference < 0) %>% group_by(patient_id) %>% slice_max(difference) %>% select(patient_id,rna_v,difference) 
dem2_rna <- full_join(dem2_rna_max,dem2_rna_min, by = c("patient_id", "difference", "rna_v"))
dem2_rna <- dem2_rna %>% group_by(patient_id) %>% slice_max(difference)
dem2_rna$log_rna <- log10(abs(dem2_rna$rna_v))
dem2_rna$viral_supression <- ifelse((dem2_rna$rna_v) > 1000, "Detectable", "Undetectable")
dem2_rna$viral_supression <- ifelse(is.na(dem2_rna$viral_supression), "Undetectable", dem2_rna$viral_supression)
dem2_rna <- dem2_rna %>% distinct(patient_id, .keep_all = TRUE) %>% select(patient_id,rna_v, log_rna, viral_supression)


## History of Other AIDS Defining Illness ##
dem2_ce <- left_join(dem2, ce, by = "patient_id")
dem2_ce$upper_limit <- as.Date(dem2_ce$baseline) + 30
dem2_ce$lower_limit <- as.Date(dem2_ce$baseline) %m-% months(6)
dem2_ce <- dem2_ce %>% filter(ce_d  <= upper_limit & ce_d>= lower_limit)
dem2_ce <- dem2_ce %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(ce_d),units ="days")))
dem2_ce <- dem2_ce %>% group_by(patient_id) %>% slice_min(abs(difference))
dem2_ce <- dem2_ce %>% select(patient_id,ce_id)



## History of Other AIDS Defining Illness ##
dem2_ce <- left_join(dem2, ce, by = "patient_id")
dem2_ce$upper_limit <- as.Date(dem2_ce$baseline) + 30
dem2_ce$lower_limit <- as.Date(dem2_ce$baseline) %m-% months(6)
dem2_ce <- dem2_ce %>% filter(ce_d  <= upper_limit & ce_d>= lower_limit)
dem2_ce <- dem2_ce %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(ce_d),units ="days")))
dem2_ce <- dem2_ce %>% group_by(patient_id) %>% slice_min(abs(difference))
dem2_ce <- dem2_ce %>% select(patient_id,ce_id)

## History of TB ##
dem2_tb <- left_join(dem2, ce_tb, by = "patient_id")
dem2_tb$upper_limit <- as.Date(dem2_tb$baseline) + 30
dem2_tb$lower_limit <- as.Date(dem2_tb$baseline) %m-% months(6)
dem2_tb <- dem2_tb %>% filter(tbdiagnosis_d  <= upper_limit & tbdiagnosis_d >= lower_limit)
dem2_tb <- dem2_tb %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(tbdiagnosis_d),units ="days")))
dem2_tb <- dem2_tb %>% group_by(patient_id) %>% slice_min(abs(difference))
dem2_tb <- dem2_tb %>% mutate(tb = "Yes")
dem2_tb <- dem2_tb %>% select(patient_id,tb)

demdata_2 <- list(demdata_2, dem2_rna, dem2_cd4, dem2_ce, dem2_tb) %>% 
  reduce(left_join, by = "patient_id")
demdata_2 <- demdata_2 %>% distinct(patient_id, .keep_all = TRUE)
demdata_2 <- demdata_2 %>% mutate(ade_type= if_else(grepl('^ade', ce_id) & ce_id != "ade_tuberculosis" , 'Yes', 'No'))
demdata_2 <- demdata_2 %>% mutate(tb = (if_else(is.na(tb), "No", "Yes")))
demdata_2 <- demdata_2 %>% mutate(regiment = ifelse(ii1 > 0 | ii2 > 0, "IINSTI-based",
                                                    ifelse(nnrti > 0 & pi==0, "NNRTI-based",
                                                           ifelse(pi > 0, "PI-based", "Other"))))


# count of regiments till baseline or Number of prior regiment
reg <- left_join(dem2,art, by = "patient_id") %>% filter(art_sd < baseline) 
reg <- reg %>% mutate(regi = ifelse(ii1 > 0 | ii2 > 0, "IINSTI-based",
                                    ifelse(nnrti > 0 & pi==0, "NNRTI-based",
                                           ifelse(pi > 0, "PI-based", "Other"))))

reg <- reg %>% group_by(patient_id) %>% summarise(num = n_distinct(regi))
demdata_2 <- left_join(demdata_2,reg, by = "patient_id")

# Cohort closing date for respective countries
demdata_2 <- left_join(demdata_2,follow, by = "patient_id")
demdata_2$l_alive_d_n <- as.numeric(as.Date(demdata_2$l_alive_d))
quantile(demdata_2$l_alive_d_n[demdata_2$study_site.factor == "Brazil"],c(0.9,0.95, 0.97))
demdata_2$cd <- NULL
demdata_2$cd[demdata_2$study_site.factor == "Brazil"] <- as.Date(18983)


quantile(demdata_2$l_alive_d_n[demdata_2$study_site.factor == "Honduras"],c(0.9,0.95, 0.97))
demdata_2$cd[demdata_2$study_site.factor == "Honduras"] <- as.Date(19216)

quantile(demdata_2$l_alive_d_n[demdata_2$study_site.factor == "Chile"],c(0.9,0.95, 0.97))
demdata_2$cd[demdata_2$study_site.factor == "Chile"] <- as.Date(18972)

quantile(demdata_2$l_alive_d_n[demdata_2$study_site.factor == "Haiti"],c(0.9,0.95, 0.97))
demdata_2$cd[demdata_2$study_site.factor == "Haiti"] <- as.Date(19072)

# Variables for cumulative probabilities and survival analysis
demdata_2$ltfu <- ifelse(as.Date(demdata_2$cd) - as.Date(demdata_2$l_alive_d) > 365, 1,0)
demdata_2$drop_d <- ifelse(demdata_2$ltfu == 1, demdata_2$l_alive_d, demdata_2$drop_d)

demdata_2$art_date <- ifelse(is.na(demdata_2$DTG_d),demdata_2$art_sd,demdata_2$DTG_d)
demdata_2$cd <- as.Date(demdata_2$cd)

demdata_2$event <- ifelse(demdata_2$death_y == 1 & demdata_2$switch == 0 , "Dead",
                          ifelse(demdata_2$death_y == 0  & demdata_2$switch == 1,"Switched to DTG",
                                 ifelse(demdata_2$death_y == 1 & demdata_2$switch == 1 & demdata_2$death_d <= demdata_2$DTG_d,"Dead",
                                        ifelse(demdata_2$death_y == 1 & demdata_2$switch == 1 & demdata_2$death_d > demdata_2$DTG_d,"Switched to DTG","Censored"))))

demdata_2$event_d <- ifelse(demdata_2$event == "Dead" ,demdata_2$death_d,
                            ifelse(demdata_2$event == "Switched to DTG", demdata_2$DTG_d,
                                   ifelse(demdata_2$event == "Censored"  ,demdata_2$l_alive_d, NA)))

demdata_2$event.num <- with(demdata_2,ifelse(event == "Switched to DTG", 1,0))
demdata_2$time <- ifelse(demdata_2$event.num ==1,(as.numeric(as.Date(demdata_2$DTG_d) - as.Date(demdata_2$baseline)))/365.25 +1,(as.numeric(as.Date(demdata_2$l_alive_d) - as.Date(demdata_2$baseline)))/365.25)

demdata_2$time_care <- (as.numeric(as.Date(demdata_2$l_alive_d) - as.Date(demdata_2$baseline)))/365.25


demdata_2$time_event <- ifelse(demdata_2$event == "Dead",(as.numeric(as.Date(demdata_2$death_d) - as.Date(demdata_2$baseline)))/365.25,
                               ifelse(demdata_2$event == "Switched to DTG",(as.numeric(as.Date(demdata_2$DTG_d) - as.Date(demdata_2$baseline) + 1))/365.25 ,
                                      ifelse(demdata_2$event == "Censored"  ,(as.numeric(as.Date(demdata_2$l_alive_d) - as.Date(demdata_2$baseline)))/365.25, NA)))



# CD4 category 
demdata_2$cd4.factor <- as.factor(ifelse(demdata_2$cd4_v > 350, ">350 cells/ mm3", 
                                         ifelse(demdata_2$cd4_v <= 350,"<350 cells/ mm3", "Missing")))

demdata_2$age_gt50 <- ifelse(demdata_2$age >= 50, 1,0)
demdata_2$age_gt50.factor <- factor(demdata_2$age_gt50,
                                    levels=c(1,0),
                                    labels=c("Greater then or equal to 50","Less then 50"))

demdata_2$switch.factor <- factor(demdata_2$switch,
                                  levels=c(1,0),
                                  labels=c("Yes","No"))

demdata_2$death_y.factor <- factor(demdata_2$death_y,
                                   levels=c(1,0),
                                   labels=c("Yes","No"))

demdata_2$drop_y.factor <- factor(demdata_2$drop_y,
                                  levels=c(1,0),
                                  labels=c("Yes","No"))
demdata_2$gap <- ifelse(demdata_2$gap.time > 1,1,0)
demdata_2$gap.factor <- factor(demdata_2$gap,
                               levels=c(1,0),
                               labels=c("Yes","No"))

demdata_2$gender.factor <- ifelse(demdata_2$gender.factor == "Male", "Male", "Female")

demdata_2$DTG_pre <- ifelse(demdata_2$DTG_warning == "Pre-warning" |demdata_2$DTG_warning == "During warning", 1,0)
demdata_2$DTG_pre.factor <- factor(demdata_2$DTG_pre,
                                   levels=c(1,0),
                                   labels=c("Pre and during DTG", "Post DTG" ))

demdata_2$age_gender_period <- with(demdata_2,ifelse(age_gt50 == 0 & DTG_pre.factor == "Pre and during DTG" & gender.factor == "Male", "Younger male pre",
                                                     ifelse(age_gt50 == 1  & DTG_pre.factor == "Pre and during DTG" & gender.factor == "Male", "Older male pre",
                                                            ifelse(age_gt50 == 0 & DTG_pre.factor == "Post DTG"           & gender.factor == "Male", "Younger male post",
                                                                   ifelse(age_gt50 == 0 & DTG_pre.factor == "Pre and during DTG" & gender.factor == "Female",  "Younger female pre",
                                                                          ifelse(age_gt50 == 1  & DTG_pre.factor == "Post DTG"           & gender.factor == "Male","Older male post",
                                                                                 ifelse(age_gt50 == 1  & DTG_pre.factor == "Pre and during DTG" & gender.factor == "Female", "Older female pre",
                                                                                        ifelse(age_gt50 == 0 & DTG_pre.factor ==  "Post DTG" & gender.factor == "Female", "Younger female post",
                                                                                               ifelse(age_gt50 == 1  & DTG_pre.factor == "Post DTG" & gender.factor == "Female","Older female post", "Missing")))))))))      




## variables Graph
demdata_2$year <- word(demdata_2$art_date, 1, sep = fixed('-'))
demdata_2$month <- word(demdata_2$art_date, 2, sep = fixed('-'))
demdata_2$month<- ifelse(demdata_2$month == "01" | demdata_2$month == "02" | demdata_2$month == "03"|
                           demdata_2$month == "04" | demdata_2$month == "05" | demdata_2$month == "06", "06", "12")

demdata_2 <- demdata_2 %>% mutate(year_semi = ifelse( month == "06",paste(year, month,30, sep = "-"),
                                                      ifelse( month == "12",paste(year, month,31, sep = "-"), NA)))

demdata_2$year_semi <- as.Date(demdata_2$year_semi,format = "%Y-%m-%d")
demdata_2$alive_year <- word(demdata_2$l_alive_d, 1, sep = fixed('-'))
demdata_2$alive_month <- word(demdata_2$l_alive_d, 2, sep = fixed('-'))
demdata_2$alive_month<- ifelse(demdata_2$alive_month == "01" | demdata_2$alive_month == "02" | demdata_2$alive_month == "03"|
                                 demdata_2$alive_month == "04" | demdata_2$alive_month == "05" | demdata_2$alive_month == "06", "06", "12")

demdata_2 <- demdata_2 %>% mutate(alive_semi = ifelse( alive_month == "06",paste(alive_year, alive_month,30, sep = "-"),
                                                       ifelse( alive_month == "12",paste(alive_year, alive_month,31, sep = "-"), NA)))

demdata_2$bs_year <- word(demdata_2$baseline, 1, sep = fixed('-'))
demdata_2$bs_month <- word(demdata_2$baseline, 2, sep = fixed('-'))
demdata_2$bs_month<- ifelse(demdata_2$bs_month == "01" | demdata_2$bs_month == "02" | demdata_2$bs_month == "03"|
                              demdata_2$bs_month == "04" | demdata_2$bs_month == "05" | demdata_2$bs_month == "06", "06", "12")

demdata_2 <- demdata_2 %>% mutate(base_semi = ifelse( bs_month == "06",paste(bs_year, bs_month,30, sep = "-"),
                                                      ifelse( bs_month == "12",paste(bs_year, bs_month,31, sep = "-"), NA)))

demdata_2$swt_year <- word(demdata_2$DTG_d, 1, sep = fixed('-'))
demdata_2$swt_month <- word(demdata_2$DTG_d, 2, sep = fixed('-'))
demdata_2$swt_month<- ifelse(demdata_2$swt_month == "01" | demdata_2$swt_month == "02" | demdata_2$swt_month == "03"|
                               demdata_2$swt_month == "04" | demdata_2$swt_month == "05" | demdata_2$swt_month == "06", "06", "12")

demdata_2 <- demdata_2 %>% mutate(dtg_semi = ifelse( swt_month == "06",paste(swt_year, swt_month,30, sep = "-"),
                                                     ifelse( swt_month == "12",paste(swt_year, swt_month,31, sep = "-"), NA)))




## Proportions graph variables
demdata_2$denom2017a <- ifelse(demdata_2$baseline <= as.Date("2017-06-30") & demdata_2$l_alive_d > as.Date("2016-12-31"),1,0)
demdata_2$num2017a <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2017-06-30") & demdata_2$denom2017a ==1,1,0))

demdata_2$denom2017b <- ifelse(demdata_2$baseline <= as.Date("2017-12-31") & demdata_2$l_alive_d > as.Date("2017-06-30"),1,0)
demdata_2$num2017b <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2017-12-31") & demdata_2$denom2017b ==1,1,0))

demdata_2$denom2018a <- ifelse(demdata_2$baseline <= as.Date("2018-06-30") & demdata_2$l_alive_d > as.Date("2017-12-31"),1,0)
demdata_2$num2018a <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2018-06-30") & demdata_2$denom2018a ==1,1,0))

demdata_2$denom2018b <- ifelse(demdata_2$baseline <= as.Date("2018-12-31") & demdata_2$l_alive_d > as.Date("2018-06-30"),1,0)
demdata_2$num2018b <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2018-12-31") & demdata_2$denom2018b ==1,1,0))
demdata_2$denom2019a <- ifelse(demdata_2$baseline <= as.Date("2019-06-30") & demdata_2$l_alive_d > as.Date("2018-12-31"),1,0)
demdata_2$num2019a <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2019-06-30") & demdata_2$denom2019a ==1,1,0))

demdata_2$denom2019b <- ifelse(demdata_2$baseline <= as.Date("2019-12-31") & demdata_2$l_alive_d > as.Date("2019-06-30"),1,0)
demdata_2$num2019b <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2019-12-31") & demdata_2$denom2019b ==1,1,0))
demdata_2$denom2020a <- ifelse(demdata_2$baseline <= as.Date("2020-06-30") & demdata_2$l_alive_d > as.Date("2019-12-31"),1,0)
demdata_2$num2020a <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2020-06-30") & demdata_2$denom2020a ==1,1,0))

demdata_2$denom2020b <- ifelse(demdata_2$baseline <= as.Date("2020-12-31") & demdata_2$l_alive_d > as.Date("2020-06-30"),1,0)
demdata_2$num2020b <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2020-12-31") & demdata_2$denom2020b ==1,1,0))

demdata_2$denom2021a <- ifelse(demdata_2$baseline <= as.Date("2021-06-30") & demdata_2$l_alive_d > as.Date("2020-12-31"),1,0)
demdata_2$num2021a <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2021-06-30") & demdata_2$denom2021a ==1,1,0))

demdata_2$denom2021b <- ifelse(demdata_2$baseline <= as.Date("2021-12-31") & demdata_2$l_alive_d > as.Date("2021-06-30"),1,0)
demdata_2$num2021b <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2021-12-31") & demdata_2$denom2021b ==1,1,0))
demdata_2$denom2022a <- ifelse(demdata_2$baseline <= as.Date("2022-06-30") & demdata_2$l_alive_d > as.Date("2021-12-31"),1,0)
demdata_2$num2022a <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2022-06-30") & demdata_2$denom2022a ==1,1,0))

demdata_2$denom2022b <- ifelse(demdata_2$baseline <= as.Date("2022-12-31") & demdata_2$l_alive_d > as.Date("2022-06-30"),1,0)
demdata_2$num2022b <- ifelse(is.na(demdata_2$DTG_d),0,
                             ifelse(demdata_2$DTG_d <= as.Date("2022-12-31") & demdata_2$denom2022b ==1,1,0))

year <-  c("June 2017","December 2017","June 2018","December 2018","June 2019","December 2019",
           "June 2020","December 2020","June 2021","December 2021","June 2022","December 2022")
num <- c(sum(demdata_2$num2017a),sum(demdata_2$num2017b),sum(demdata_2$num2018a),sum(demdata_2$num2018b),sum(demdata_2$num2019a),sum(demdata_2$num2019b),
         sum(demdata_2$num2020a),sum(demdata_2$num2020b),sum(demdata_2$num2021a),sum(demdata_2$num2021b),sum(demdata_2$num2022a),sum(demdata_2$num2022b))

denom <- c(sum(demdata_2$denom2017a),sum(demdata_2$denom2017b),sum(demdata_2$denom2018a),sum(demdata_2$denom2018b),sum(demdata_2$denom2019a),sum(demdata_2$denom2019b),
           sum(demdata_2$denom2020a),sum(demdata_2$denom2020b),sum(demdata_2$denom2021a),sum(demdata_2$denom2021b),sum(demdata_2$denom2022a),sum(demdata_2$denom2022b))


prop <- c(sum(demdata_2$num2017a)/sum(demdata_2$denom2017a),sum(demdata_2$num2017b)/sum(demdata_2$denom2017b),sum(demdata_2$num2018a)/sum(demdata_2$denom2018a),sum(demdata_2$num2018b)/sum(demdata_2$denom2018b),
          sum(demdata_2$num2019a)/sum(demdata_2$denom2019a),sum(demdata_2$num2019b)/sum(demdata_2$denom2019b),sum(demdata_2$num2020a)/sum(demdata_2$denom2020a),sum(demdata_2$num2020b)/sum(demdata_2$denom2020b),
          sum(demdata_2$num2021a)/sum(demdata_2$denom2021a),sum(demdata_2$num2021b)/sum(demdata_2$denom2021b),sum(demdata_2$num2022a)/sum(demdata_2$denom2022a),sum(demdata_2$num2022b)/sum(demdata_2$denom2022b))

proportion <- data.frame(year,num,denom,prop)


t <- sum(demdata_2$num2017b)/sum(demdata_2$denom2017b)
t <- sum(demdata_2$num2017b[demdata_2$study_site.factor == "Chile"])/sum(demdata_2$denom2017b[demdata_2$study_site.factor == "Chile"])
t <- sum(demdata_2$num2017b[demdata_2$study_site.factor == "Haiti"])/sum(demdata_2$denom2017b[demdata_2$study_site.factor == "Haiti"])

t

# greater then 1 means atleast 1 switch as num is count of unique regiments
demdata_2$num_cat <- ifelse(demdata_2$num > 1, 1,0)
demdata_2$num_cat.factor <- factor(demdata_2$num_cat,
                                   levels=c(0,1),
                                   labels=c("One ART regiment before baseline","More then one ART regiment before baseline"))


## ----------------------  ##
##          Models         ##
## ----------------------  ##

# Re level for cox regression model

demdata_2 <- demdata_2 %>% select(-DTG_start,-DTG0,-art1,-recart_y , -art_rs3,-ii2,-art_rs4,-check,-naive)
demdata_2$viral_supression1 <- ifelse(is.na(demdata_2$viral_supression), "Undetectable", demdata_2$viral_supression)
demdata_2 <- within(demdata_2, age_gender_period <- relevel(factor(age_gender_period), ref = "Younger male post"))
demdata_2 <- within(demdata_2, DTG_pre.factor <- relevel(factor(DTG_pre.factor), ref = "Post DTG"))
demdata_2 <- within(demdata_2, gender.factor <- relevel(factor(gender.factor), ref = "Male"))
demdata_2 <- within(demdata_2, age_gt50.factor <- relevel(factor(age_gt50.factor), ref = "Less then 50"))
demdata_2 <- within(demdata_2, viral_supression <- relevel(factor(viral_supression), ref = "Undetectable"))
demdata_2 <- within(demdata_2, viral_supression1 <- relevel(factor(viral_supression1), ref = "Undetectable"))
demdata_2 <- within(demdata_2, num_cat.factor <- relevel(factor(num_cat.factor), ref = "One ART regiment before baseline"))


dd <- datadist(demdata_2)
options(datadist= 'dd')

dd$limits$age_gt50.factor[2] <- "Less then 50"
dd$limits$gender.factor[2] <- "Male"
dd$limits$tb[2] <- "No"
dd$limits$ade_type[2] <- "No"
dd$limits$DTG_pre.factor[2] <- "Post DTG"
dd$limits$age_gender_period[2] <- "Younger male post"
dd$limits$num_cat.factor[2] <- "One ART regiment before baseline"
dd$limits$viral_supression1[2] <- "Undetectable"


## Model ##

mi  <- cph(Surv(time, event.num) ~   age_gender_period + strat(study_site.factor), data = demdata_2, surv = TRUE, x = TRUE, y =TRUE)
s <- as.data.frame(summary(mi))
a <- as.data.frame(anova(mi))

variable <- c("older female post (ref = younger male post)", "older female pre (ref = younger male post)", "older male post (ref = younger male post)", "older male pre (ref = younger male post)","younger female post (ref = younger male post)", "younger female pre (ref = younger male post)","younger male pre(ref = younger male post)")
hr <- c(s$Effect[2],s$Effect[4],s$Effect[6],s$Effect[8],s$Effect[10],s$Effect[12],s$Effect[14])
ll <- c(s$`Lower 0.95`[2],s$`Lower 0.95`[4],s$`Lower 0.95`[6],s$`Lower 0.95`[8],s$`Lower 0.95`[10],s$`Lower 0.95`[12],s$`Lower 0.95`[14])
ul <- c(s$`Upper 0.95`[2],s$`Upper 0.95`[4],s$`Upper 0.95`[6],s$`Upper 0.95`[8],s$`Upper 0.95`[10],s$`Upper 0.95`[12],s$`Upper 0.95`[14])
p <- c(round(a$P[1],digits = 3 ), "","","","","","")

cox_i <- data.frame(variable,hr,ll,ul,p)

## Interaction terms model


m  <- cph(Surv(time, event.num) ~ age_gt50.factor*gender.factor*DTG_pre.factor + strat(study_site.factor), data= demdata_2, surv = TRUE, x = TRUE, y =TRUE)

overall_a <- data.matrix(m$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL


#OMP
ind_test = 1
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))



#OFP
ind_test = c(1,2,4)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(1,2,3,4,5,6,7)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(1,3,5)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 2
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))

#YFPr
ind_test = c(2,3,6)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))


## young female pre vs young female post
ind_test = c(2,6)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


interaction <-  data.frame(variable, rr,ll,ul,p)

m2  <- cph(Surv(time, event.num) ~  tb  + viral_supression1 + age_gender_period + strat(study_site.factor)+ num_cat, data = demdata_2, surv = TRUE, x = TRUE, y =TRUE)
s <- as.data.frame(summary(m2))
a <- as.data.frame(anova(m2))
#Model summary table
variable <- c(" TB (ref = no)", "Viral supression (ref = undetectable)" , "Number of ART regiments before baseline (ref = One)", 
              "older female post (ref = younger male post)", "older female pre (ref = younger male post)", "older male post (ref = younger male post)", "older male pre (ref = younger male post)","younger female post (ref = younger male post)", "younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")


hr <- c(s$Effect[4],s$Effect[6],s$Effect[2],s$Effect[8],s$Effect[10],s$Effect[12],s$Effect[14],s$Effect[16],s$Effect[18],s$Effect[20])
ll <- c(s$`Lower 0.95`[4],s$`Lower 0.95`[6],s$`Lower 0.95`[2],s$`Lower 0.95`[8],s$`Lower 0.95`[10],s$`Lower 0.95`[12],s$`Lower 0.95`[14],s$`Lower 0.95`[16],s$`Lower 0.95`[18],s$`Lower 0.95`[20])
ul <- c(s$`Upper 0.95`[4],s$`Upper 0.95`[6],s$`Upper 0.95`[2],s$`Upper 0.95`[8],s$`Upper 0.95`[10],s$`Upper 0.95`[12],s$`Upper 0.95`[14],s$`Upper 0.95`[16],s$`Upper 0.95`[18],s$`Upper 0.95`[20])
p <- c(round(a$P[1],digits = 3 ), round(a$P[2], digits = 3), round(a$P[4],digits=3),round(a$P[4], digits =3 ), "","","","","","")

cox <- data.frame(variable,hr,ll,ul,p)

m  <- cph(Surv(time, event.num) ~ age_gt50.factor + gender.factor + DTG_pre.factor + tb +  viral_supression1 
          + strat(study_site.factor) + age_gt50.factor*gender.factor*DTG_pre.factor + num_cat ,
          data = demdata_2, surv = TRUE, x = TRUE, y =TRUE) 

overall_a <- data.matrix(m$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL



var_overall <- c("TB (ref=no)",
                 "Viral supression (ref= undetectable)",
                 "older male post DTG warning (ref = younger male post)", "older female postDTG warning (ref = younger male post)",
                 "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
                 "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
                 "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                 "older female pre  & during DTG warning: older male pre  & during DTG warning",
                 "younger female pre & during DTG warning : older female pre & during DTG warning",
                 "Number of ART regiments before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m1$`m$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_a[ind,ind]))



#OFP
ind = c(1,2,7)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,7,8,9,10)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,8)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m1$`m$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

#YFPr
ind = c(2,3,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m1$`m$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_a[ind,ind]))


## young female pre vs young female post
ind = c(3,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,7)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,9,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,7,8,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m1$`m$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_a[ind,ind]))


ind = 5
beta = m1$`m$coefficients`[ind]
vl_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
vl_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
vl_point = exp(beta)
vl <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 6
beta = m1$`m$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_a[ind,ind]))

rr <- c(tb_point,vl_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,vl_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,vl_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(vl))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))


cox_a <-  data.frame(var_overall, rr,ll,ul,p)


## Imputed model 

set.seed(7)
cph_multi_i <- aregImpute(~ age_gt50.factor + gender.factor + DTG_pre.factor + tb + ade_type + viral_supression 
                          + strat(study_site.factor) + age_gt50.factor*gender.factor*DTG_pre.factor + num_cat.factor + time*event.num,
                          data = demdata_2, n.impute=20)

cph_multi <- fit.mult.impute(Surv(time, event.num) ~ age_gt50.factor + gender.factor + DTG_pre.factor + tb + ade_type + viral_supression 
                             + strat(study_site.factor) + age_gt50.factor*gender.factor*DTG_pre.factor + num_cat.factor, 
                             data = demdata_2, fitter=cph, fitargs=list(x=TRUE, y=TRUE,surv=TRUE), xtrans=cph_multi_i, n.impute=20)

overall_a <- data.matrix(cph_multi$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m3 <- as.data.frame(cph_multi$coefficients)
rownames(m3) <- NULL



# var_overall <- c("TB (ref=no)","Viral supression (ref= undetectable)",
#                  "older male post DTG warning", "older female postDTG warning",
#                  "older female pre  & during DTG warning","older male pre  & during DTG warning","younger female post DTG warning","younger female pre  & during DTG warning","younger male pre  & during DTG warning",
#                  "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
#                  "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
#                  "older female pre  & during DTG warning: older male pre  & during DTG warning",
#                  "younger female pre & during DTG warning : older female pre & during DTG warning",
#                  "Number of ART regiments before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_a[ind,ind]))



#OFP
ind = c(1,2,7)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,7,8,9,10)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,8)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m3$`cph_multi$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

#YFPr
ind = c(2,3,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m3$`cph_multi$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_a[ind,ind]))


## young female pre vs young female post
ind = c(3,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,7)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m3$`cph_multi$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,9,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,7,8,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m3$`cph_multi$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m3$`cph_multi$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_a[ind,ind]))

# ind = 5
# beta = m3$`cph_multi$coefficients`[ind]
# ade_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
# ade_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
# ade_point = exp(beta)
# ade <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 5
beta = m3$`cph_multi$coefficients`[ind]
vl_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
vl_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
vl_point = exp(beta)
vl <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 6
beta = m3$`cph_multi$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_a[ind,ind]))

rr <- c(tb_point,vl_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,vl_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,vl_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(vl))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))


cox_imp <-  data.frame(var_overall, rr,ll,ul,p)
demdata_2 <- demdata_2 %>% ungroup(patient_id)
demdata_2$clinic_year <- as.numeric(word(demdata_2$l_alive_d, 1, sep = fixed('-')))

## ----------------------  ##
##          Labels         ##
## ----------------------  ##
label(demdata_2$age)   <- "Age at baseline"
label(demdata_2$study_site.factor)   <- "Study site"
label(demdata_2$baseline_year) <- "Year of baseline"
label(demdata_2$enrol_year)   <- "Year of clinic enrollment"
label(demdata_2$DTG_warning) <- "DTG warning calendar era at baseline"
label(demdata_2$clinic_year) <- "Year of last clinic visit"
label(demdata_2$cd4_v)  <- "CD4 at baseline "
label(demdata_2$viral_supression1) <-  "HIV-1 RNA at baseline categorized as detectable vs undetectable"
label(demdata_2$viral_supression) <-  "HIV-1 RNA at baseline categorized as detectable vs undetectable (missing not categorised)"
label(demdata_2$log_rna) <- "HIV-1 RNA at baseline (log10-transformed) "
label(demdata_2$ade_type) <- "History of other AIDS-Defining illness (not TB)"
#label(demdata_2$DTG.factor) <- "DTG Given Yes/No"
label(demdata_2$tb) <- "TB Yes/No"
#label(demdata_2$base) <- "ART Class"
label(demdata_2$gender.factor) <- "Gender  "
label(demdata_2$duration_art)  <- "Time from ART initiation(in years) to baseline"
label(demdata_2$cd4.factor) <- "CD4 at baseline "
label(demdata_2$switch.factor) <- "Switched to DTG "
label(demdata_2$dtg_time) <- "Time from baseline to DTG switch(in years)"
label(demdata_2$age_gt50.factor) <- "Age(category) at baseline"
label(demdata_2$gap.time) <- "Gaps in care(in years)"
label(demdata_2$drop_y.factor) <- "End of follow-up (administrative censoring)"
label(demdata_2$death_y.factor) <- "Death"
label(demdata_2$gap.factor) <- "Gaps in care"
label(demdata_2$time_event) <- "Time from baseline to event or censoring"
label(demdata_2$time_care) <- "Time from baseline to last clinical activity"
label(demdata_2$num) <- "Number of regiments prior to baseline"
label(demdata_2$num_cat.factor) <- "Number of regiments before baseline"
label(demdata_2$regiment) <- "Type of regiment at baseline"





### Sensitivity analysis

m  <- cph(Surv(time, event.num) ~ age_gt50.factor + gender.factor + DTG_pre.factor + tb   
          + strat(study_site.factor) + age_gt50.factor*gender.factor*DTG_pre.factor + num_cat ,
          data = demdata_2, surv = TRUE, x = TRUE, y =TRUE) 

overall_a <- data.matrix(m$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL



var_vl <- c("TB (ref=no)",
            "older male post DTG warning (ref = younger male post)", "older female postDTG warning (ref = younger male post)",
            "older female pre  & during DTG warning(ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
            "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
            "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
            "older female pre  & during DTG warning: older male pre  & during DTG warning",
            "younger female pre & during DTG warning : older female pre & during DTG warning",
            "Number of ART regiments before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m1$`m$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_a[ind,ind]))



#OFP
ind = c(1,2,6)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,6,7,8,9)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,7)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m1$`m$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

#YFPr
ind = c(2,3,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m1$`m$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_a[ind,ind]))


## young female pre vs young female post
ind = c(3,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,6)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,6)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,6,8,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,6,7,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m1$`m$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_a[ind,ind]))

# ind = 5
# beta = m1$`m$coefficients`[ind]
# ade_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
# ade_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
# ade_point = exp(beta)
# ade <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 5
beta = m1$`m$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_a[ind,ind]))

rr <- c(tb_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))


cox_vl <-  data.frame(var_vl, rr,ll,ul,p)


save(demdata_2, file="demdata_2.Rdata")
save(cox, file="cox.Rdata")
save(cox_i, file="cox_i.Rdata")
save(cox_a, file="cox_a.Rdata")
save(cox_imp, file="cox_imp.Rdata")
save(cox_vl, file="cox_vl.Rdata")
save(interaction, file="interaction.Rdata")


## Time updated analysis ##


dd0 <- demdata_2 %>% select(patient_id,baseline,DTG_d,study_site.factor,gender.factor,age,ade_type,age_gt50,age_gt50.factor,l_alive_d,birth_d,num_cat.factor,num_cat,switch.factor)
dd1 <- dd0 %>% mutate(warning_d = as.Date("2019-07-31"))
# To avoid warning messages in case baseline is same as rna_d(baseline defined using rna_d,visit_d,cd4_d etc.)
dd1$start_d <- as.Date(dd1$baseline) -1 
dd1$end_d <- ifelse(is.na(dd1$DTG_d),dd1$l_alive_d,dd1$DTG_d)
dd1$start_d <- as.Date(dd1$start_d)
dd1$end_d <- as.Date(dd1$end_d)


## Add Time updated TB analysis to the AIM 2
dd.tb <- left_join(dd1,ce_tb,by = "patient_id") %>% select(tbdiagnosis_d,patient_id,start_d,end_d,baseline,study_site.factor) %>% filter(!is.na(tbdiagnosis_d))
dd.tb$baseline_6 <- as.Date(dd.tb$baseline) %m-% months(6)
dd.tb$tb <- with(dd.tb, ifelse(tbdiagnosis_d >= baseline_6, "Yes","No"))
dd.tb <- dd.tb %>% filter(tbdiagnosis_d >= baseline_6) %>% filter(tbdiagnosis_d <= end_d) %>% select(patient_id,tbdiagnosis_d,baseline,end_d,tb)
#dd.tb$start_date <- with(dd.tb,ifelse(tbdiagnosis_d == baseline | tbdiagnosis_d == end_d , as.Date(baseline) - 1,as.Date(tbdiagnosis_d)))
dd.tb$start_date <- as.Date(dd.tb$tbdiagnosis_d)
dd.tb$end_date <-  as.Date(dd.tb$tbdiagnosis_d)
#dd.tb$start_date <- as.Date(dd.tb$start_date)
#dd.tb$end_date <- as.Date(dd.tb$end_date)
dd.tb <- dd.tb %>% select(-baseline)
#dd.tb$tb <- "Yes"

## Creating long format for "period" variable 
dd1 <- dd1 %>% mutate(period0 = ifelse(start_d < warning_d | (start_d < end_d & end_d < warning_d), "Pre and during DTG",
                                       "Post DTG"))
dd1 <- dd1 %>% mutate(period1 = ifelse(end_d < warning_d, "Pre and during DTG",
                                       ifelse(end_d >= warning_d, "Post DTG",NA)))
dd1 <- dd1 %>% mutate(period1 = ifelse(period0 == period1, NA,period1))

dd1 <- dd1 %>% mutate(period0 = ifelse(start_d < warning_d | (start_d < end_d & end_d < warning_d), "Pre and during DTG",
                                       "Post DTG"))
dd1 <- dd1 %>% mutate(period1 = ifelse(end_d < warning_d, "Pre and during DTG",
                                       ifelse(end_d >= warning_d, "Post DTG",NA)))
dd1 <- dd1 %>% mutate(period1 = ifelse(period0 == period1, NA,period1))

dd.l <- gather(dd1, condition, period, period0:period1, factor_key=TRUE) 
dd.l <- dd.l %>% filter(!is.na(period))  
dd.l <- dd.l %>% mutate(event.num.up = ifelse(is.na(DTG_d),0,1))
dd.l <- dd.l %>% group_by(patient_id) %>% add_count()
dd.l <- dd.l %>% mutate(event.num.up =ifelse(n ==2 & !is.na(DTG_d) & condition == "period0" & end_d > warning_d,0,event.num.up))
dd.l <- dd.l %>% mutate(start_d = if_else(n == 2 & condition == "period1", as.Date(warning_d),as.Date(start_d )))
dd.l <- dd.l %>% mutate(end_d = if_else(n == 2 & condition == "period0", as.Date(warning_d),as.Date(end_d )))
dd.l <- dd.l %>% mutate(age.new = difftime(as.Date(baseline),as.Date(birth_d),units = "days")/365.25)
dd.l <- dd.l %>% select(-age,-condition,-period,-event.num.up,-n,-age,-birth_d,-age_gt50,-age_gt50.factor)


#RNA Long format
rna <- left_join(dd.l,lab_rna, by = "patient_id") 
# rna$lower_limit <- as.Date(rna$baseline) %m-% months(6)
#rna <- rna %>% filter(cd4_d   <= upper_limit & cd4_d >= lower_limit)
rna$flag <- with(rna, ifelse(is.na(DTG_d) | DTG_d>= rna_d  ,1,0))
rna <- rna %>% filter(flag == 1)

all =list(dd.l,rna)
all   <- do.call(bind_rows, all)
all$rna_d <- as.Date(all$rna_d)

# subtracting one date from rna_d when the start date is same as rna_d and if rna_d == DTG_d or rna_d == l_alive_d
all$rna_d <- with(all,ifelse(!is.na(rna_d) & as.Date(rna_d) == as.Date(baseline) |!is.na(DTG_d) & !is.na(rna_d) & as.Date(rna_d) == as.Date(DTG_d) | !is.na(rna_d) & as.Date(rna_d) == as.Date(l_alive_d), as.Date(rna_d)-1,as.Date(rna_d))) 
all$rna_d <- as.Date(all$rna_d)
# all$rna_d <- with(all,ifelse(rna_d == baseline , rna_d-1,rna_d))
# all$rna_d <- as.Date(all$rna_d)
all <- all %>% select(patient_id,start_d,end_d,DTG_d,rna_v,rna_d,ade_type,gender.factor,num_cat.factor,study_site.factor,num_cat,age.new,baseline)
all$start_date <- if_else(is.na(all$rna_d),all$start_d,all$rna_d)
all$end_date <- if_else(is.na(all$rna_d),all$end_d,all$rna_d)
all <- all %>% group_by(patient_id) %>% distinct(start_date, .keep_all = TRUE)
check <- all %>% select(patient_id,DTG_d,start_date,end_date,rna_d,baseline,start_d,end_d)


# Event
#all$event <- with(all,ifelse(is.na(DTG_d),0,1))

#all$event <- ifelse(!is.na(DTG_d) & start_date <= DTG_d,0,all$event)

all <- full_join(all,dd.tb, by = c("patient_id","start_date","end_date"))
all  <-  all %>% group_by(patient_id) %>% fill(DTG_d, .direction = "updown") 
all  <-  all %>% group_by(patient_id) %>% fill(tb, .direction = "updown") 
all  <-  all %>% group_by(patient_id) %>% fill(tbdiagnosis_d, .direction = "updown") 
all  <-  all %>% group_by(patient_id) %>% fill(baseline, .direction = "updown")
all  <-  all %>% group_by(patient_id) %>% fill(study_site.factor, .direction = "updown")
check <- all %>% filter( patient_id == 70232 ) %>% select(patient_id,DTG_d,start_date,end_date,rna_d,baseline,tbdiagnosis_d)
all$check <- with(all,ifelse(is.na(DTG_d),0,1))
all$check <- with(all,ifelse(DTG_d == tbdiagnosis_d &  tbdiagnosis_d == end_date & end_date == start_date,1,0))

all <- all %>% filter(check == 0 | is.na(check) )           
check7 <- all %>% filter(patient_id == 70232 ) %>% select(patient_id,DTG_d,start_date,end_date,rna_d,baseline,tbdiagnosis_d)

# check <- check %>% select(patient_id,DTG_d,start_date,end_date,rna_d,baseline,tbdiagnosis_d,tb)
# time updated viral failure
#all <- all %>% select(-start_d,-end_d)
all$viral_failure <- ifelse((all$rna_v) > 1000, "Yes", "No")
all <- all %>% group_by(patient_id) %>% arrange(start_date,.by_group=TRUE)
all  <-  all %>% group_by(patient_id) %>% fill(rna_d, .direction = "updown")
all  <-  all %>% group_by(patient_id) %>% fill(viral_failure, .direction = "down") ## changed the direction of filling the NA as we need the values for RNA 
all  <-  all %>% group_by(patient_id) %>% mutate(diff1 = start_date - lag(start_date))
all <- all %>% group_by(patient_id) %>% mutate(viral_failure.n = ifelse(start_date == rna_d,viral_failure,
                                                                        ifelse(diff1 >= 0 & diff1 <= 365, viral_failure,NA)))
all$bs2 <- as.Date(all$baseline) -1
#check8 <- all %>% filter(patient_id == 70232 ) %>% select(patient_id,DTG_d,start_date,end_date,rna_d,baseline,tbdiagnosis_d,bs2)

all <- all %>% filter(start_date >= bs2 )

check8 <- all %>% filter(patient_id == 70232 ) %>% select(patient_id,DTG_d,start_date,end_date,rna_d,baseline,tbdiagnosis_d)

#  Time variables 
all <- all %>% group_by(patient_id) %>% mutate(end_date = sort(end_date))
all$start_date <- as.Date(all$start_date, '%Y-%m-%d')
all$end_date <- as.Date(all$end_date, '%Y-%m-%d')
all <- all %>% group_by(patient_id) %>% mutate(diffDate = difftime(start_date, lag(start_date), units = "days"))
all <- all %>% group_by(patient_id) %>% mutate(diff = difftime(end_date, start_date, units = "days"))
all$diffDate <- as.numeric(all$diffDate)
all$diff <- as.numeric(all$diff)
all$diffDate <- ifelse(is.na(all$diffDate),0,all$diffDate)
all <- all %>% group_by(patient_id) %>% mutate(csum = cumsum(diffDate))
all <- all %>% group_by(patient_id) %>% mutate(csum1 = cumsum(diff))





# all$diff <- all$end_date -all$start_date
# shift <- function(x, n){
#   c(x[-(seq(n))], rep(NA, n))
# }
# all <- all %>% group_by(patient_id) %>% mutate(last = last(diff), first.indi = as.integer(row_number() == 1L),last.indi = as.integer(row_number() == n()))
#all$time0 <- with(all, ifelse(first.indi == 1 , 0,csum))
all$time0 <- all$csum
all$time1 <- all$csum1
#  all$time1 <- shift(all$csum, 1)
# all$time1 <- with(all, ifelse(last.indi == 1 , time0+ last,time1))
# 

# time updated Period 
all$period <- with(all,ifelse(start_date >= "2019-07-31","Post DTG warning","Pre and during DTG warning")) 
## TB Time updated 
all$tb <- with(all,ifelse(is.na(tbdiagnosis_d) | start_date < tbdiagnosis_d, "No","Yes"))
#all$tb <- ifelse(is.na(all$tb),"No",all$tb)
all <- all %>% select(patient_id,time0,time1,period,viral_failure.n,gender.factor,tb,ade_type,num_cat.factor,study_site.factor,num_cat,start_date,rna_d,start_date,end_date,DTG_d,age.new,tbdiagnosis_d,rna_v,baseline)

#age New
all$age_gt50  <- with(all,ifelse(age.new>= 50, 1,0))
all$age_gt50.factor <- factor(all$age_gt50,
                              levels=c(1,0),
                              labels=c("Greater then or equal to 50","Less then 50"))

all$age_gender_period <- with(all,ifelse(age_gt50 == 0 & period == "Pre and during DTG warning" & gender.factor == "Male", "Younger male pre",
                                         ifelse(age_gt50 == 1  & period == "Pre and during DTG warning" & gender.factor == "Male", "Older male pre",
                                                ifelse(age_gt50 == 0 & period == "Post DTG warning"           & gender.factor == "Male", "Younger male post",
                                                       ifelse(age_gt50 == 0 & period == "Pre and during DTG warning" & gender.factor == "Female",  "Younger female pre",
                                                              ifelse(age_gt50 == 1  & period == "Post DTG warning"           & gender.factor == "Male","Older male post",
                                                                     ifelse(age_gt50 == 1  & period == "Pre and during DTG warning" & gender.factor == "Female", "Older female pre",
                                                                            ifelse(age_gt50 == 0 & period ==  "Post DTG warning" & gender.factor == "Female", "Younger female post",
                                                                                   ifelse(age_gt50 == 1  & period == "Post DTG warning" & gender.factor == "Female","Older female post", "Missing")))))))))      



# Updating event occurring when the viral lab date and DTG date are same.
all <- all %>% group_by(patient_id) %>% mutate(flag_last = ifelse(row_number() == n(),1,0))
all$event <- with(all,ifelse(flag_last ==1 & !is.na(DTG_d) ==1 ,1,0))
all$event <- with(all,ifelse(!is.na(DTG_d) & DTG_d == end_date,1,0))
#all$inc <- with(all,ifelse(!is.na(DTG_d) & DTG_d == end_date & end_date == start_date,1,0))
all <- all %>% group_by(patient_id) %>% mutate(csum3 = cumsum(event))
all <- all %>% filter(csum3 <=1)
check0 <- all %>% filter(start_date == end_date)
all <- all %>% filter(start_date < end_date)
# check <- all %>% filter(start_date == end_date)
# test <- all %>% select(patient_id,start_date,end_date,baseline,DTG_d,rna_d,time0,time1)
# 
# t0 <- all %>% filter(patient_id == 2040)
# t1 <- all %>% filter(patient_id == 5046)
# t3 <- all %>% filter(patient_id == 5038)
# 
# t2 <- all %>% filter(start_date == end_date)
dd <- datadist(all)


all <- within(all, age_gender_period <- relevel(factor(age_gender_period), ref = "Younger male post"))
all <- within(all, period <- relevel(factor(period), ref = "Post DTG warning"))
all <- within(all, tb <- relevel(factor(tb), ref = "No"))
all <- within(all, ade_type <- relevel(factor(ade_type), ref = "No"))
all <- within(all, gender.factor <- relevel(factor(gender.factor), ref = "Male"))
all <- within(all, age_gt50.factor <- relevel(factor(age_gt50.factor), ref = "Less then 50"))
all <- within(all, viral_failure.n <- relevel(factor(viral_failure.n), ref = "No"))
all <- within(all, num_cat.factor <- relevel(factor(num_cat.factor), ref = "One ART regiment before baseline"))
#test <- all %>% select(patient_id,DTG_d,start_date,end_date,rna_d,event,time0,time1,viral_failure.n,baseline)

## Model ##
m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor*gender.factor*period + strata(study_site.factor), data= all, cluster = patient_id)

overall_a <- data.matrix(m$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL


#OMP
ind_test = 1
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))



#OFP
ind_test = c(1,2,4)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(1,2,3,4,5,6,7)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(1,3,5)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 2
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))

#YFPr
ind_test = c(2,3,6)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))


## young female pre vs young female post
ind_test = c(2,6)
test_var = overall_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


t.up.interaction  <-  data.frame(variable, rr,ll,ul,p)


mi  <- coxph(Surv(time0, time1,event) ~ age_gender_period + strata(study_site.factor), data= all, cluster = patient_id)

overall_a <- data.matrix(mi$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m1 <- as.data.frame(mi$coefficients)
rownames(m1) <- NULL
#OMP

ind_test = 3
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))

#OFP
ind_test = 1
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))


#OFPr
ind_test = 2
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))



#OMPr
ind_test = 4
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))

#YFP
ind_test = 5
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))


#YFPr
ind_test = 6
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))

#YMPr
ind_test = 7
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overall_a[ind_test,ind_test]))

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3))


t.up_i  <-  data.frame(variable, rr,ll,ul,p)



m2  <- coxph(Surv(time0, time1,event) ~  tb +  num_cat.factor  + viral_failure.n  + age_gender_period + strata(study_site.factor), data = all, cluster = patient_id)

overall_a <- data.matrix(m2$var)
m1 <- as.data.frame(m2$coefficients)

s <- summary(m2)
rownames(m1) <- NULL
s1 <- as.data.frame(s$conf.int)


#TB
ind = 1
beta = as.numeric(m1$`m2$coefficients`[ind])
z1 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 2
beta = as.numeric(m1$`m2$coefficients`[ind])
z2 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 3
beta = as.numeric(m1$`m2$coefficients`[ind])
z3 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 4
beta = as.numeric(m1$`m2$coefficients`[ind])
z4 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 5
beta = as.numeric(m1$`m2$coefficients`[ind])
z5 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 6
beta = as.numeric(m1$`m2$coefficients`[ind])
z6 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 7
beta = as.numeric(m1$`m2$coefficients`[ind])
z7 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 8
beta = as.numeric(m1$`m2$coefficients`[ind])
z8 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 9
beta = as.numeric(m1$`m2$coefficients`[ind])
z9 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 10
beta = as.numeric(m1$`m2$coefficients`[ind])
z10 <- beta/sqrt(as.numeric(overall_a[ind,ind]))


#Model summary table
variable <- c(" TB (ref = no)",  "Number of ART regiment before baseline (ref = One)","Viral failure (ref = No)" , 
              "older female post (ref = younger male post)", "older female pre (ref = younger male post)", "older male post (ref = younger male post)", "older male pre (ref = younger male post)","younger female post (ref = younger male post)", "younger female pre","younger male pre ")


hr <- c(s1$`exp(coef)`[1],s1$`exp(coef)`[2],s1$`exp(coef)`[3],s1$`exp(coef)`[4],s1$`exp(coef)`[5],s1$`exp(coef)`[6],s1$`exp(coef)`[7],s1$`exp(coef)`[8],s1$`exp(coef)`[9],s1$`exp(coef)`[10])
ll <- c(s1$`lower .95`[1],s1$`lower .95`[2],s1$`lower .95`[3],s1$`lower .95`[4],s1$`lower .95`[5],s1$`lower .95`[6],s1$`lower .95`[7],s1$`lower .95`[8],s1$`lower .95`[9],s1$`lower .95`[10])
ul <- c(s1$`upper .95`[1],s1$`upper .95`[2],s1$`upper .95`[3],s1$`upper .95`[4],s1$`upper .95`[5],s1$`upper .95`[6],s1$`upper .95`[7],s1$`upper .95`[8],s1$`upper .95`[9],s1$`upper .95`[10])
p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3))

t.up <- data.frame(variable,hr,ll,ul,p)

m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb  + viral_failure.n 
            + age_gt50.factor*gender.factor*period + num_cat.factor + strata(study_site.factor) ,
            data = all, cluster=patient_id) 

overall_a <- data.matrix(m$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL



var_overall <- c("TB (ref=no)","Viral failure (ref= No)",
                 "older male post DTG warning", "older female postDTG warning",
                 "older female pre  & during DTG warning","older male pre  & during DTG warning","younger female post DTG warning","younger female pre  & during DTG warning","younger male pre  & during DTG warning",
                 "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
                 "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                 "older female pre  & during DTG warning: older male pre  & during DTG warning",
                 "younger female pre & during DTG warning : older female pre & during DTG warning",
                 "Number of ART regiment before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m1$`m$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_a[ind,ind]))



#OFP
ind = c(1,2,7)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,7,8,9,10)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,8)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m1$`m$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

#YFPr
ind = c(2,3,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m1$`m$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_a[ind,ind]))


## young female pre vs young female post
ind = c(3,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,7)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m1$`m$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,9,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,7,8,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m1$`m$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_a[ind,ind]))

# ind = 5
# beta = m1$`m$coefficients`[ind]
# ade_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
# ade_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
# ade_point = exp(beta)
# ade <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 5
beta = m1$`m$coefficients`[ind]
vl_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
vl_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
vl_point = exp(beta)
vl <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 6
beta = m1$`m$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_a[ind,ind]))

rr <- c(tb_point,vl_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,vl_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,vl_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(vl))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))


t.up_a <-  data.frame(var_overall, rr,ll,ul,p)


## Imputed model 

set.seed(7)
cph_multi_i <- aregImpute(~ age_gt50.factor + gender.factor + period + tb  + viral_failure.n 
                          + strata(study_site.factor) + age_gt50.factor*gender.factor*period + num_cat.factor + time1*event,
                          data = all, n.impute=20)

cph_multi <- fit.mult.impute(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb + viral_failure.n 
                             + strat(study_site.factor) + age_gt50.factor*gender.factor*period + num_cat.factor, data = all,cluster=all$patient_id,
                             fitter=cph, fitargs=list(x=TRUE, y=TRUE,surv=TRUE), xtrans=cph_multi_i, n.impute=20)

overall_a <- data.matrix(cph_multi$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m3 <- as.data.frame(cph_multi$coefficients)
rownames(m3) <- NULL



var_overall <- c("TB (ref=no)",
                 "Viral failure (ref= No)",
                 "Younger male post DTG warning (ref = Older male post)", "older female postDTG warning (ref = younger male post)",
                 "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
                 "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
                 "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                 "older female pre  & during DTG warning: older male pre  & during DTG warning",
                 "younger female pre & during DTG warning : older female pre & during DTG warning",
                 "Number of ART regiment before baseline (ref = One)" , "Younger male pre : Older male pre","Older female Pre: Older female Post","Older Male Pre: Older Female Post")

#YMP_OMP
ind = 1
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
YMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMP_point = exp(beta*(-1))
z1 <- beta/sqrt(as.numeric(overall_a[ind,ind]))



#OFP
ind = c(1,2,7)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,7,8,9,10)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,8)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m3$`cph_multi$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

#YFPr
ind = c(2,3,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m3$`cph_multi$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_a[ind,ind]))


## young female pre vs young female post
ind = c(3,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#YFP vs OFP
ind = c(1,7)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m3$`cph_multi$coefficients`[ind])
YFP_OFP_UB = exp(beta + 1.96*sqrt(se))
YFP_OFP_LB = exp(beta - 1.96*sqrt(se))
YFP_OFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m3$`cph_multi$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,9,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference is opposite)
ind = c(1,7,8,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m3$`cph_multi$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m3$`cph_multi$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 5
beta = m3$`cph_multi$coefficients`[ind]
vl_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
vl_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
vl_point = exp(beta)
vl <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 6
beta = m3$`cph_multi$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = c(1,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m3$`cph_multi$coefficients`[ind])
YMPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
YMPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
YMPr_OMPr_point = exp(beta)
z14 <- beta/sqrt(se)

ind = c(3,8,9,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_OFP_UB = exp(beta + 1.96*sqrt(se))
OFPr_OFP_LB = exp(beta - 1.96*sqrt(se))
OFPr_OFP_point = exp(beta)
z15 <- beta/sqrt(se)

ind = c(3,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OMPr_OMP_UB = exp(beta + 1.96*sqrt(se))
OMPr_OMP_LB = exp(beta - 1.96*sqrt(se))
OMPr_OMP_point = exp(beta)
z16 <- beta/sqrt(se)

rr <- c(tb_point,vl_point,YMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,YFP_OFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point,YMPr_OMPr_point,OFPr_OFP_point,OMPr_OMP_point)

ll <- c(tb_LB,vl_LB,YMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,YFP_OFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB,YMPr_OMPr_LB,OFPr_OFP_LB,OMPr_OMP_LB)
ul <- c(tb_UB,vl_UB,YMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,YFP_OFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB,YMPr_OMPr_UB,OFPr_OFP_UB,OMPr_OMP_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(vl))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3),
       round((1 - pnorm(abs(z14))) * 2,digits=3),
       round((1 - pnorm(abs(z15))) * 2,digits=3),
       round((1 - pnorm(abs(z16))) * 2,digits=3))


t.up_imp <-  data.frame(var_overall, rr,ll,ul,p)
#all <- all %>% ungroup(patient_id)



### Sensitivity analysis

m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb   
            + strata(study_site.factor) + age_gt50.factor*gender.factor*period + num_cat.factor , 
            data = all, cluster= patient_id) 

overall_a <- data.matrix(m$var)
rownames(overall_a) <- NULL
colnames(overall_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL

var_vl <- c("TB (ref=no)",
            "older male post DTG warning (ref = younger male post)", "older female postDTG warning (ref = younger male post)",
            "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
            "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
            "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
            "older female pre  & during DTG warning: older male pre  & during DTG warning",
            "younger female pre & during DTG warning : older female pre & during DTG warning",
            "Number of ART regiment before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m1$`m$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_a[ind,ind]))



#OFP
ind = c(1,2,6)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,6,7,8,9)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,7)
var = overall_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m1$`m$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

#YFPr
ind = c(2,3,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m1$`m$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_a[ind,ind]))


## young female pre vs young female post
ind = c(3,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,6)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,6)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m1$`m$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,6,8,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference is opposite)
ind = c(1,6,7,9)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m1$`m$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 5
beta = m1$`m$coefficients`[ind]
ade_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
ade_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
ade_point = exp(beta)
ade <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 6
beta = m1$`m$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_a[ind,ind]))

rr <- c(tb_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))

t.up_vl <-  data.frame(var_vl, rr,ll,ul,p)

##   Labels   ##
label(all$study_site.factor)   <- "Study site"
label(all$gender.factor)   <- "Gender"
label(all$age.new)   <- "Age at baseline"
label(all$age_gt50.factor)   <- "Age at baseline(chategorized) "
label(all$period) <- "Indicator for DTG warning: Pre/Post DTG warning"
label(all$viral_failure.n) <-  "Virological failure: >1000(Yes)"
label(all$ade_type) <- "History of other AIDS-Defining illness (not TB)"
label(all$tb) <- "TB Yes/No"
label(all$rna_v) <- "Viral load value"
label(all$event) <- "Event"



## Save files ##
save(all, file="all.Rdata")
save(t.up, file="t.up.Rdata")
save(t.up_a, file="t.up_a.Rdata")
save(t.up_imp, file="t.up_imp.Rdata")
save(t.up_vl, file="t.up_vl.Rdata")
save(t.up.interaction, file="t.up.interaction.Rdata")
save(t.up_i, file="t.up_i.Rdata")
save(all, file="all.Rdata")


### -------------------------- ###
### Separate analysis for site ###
### -------------------------- ###

## Since the data is highly dominated with observations from Haiti and stratifying by site assumes only
## same HR Performing time-updated survival analysis for Haiti and other sites separately  

all$haiti <- with(all,ifelse(study_site.factor == "Haiti", 1,0))
all_haiti <- all %>% filter(haiti ==1)

## Re-leveling ##
all_haiti <- within(all_haiti, age_gender_period <- relevel(factor(age_gender_period), ref = "Younger male post"))
all_haiti <- within(all_haiti, period <- relevel(factor(period), ref = "Post DTG warning"))
all_haiti <- within(all_haiti, gender.factor <- relevel(factor(gender.factor), ref = "Male"))
all_haiti <- within(all_haiti, age_gt50.factor <- relevel(factor(age_gt50.factor), ref = "Less then 50"))
all_haiti <- within(all_haiti, viral_failure.n <- relevel(factor(viral_failure.n), ref = "No"))
all_haiti <- within(all_haiti, num_cat.factor <- relevel(factor(num_cat.factor), ref = "One ART regiment before baseline"))
# There were no cases of ADE in Haiti therefore removing it from the RHS of analysis

## Model ##
m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor*gender.factor*period, data= all_haiti, cluster = patient_id)

overall_haiti_a <- data.matrix(m$var)
rownames(overall_haiti_a) <- NULL
colnames(overall_haiti_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL


#OMP
ind_test = 1
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))



#OFP
ind_test = c(1,2,4)
test_var = overall_haiti_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(1,2,3,4,5,6,7)
test_var = overall_haiti_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(1,3,5)
test_var = overall_haiti_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 2
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))

#YFPr
ind_test = c(2,3,6)
test_var = overall_haiti_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))


## young female pre vs young female post
ind_test = c(2,6)
test_var = overall_haiti_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


t.hup.interaction  <-  data.frame(variable, rr,ll,ul,p)


mi  <- coxph(Surv(time0, time1,event) ~ age_gender_period, data= all_haiti, cluster = patient_id)

overall_haiti_a <- data.matrix(mi$var)
rownames(overall_haiti_a) <- NULL
colnames(overall_haiti_a) <- NULL
m1 <- as.data.frame(mi$coefficients)
rownames(m1) <- NULL
#OMP

ind_test = 3
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))

#OFP
ind_test = 1
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))


#OFPr
ind_test = 2
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))



#OMPr
ind_test = 4
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))

#YFP
ind_test = 5
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))


#YFPr
ind_test = 6
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))

#YMPr
ind_test = 7
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overall_haiti_a[ind_test,ind_test]))

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3))


t.hup_i  <-  data.frame(variable, rr,ll,ul,p)


m2  <- coxph(Surv(time0, time1,event) ~  tb + num_cat.factor  + viral_failure.n  + age_gender_period, data = all_haiti, cluster = patient_id)

s <- summary(m2)
rownames(m1) <- NULL
s1 <- as.data.frame(s$conf.int)
a1 <- round(s$coefficients[51:60],digits=3)
colnames(a1) <- NULL

#Model summary table
variable <- c(" TB (ref = no)", "Number of ART regiment before baseline (ref = One)","Viral failure (ref = No)" , 
              "older female post (ref = younger male post)", "older female pre (ref = younger male post)", "older male post (ref = younger male post)", "older male pre (ref = younger male post)","younger female post (ref = younger male post)", "younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")


hr <- c(s1$`exp(coef)`[1],s1$`exp(coef)`[2],s1$`exp(coef)`[3],s1$`exp(coef)`[4],s1$`exp(coef)`[5],s1$`exp(coef)`[6],s1$`exp(coef)`[7],s1$`exp(coef)`[8],s1$`exp(coef)`[9],s1$`exp(coef)`[10])
ll <- c(s1$`lower .95`[1],s1$`lower .95`[2],s1$`lower .95`[3],s1$`lower .95`[4],s1$`lower .95`[5],s1$`lower .95`[6],s1$`lower .95`[7],s1$`lower .95`[8],s1$`lower .95`[9],s1$`lower .95`[10])
ul <- c(s1$`upper .95`[1],s1$`upper .95`[2],s1$`upper .95`[3],s1$`upper .95`[4],s1$`upper .95`[5],s1$`upper .95`[6],s1$`upper .95`[7],s1$`upper .95`[8],s1$`upper .95`[9],s1$`upper .95`[10])
p <- c(a1[1],a1[2],a1[3],a1[4],a1[5],a1[6],a1[7],a1[8],a1[9],a1[10])

t.hup <- data.frame(variable,hr,ll,ul,p)

m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb + viral_failure.n 
            + age_gt50.factor*gender.factor*period + num_cat ,
            data = all_haiti, cluster=patient_id) 

overall_haiti_a <- data.matrix(m$var)
rownames(overall_haiti_a) <- NULL
colnames(overall_haiti_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL



var_overall_haiti <- c("TB (ref=no)",
                       "Viral failure (ref= No)",
                       "older male post DTG warning (ref = younger male post)", "older female postDTG warning (ref = younger male post)",
                       "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
                       "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
                       "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                       "older female pre  & during DTG warning: older male pre  & during DTG warning",
                       "younger female pre & during DTG warning : older female pre & during DTG warning",
                       "Number of ART regiment before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m1$`m$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))



#OFP
ind = c(1,2,7)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,7,8,9,10)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,8)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m1$`m$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))

#YFPr
ind = c(2,3,9)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m1$`m$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))


## young female pre vs young female post
ind = c(3,9)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,7)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m1$`m$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m1$`m$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,9)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,9,10)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,7,8,10)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m1$`m$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))


ind = 5
beta = m1$`m$coefficients`[ind]
vl_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
vl_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
vl_point = exp(beta)
vl <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))

ind = 6
beta = m1$`m$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))

rr <- c(tb_point,vl_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,vl_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,vl_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(vl))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))


t.hup_a <-  data.frame(var_overall_haiti, rr,ll,ul,p)


## Imputed model 

set.seed(7)
cph_multi_i <- aregImpute(~ age_gt50.factor + gender.factor + period + tb  + viral_failure.n 
                          + age_gt50.factor*gender.factor*period + num_cat.factor + time1*event,
                          data = all_haiti, n.impute=20)

cph_multi <- fit.mult.impute(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb  + viral_failure.n 
                             + age_gt50.factor*gender.factor*period + num_cat.factor, cluster = all_haiti$patient_id,
                             data = all_haiti, fitter=cph, fitargs=list(x=TRUE, y=TRUE,surv=TRUE), xtrans=cph_multi_i, n.impute=20)

overall_haiti_a <- data.matrix(cph_multi$var)
rownames(overall_haiti_a) <- NULL
colnames(overall_haiti_a) <- NULL
m3 <- as.data.frame(cph_multi$coefficients)
rownames(m3) <- NULL



var_overall_haiti <- c("TB (ref=no)",
                       "Viral failure (ref= No)",
                       "Younger male post DTG warning (ref = Older male post)", "older female postDTG warning (ref = younger male post)",
                       "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
                       "young female pre : young female post DTG warning", "Younger Female post DTG warning : older female post DTG warning",
                       "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                       "older female pre  & during DTG warning: older male pre  & during DTG warning",
                       "younger female pre & during DTG warning : older female pre & during DTG warning",
                       "Number of ART regiment before baseline (ref = One)","Younger male pre & during : Older male pre & during","Older female Pre: Older female Post","Older Male Pre: Older Female Post")

#OMP
ind = 1
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
YMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YMP_point = exp(beta*(-1))
z1 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))



#OFP
ind = c(1,2,7)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,7,8,9,10)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,8)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))

#YFPr
ind = c(2,3,9)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))


## young female pre vs young female post
ind = c(3,9)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#YFP vs OFP
ind = c(1,7)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m3$`cph_multi$coefficients`[ind])
YFP_OFP_UB = exp(beta + 1.96*sqrt(se))
YFP_OFP_LB = exp(beta - 1.96*sqrt(se))
YFP_OFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m3$`cph_multi$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,9)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,9,10)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,7,8,10)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m3$`cph_multi$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))


ind = 5
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
vl_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
vl_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
vl_point = exp(beta)
vl <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))

ind = 6
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))

ind = c(1,8)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m3$`cph_multi$coefficients`[ind])
YMPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
YMPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
YMPr_OMPr_point = exp(beta)
z14 <- beta/sqrt(se)

ind = c(3,8,9,10)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_OFP_UB = exp(beta + 1.96*sqrt(se))
OFPr_OFP_LB = exp(beta - 1.96*sqrt(se))
OFPr_OFP_point = exp(beta)
z15 <- beta/sqrt(se)

ind = c(3,8)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OMPr_OMP_UB = exp(beta + 1.96*sqrt(se))
OMPr_OMP_LB = exp(beta - 1.96*sqrt(se))
OMPr_OMP_point = exp(beta)
z16 <- beta/sqrt(se)

rr <- c(tb_point,vl_point,YMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,YFP_OFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point,YMPr_OMPr_point,OFPr_OFP_point,OMPr_OMP_point)

ll <- c(tb_LB,vl_LB,YMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,YFP_OFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB,YMPr_OMPr_LB,OFPr_OFP_LB,OMPr_OMP_LB)
ul <- c(tb_UB,vl_UB,YMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,YFP_OFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB,YMPr_OMPr_UB,OFPr_OFP_UB,OMPr_OMP_UB)


p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(vl))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3),
       round((1 - pnorm(abs(z14))) * 2,digits=3),
       round((1 - pnorm(abs(z15))) * 2,digits=3),
       round((1 - pnorm(abs(z16))) * 2,digits=3))


t.hup_imp <-  data.frame(var_overall_haiti, rr,ll,ul,p)


### Sensitivity analysis

m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb   
            + age_gt50.factor*gender.factor*period + num_cat , 
            data = all_haiti, cluster= patient_id) 

overall_haiti_a <- data.matrix(m$var)
rownames(overall_haiti_a) <- NULL
colnames(overall_haiti_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL



var_vl <- c("TB (ref=no)",
            "older male post DTG warning (ref = younger male post)", "older female postDTG warning (ref = younger male post)",
            "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning",
            "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
            "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
            "older female pre  & during DTG warning: older male pre  & during DTG warning",
            "younger female pre & during DTG warning : older female pre & during DTG warning",
            "Number of ART regiment before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m1$`m$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))



#OFP
ind = c(1,2,6)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,6,7,8,9)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,7)
var = overall_haiti_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m1$`m$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))

#YFPr
ind = c(2,3,8)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m1$`m$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))


## young female pre vs young female post
ind = c(3,8)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,6)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,6)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m1$`m$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,8)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,6,8,9)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,6,7,9)
var = as.numeric(overall_haiti_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m1$`m$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))


ind = 5
beta = m1$`m$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_haiti_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_haiti_a[ind,ind]))

rr <- c(tb_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))

t.hup_vl <-  data.frame(var_vl, rr,ll,ul,p)

##   Labels   ##
label(all_haiti$study_site.factor)   <- "Study site"
label(all_haiti$gender.factor)   <- "Gender"
label(all_haiti$age.new)   <- "Age at baseline"
label(all_haiti$age_gt50.factor)   <- "Age at baseline(chategorized) "
label(all_haiti$period) <- "Indicator for DTG warning: Pre/Post DTG warning"
label(all_haiti$viral_failure.n) <-  "Virological failure: >1000(Yes)"
label(all_haiti$ade_type) <- "History of other AIDS-Defining illness (not TB)"
label(all_haiti$tb) <- "TB Yes/No"
label(all_haiti$rna_v) <- "Viral load value"
label(all_haiti$event) <- "Event"



## Other Sites (This long nested ifelse bacause study_site is a factor filtering was not working properly)
all$study_site <- with(all, ifelse(study_site.factor == "Haiti",1,
                                   ifelse(study_site.factor == "Brazil",2,
                                          ifelse(study_site.factor == "Chile",3,4))))

all_other <- all %>% filter(study_site != 1)
all_other$study_site.factor <- factor(all_other$study_site,
                                      levels=c(2,3,4),
                                      labels=c("Brazil", "Chile","Honduras" ))

## Model ##
m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor*gender.factor*period + strata(study_site.factor), data= all_other, cluster = patient_id)

overall_other_a <- data.matrix(m$var)
rownames(overall_other_a) <- NULL
colnames(overall_other_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL


#OMP
ind_test = 1
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))



#OFP
ind_test = c(1,2,4)
test_var = overall_other_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(1,2,3,4,5,6,7)
test_var = overall_other_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(1,3,5)
test_var = overall_other_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 2
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))

#YFPr
ind_test = c(2,3,6)
test_var = overall_other_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))


## young female pre vs young female post
ind_test = c(2,6)
test_var = overall_other_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


t.other.up.interaction  <-  data.frame(variable, rr,ll,ul,p)


mi  <- coxph(Surv(time0, time1,event) ~ age_gender_period + strata(study_site.factor), data= all_other, cluster = patient_id)

overall_other_a <- data.matrix(mi$var)
rownames(overall_other_a) <- NULL
colnames(overall_other_a) <- NULL
m1 <- as.data.frame(mi$coefficients)
rownames(m1) <- NULL
#OMP

ind_test = 3
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))

#OFP
ind_test = 1
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))


#OFPr
ind_test = 2
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))



#OMPr
ind_test = 4
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))

#YFP
ind_test = 5
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))


#YFPr
ind_test = 6
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))

#YMPr
ind_test = 7
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overall_other_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overall_other_a[ind_test,ind_test]))

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3))


t.other.up_i  <-  data.frame(variable, rr,ll,ul,p)



m2  <- coxph(Surv(time0, time1,event) ~  tb  + num_cat.factor  + viral_failure.n  + age_gender_period + strata(study_site.factor), data = all_other, cluster = patient_id)

overall_a <- data.matrix(m2$var)
m1 <- as.data.frame(m2$coefficients)

s <- summary(m2)
rownames(m1) <- NULL
s1 <- as.data.frame(s$conf.int)

#TB
ind = 1
beta = as.numeric(m1$`m2$coefficients`[ind])
z1 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 2
beta = as.numeric(m1$`m2$coefficients`[ind])
z2 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 3
beta = as.numeric(m1$`m2$coefficients`[ind])
z3 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 4
beta = as.numeric(m1$`m2$coefficients`[ind])
z4 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 5
beta = as.numeric(m1$`m2$coefficients`[ind])
z5 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 6
beta = as.numeric(m1$`m2$coefficients`[ind])
z6 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 7
beta = as.numeric(m1$`m2$coefficients`[ind])
z7 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 8
beta = as.numeric(m1$`m2$coefficients`[ind])
z8 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 9
beta = as.numeric(m1$`m2$coefficients`[ind])
z9 <- beta/sqrt(as.numeric(overall_a[ind,ind]))

ind = 10
beta = as.numeric(m1$`m2$coefficients`[ind])
z10 <- beta/sqrt(as.numeric(overall_a[ind,ind]))


#Model summary table
variable <- c(" TB (ref = no)", "Number of ART regiment before baseline (ref = One)","Viral failure (ref = No)" , 
              "older female post (ref = younger male post)", "older female pre (ref = younger male post)", "older male post (ref = younger male post)", "older male pre (ref = younger male post)","younger female post", "younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")


hr <- c(s1$`exp(coef)`[1],s1$`exp(coef)`[2],s1$`exp(coef)`[3],s1$`exp(coef)`[4],s1$`exp(coef)`[5],s1$`exp(coef)`[6],s1$`exp(coef)`[7],s1$`exp(coef)`[8],s1$`exp(coef)`[9],s1$`exp(coef)`[10])
ll <- c(s1$`lower .95`[1],s1$`lower .95`[2],s1$`lower .95`[3],s1$`lower .95`[4],s1$`lower .95`[5],s1$`lower .95`[6],s1$`lower .95`[7],s1$`lower .95`[8],s1$`lower .95`[9],s1$`lower .95`[10])
ul <- c(s1$`upper .95`[1],s1$`upper .95`[2],s1$`upper .95`[3],s1$`upper .95`[4],s1$`upper .95`[5],s1$`upper .95`[6],s1$`upper .95`[7],s1$`upper .95`[8],s1$`upper .95`[9],s1$`upper .95`[10])
p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3))

t.other.up <- data.frame(variable,hr,ll,ul,p)

m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb + viral_failure.n 
            + age_gt50.factor*gender.factor*period + num_cat.factor + strata(study_site.factor) ,
            data = all_other, cluster=patient_id) 

overall_other_a <- data.matrix(m$var)
rownames(overall_other_a) <- NULL
colnames(overall_other_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL



var_overall_other <- c("TB (ref=no)","Viral failure (ref= No)",
                       "older male post DTG warning (ref = younger male post) ", "older female postDTG warning (ref = younger male post)",
                       "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
                       "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
                       "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                       "older female pre  & during DTG warning: older male pre  & during DTG warning",
                       "younger female pre & during DTG warning : older female pre & during DTG warning",
                       "Number of ART regiment before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m1$`m$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))



#OFP
ind = c(1,2,7)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,7,8,9,10)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,8)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m1$`m$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

#YFPr
ind = c(2,3,9)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m1$`m$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))


## young female pre vs young female post
ind = c(3,9)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,7)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m1$`m$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,9)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,9,10)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,7,8,10)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m1$`m$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

# ind = 5
# beta = m1$`m$coefficients`[ind]
# ade_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
# ade_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
# ade_point = exp(beta)
# ade <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

ind = 5
beta = m1$`m$coefficients`[ind]
vl_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
vl_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
vl_point = exp(beta)
vl <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

ind = 6
beta = m1$`m$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

rr <- c(tb_point,vl_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,vl_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,vl_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(vl))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))


t.other.up_a <-  data.frame(var_overall_other, rr,ll,ul,p)


## Imputed model 

dd <- datadist(all_other)
set.seed(7)
cph_multi_i <- aregImpute(~ age_gt50.factor + gender.factor + period + tb  + viral_failure.n 
                          + strata(study_site.factor) + age_gt50.factor*gender.factor*period + num_cat.factor + time1*event,
                          data = all_other, n.impute=20)

cph_multi <- fit.mult.impute(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb  + viral_failure.n 
                             + strat(study_site.factor) + age_gt50.factor*gender.factor*period + num_cat.factor, cluster = all_other$patient_id,
                             data = all_other,fitter=cph, fitargs=list(x=TRUE, y=TRUE,surv=TRUE), xtrans=cph_multi_i, n.impute=20)

overall_other_a <- data.matrix(cph_multi$var)
rownames(overall_other_a) <- NULL
colnames(overall_other_a) <- NULL
m3 <- as.data.frame(cph_multi$coefficients)
rownames(m3) <- NULL



var_overall_other <- c("TB (ref=no)",
                       "Viral failure (ref= No)",
                       "Younger male post :older male post DTG warning", "older female postDTG warning",
                       "older female pre  & during DTG warning","older male pre  & during DTG warning","younger female post DTG warning","younger female pre  & during DTG warning","younger male pre  & during DTG warning",
                       "young female pre : young female post DTG warning", "Younger Female post DTG warning : older female post DTG warning",
                       "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
                       "older female pre  & during DTG warning: older male pre  & during DTG warning",
                       "younger female pre & during DTG warning : older female pre & during DTG warning",
                       "Number of ART regiment before baseline (ref = One)","Younger Male Pre: Older Male Pre"
                       ,"Older female Pre: Older female Post","Older Male Pre: Older Female Post" )

#YMP
ind = 1
beta = as.numeric(m3$`cph_multi$coefficients`[ind])
YMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YMP_point = exp(beta*(-1))
z1 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))



#OFP
ind = c(1,2,7)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,7,8,9,10)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,8)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m3$`cph_multi$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

#YFPr
ind = c(2,3,9)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m3$`cph_multi$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))


## young female pre vs young female post
ind = c(3,9)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,7)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m3$`cph_multi$coefficients`[ind])
YFP_OFP_UB = exp(beta + 1.96*sqrt(se))
YFP_OFP_LB = exp(beta - 1.96*sqrt(se))
YFP_OFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,7)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m3$`cph_multi$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,9)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,7,9,10)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,7,8,10)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m3$`cph_multi$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m3$`cph_multi$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

ind = 5
beta = m3$`cph_multi$coefficients`[ind]
vl_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
vl_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
vl_point = exp(beta)
vl <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

ind = 6
beta = m3$`cph_multi$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

ind = c(1,8)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m3$`cph_multi$coefficients`[ind])
YMPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
YMPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
YMPr_OMPr_point = exp(beta)
z14 <- beta/sqrt(se)

ind = c(3,8,9,10)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OFPr_OFP_UB = exp(beta + 1.96*sqrt(se))
OFPr_OFP_LB = exp(beta - 1.96*sqrt(se))
OFPr_OFP_point = exp(beta)
z15 <- beta/sqrt(se)

ind = c(3,8)
var = as.numeric(overall_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m3$`cph_multi$coefficients`[ind])
OMPr_OMP_UB = exp(beta + 1.96*sqrt(se))
OMPr_OMP_LB = exp(beta - 1.96*sqrt(se))
OMPr_OMP_point = exp(beta)
z16 <- beta/sqrt(se)

rr <- c(tb_point,vl_point,YMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,YFP_OFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point,YMPr_OMPr_point,OFPr_OFP_point,OMPr_OMP_point)

ll <- c(tb_LB,vl_LB,YMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,YFP_OFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB,YMPr_OMPr_LB,OFPr_OFP_LB,OMPr_OMP_LB)
ul <- c(tb_UB,vl_UB,YMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,YFP_OFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB,YMPr_OMPr_UB,OFPr_OFP_UB,OMPr_OMP_UB)


p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(vl))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3),
       round((1 - pnorm(abs(z14))) * 2,digits=3),
       round((1 - pnorm(abs(z15))) * 2,digits=3),
       round((1 - pnorm(abs(z16))) * 2,digits=3))


t.other.up_imp <-  data.frame(var_overall_other, rr,ll,ul,p)


### Sensitivity analysis

m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor + gender.factor + period + tb   
            + strata(study_site.factor) + age_gt50.factor*gender.factor*period + num_cat.factor , 
            data = all_other, cluster= patient_id) 

overall_other_a <- data.matrix(m$var)
rownames(overall_other_a) <- NULL
colnames(overall_other_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL



var_vl <- c("TB (ref=no)",
            "older male post DTG warning (ref = younger male post)", "older female postDTG warning (ref = younger male post)",
            "older female pre  & during DTG warning (ref = younger male post)","older male pre  & during DTG warning (ref = younger male post)","younger female post DTG warning (ref = younger male post)","younger female pre  & during DTG warning (ref = younger male post)","younger male pre  & during DTG warning (ref = younger male post)",
            "young female pre : young female post DTG warning", "older Female post DTG warning : younger female post DTG warning",
            "older male post DTG warning : older female post DTG warning", "younger female pre & during DTG warning: younger male pre & during DTG warning ",
            "older female pre  & during DTG warning: older male pre  & during DTG warning",
            "younger female pre & during DTG warning : older female pre & during DTG warning",
            "Number of ART regiment before baseline (ref = One)" )

#OMP
ind = 1
beta = as.numeric(m1$`m$coefficients`[ind])
OMP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
OMP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
OMP_point = exp(beta)
z1 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))



#OFP
ind = c(1,2,6)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_UB = exp(beta + 1.96*sqrt(se))
OFP_LB = exp(beta - 1.96*sqrt(se))
OFP_point = exp(beta)
z2 <- beta/sqrt(se)


#OFPr
ind = c(1,2,3,6,7,8,9)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_point = exp(beta)
z3 <- beta/sqrt(se)


#OMPr
ind = c(1,3,7)
var = overall_other_a[ind, ind]
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])  
OMPr_UB = exp(beta + 1.96*sqrt(se))
OMPr_LB = exp(beta - 1.96*sqrt(se))
OMPr_point = exp(beta)
z4 <- beta/sqrt(se)

#YFP
ind = 2
beta = m1$`m$coefficients`[ind]
YFP_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YFP_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YFP_point = exp(beta)
z5 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

#YFPr
ind = c(2,3,8)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_point = exp(beta)
z6 <- beta/sqrt(se)

#YMPr
ind = 3
beta = m1$`m$coefficients`[ind]
YMPr_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YMPr_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
YMPr_point = exp(beta)
z7 <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))


## young female pre vs young female post
ind = c(3,8)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YFP_UB = exp(beta + 1.96*sqrt(se))
YFPr_YFP_LB = exp(beta - 1.96*sqrt(se))
YFPr_YFP_point = exp(beta)
YFPr_YFP_LB
z8 <- beta/sqrt(se)

#OFP vs YFP
ind = c(1,6)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFP_YFP_UB = exp(beta + 1.96*sqrt(se))
OFP_YFP_LB = exp(beta - 1.96*sqrt(se))
OFP_YFP_point = exp(beta)
z9 <- beta/sqrt(se)

#OMP vs OFP
ind = c(2,6)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = (-1)*sum(m1$`m$coefficients`[ind])
OMP_OFP_UB = exp(beta + 1.96*sqrt(se))
OMP_OFP_LB = exp(beta - 1.96*sqrt(se))
OMP_OFP_point = exp(beta)
z10 <- beta/sqrt(se)



#YFPr vs YMPr
ind = c(2,8)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
YFPr_YMPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_YMPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_YMPr_point = exp(beta)
z11 <- beta/sqrt(se)


#OFPr_OMPr
ind = c(2,6,8,9)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum(m1$`m$coefficients`[ind])
OFPr_OMPr_UB = exp(beta + 1.96*sqrt(se))
OFPr_OMPr_LB = exp(beta - 1.96*sqrt(se))
OFPr_OMPr_point = exp(beta)
z12 <- beta/sqrt(se)


#OFPr YFPr (-1 because the reference id opposite)
ind = c(1,6,7,9)
var = as.numeric(overall_other_a[ind, ind])
se = sum(diag(var)) + 2*sum(var[upper.tri(var, diag = FALSE)])
beta = sum((-1)*m1$`m$coefficients`[ind])
YFPr_OFPr_UB = exp(beta + 1.96*sqrt(se))
YFPr_OFPr_LB = exp(beta - 1.96*sqrt(se))
YFPr_OFPr_point = exp(beta)
z13 <- beta/sqrt(se)

ind = 4
beta = m1$`m$coefficients`[ind]
tb_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
tb_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
tb_point = exp(beta)
tb <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

# ind = 5
# beta = m1$`m$coefficients`[ind]
# ade_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
# ade_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
# ade_point = exp(beta)
# ade <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

ind = 5
beta = m1$`m$coefficients`[ind]
rg_UB = exp(beta + 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
rg_LB = exp(beta - 1.96*sqrt(as.numeric(overall_other_a[ind,ind])))
rg_point = exp(beta)
rg <- beta/sqrt(as.numeric(overall_other_a[ind,ind]))

rr <- c(tb_point,OMP_point,OFP_point,OFPr_point,OMPr_point,YFP_point,YFPr_point,
        YMPr_point,YFPr_YFP_point,OFP_YFP_point,OMP_OFP_point,YFPr_YMPr_point,OFPr_OMPr_point,YFPr_OFPr_point,rg_point)

ll <- c(tb_LB,OMP_LB,OFP_LB,OFPr_LB,OMPr_LB,YFP_LB,YFPr_LB,
        YMPr_LB,YFPr_YFP_LB,OFP_YFP_LB,OMP_OFP_LB,YFPr_YMPr_LB,OFPr_OMPr_LB, YFPr_OFPr_LB,rg_LB)
ul <- c(tb_UB,OMP_UB,OFP_UB,OFPr_UB,OMPr_UB,YFP_UB,YFPr_UB,
        YMPr_UB,YFPr_YFP_UB,OFP_YFP_UB,OMP_OFP_UB,YFPr_YMPr_UB,OFPr_OMPr_UB, YFPr_OFPr_UB,rg_UB)

p <- c(round((1 - pnorm(abs(tb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),
       round((1 - pnorm(abs(z2))) * 2,digits =3),
       round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),
       round((1 - pnorm(abs(z5))) * 2,digits=3),
       round((1 - pnorm(abs(z6))) * 2,digits =3),
       round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits=3),
       round((1 - pnorm(abs(z9))) * 2,digits=3),
       round((1 - pnorm(abs(z10))) * 2,digits =3),
       round((1 - pnorm(abs(z11))) * 2,digits=3),
       round((1 - pnorm(abs(z12))) * 2,digits=3),
       round((1 - pnorm(abs(z13))) * 2,digits=3),
       round((1 - pnorm(abs(rg))) * 2,digits=3))

t.other.up_vl <-  data.frame(var_vl, rr,ll,ul,p)


## Save files ##
save(all_other, file="all_other.Rdata")
save(t.other.up, file="t.other.up.Rdata")
save(t.other.up_a, file="t.other.up_a.Rdata")
save(t.other.up_imp, file="t.other.up_imp.Rdata")
save(t.other.up_vl, file="t.other.up_vl.Rdata")
save(t.other.up.interaction, file="t.other.up.interaction.Rdata")
save(t.other.up_i, file="t.other.up_i.Rdata")
save(all_other, file="all_other.Rdata")

## Males only 


data_m <- demdata_1 %>% filter(gender.factor == "Male")
dem_male <- demdata_1 %>% filter(gender.factor == "Male")




data_m$age_period.factor <- with(data_m,ifelse(age_gt50.factor == "Greater then or equal to 50" & DTG_pre.factor == "Pre and during DTG", "older pre",
                                               ifelse(age_gt50.factor == "Less then 50"                & DTG_pre.factor == "Pre and during DTG", "younger pre",
                                                      ifelse(age_gt50.factor == "Greater then or equal to 50" & DTG_pre.factor == "Post DTG"           , "older post",
                                                             ifelse(age_gt50.factor == "Less then 50" & DTG_pre.factor == "Post DTG"  ,  "younger post",
                                                                    "99")))))     

#data_m$DTG.factor <- ifelse(data_m$DTG.factor == "No", 0, 1)

dd <- datadist(data_m)
options(datadist="dd")
dd$limits$age_period.factor <- "younger post"
dd$limits$tb <- "No"
dd$limits$ade_type <- "No"

# Setting levels  
data_m <- within(data_m, age_period.factor <- relevel(factor(age_period.factor), ref = "younger post"))
data_m <- within(data_m, tb <- relevel(factor(tb), ref = "No"))
data_m <- within(data_m, ade_type <- relevel(factor(ade_type), ref = "No"))




# Creating data frame
unmm <- data.frame(Variable= c(c("Age (ref =Less then 50)", "Greater then or equal to 50"),
                               c("Site (ref= Haiti)", "Brazil", "Chile", "Honduras"),
                               c("DTG: post vs. Pre and during DTG (ref = Post DTG)","Pre and during DTG"),
                               c("Other AIDS defining illness (ref =No)", "Other AIDS defining illness: Yes"),
                               c("TB (ref = No)", "TB: Yes")),
                   RR= NA, Lower=NA, Upper=NA,  p=NA)

# un-adjusted models
#age_gt50
m1_un_age_gt50 <- glm(DTG ~   age_gt50.factor, data = data_m, family=poisson(link="log"))
m1_un_Age_gt50 <- Glm(DTG ~   age_gt50.factor, data = data_m, family=poisson(link="log"))
summ_un_age_gt50 <- as.data.frame(my_summary.rms(object=m1_un_Age_gt50, object2 = m1_un_age_gt50))
exposures <- " "
anova_un_age_gt50 <- as.data.frame(anova_w_sandwich(object = m1_un_Age_gt50, objectglm = m1_un_age_gt50,coVars_woInt = c("age_gt50")))

#TB
m1_un_tb <- glm(DTG ~   tb, data = data_m, family=poisson(link="log"))
m1_un_Tb <- Glm(DTG ~   tb, data = data_m, family=poisson(link="log"))
summ_un_tb <- as.data.frame(my_summary.rms(object=m1_un_Tb, object2 = m1_un_tb))
exposures <- " "
anova_un_tb <- as.data.frame(anova_w_sandwich(object = m1_un_Tb, objectglm = m1_un_tb,coVars_woInt = c("tb")))

#ADE
m1_un_ade <- glm(DTG ~   ade_type, data = data_m, family=poisson(link="log"))
m1_un_Ade <- Glm(DTG ~   ade_type, data = data_m, family=poisson(link="log"))
summ_un_ade <- as.data.frame(my_summary.rms(object=m1_un_Ade, object2 = m1_un_ade))
exposures <- " "
anova_un_ade <- as.data.frame(anova_w_sandwich(object = m1_un_Ade, objectglm = m1_un_ade,coVars_woInt = c("ade_type")))


#site
m1_un_site <- glm(DTG ~   study_site.factor, data = data_m, family=poisson(link="log"))
m1_un_Site <- Glm(DTG ~   study_site.factor, data = data_m, family=poisson(link="log"))
summ_un_site <- as.data.frame(my_summary.rms(object=m1_un_Site, object2 = m1_un_site))
exposures <- " "
anova_un_site <- as.data.frame(anova_w_sandwich(object = m1_un_Site, objectglm = m1_un_site,coVars_woInt = c("study_site.factor")))


#DTG
m1_un_dtg <- glm(DTG ~   DTG_pre.factor, data = data_m, family=poisson(link="log"))
m1_un_Dtg <- Glm(DTG ~   DTG_pre.factor, data = data_m, family=poisson(link="log"))
summ_un_dtg <- as.data.frame(my_summary.rms(object=m1_un_Dtg, object2 = m1_un_dtg))
exposures <- " "
anova_un_dtg <- as.data.frame(anova_w_sandwich(object = m1_un_Dtg, objectglm = m1_un_dtg,coVars_woInt = c("DTG_pre")))


##Unadjusted Model
unmm$RR[unmm$Variable %in% "Pre and during DTG"] <- round(exp(summ_un_dtg$Effect[1]), digits =3)
unmm$RR[unmm$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(summ_un_ade$Effect[1]), digits =3)
unmm$RR[unmm$Variable %in% "TB: Yes"] <- round(exp(summ_un_tb$Effect[1]), digits =3)
unmm$RR[unmm$Variable %in% "Brazil"] <- round(exp(summ_un_site$Effect[1]), digits =3)
unmm$RR[unmm$Variable %in% "Chile"] <- round(exp(summ_un_site$Effect[2]), digits =3)
unmm$RR[unmm$Variable %in% "Honduras"] <- round(exp(summ_un_site$Effect[3]), digits =3)
unmm$RR[unmm$Variable %in% "Greater then or equal to 50"] <- round(exp(summ_un_age_gt50$Effect[1]), digits =3)
unmm$Lower[unmm$Variable %in% "Pre and during DTG"] <- round(exp(summ_un_dtg$`Lower 0.95`[1]), digits =3)
unmm$Lower[unmm$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(summ_un_ade$`Lower 0.95`[1]), digits =3)
unmm$Lower[unmm$Variable %in% "TB: Yes"] <- round(exp(summ_un_tb$`Lower 0.95`[1]), digits =3)
unmm$Lower[unmm$Variable %in% "Brazil"] <- round(exp(summ_un_site$`Lower 0.95`[1]), digits =3)
unmm$Lower[unmm$Variable %in% "Chile"] <- round(exp(summ_un_site$`Lower 0.95`[2]), digits =3)
unmm$Lower[unmm$Variable %in% "Honduras"] <- round(exp(summ_un_site$`Lower 0.95`[3]), digits =3)
unmm$Lower[unmm$Variable %in% "Greater then or equal to 50"] <- round(exp(summ_un_age_gt50$`Lower 0.95`[1]), digits =3)
unmm$Upper[unmm$Variable %in% "Pre and during DTG"] <- round(exp(summ_un_dtg$`Upper 0.95`[1]), digits =3)
unmm$Upper[unmm$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(summ_un_ade$`Upper 0.95`[1]), digits =3)
unmm$Upper[unmm$Variable %in% "TB: Yes"] <- round(exp(summ_un_tb$`Upper 0.95`[1]), digits =3)
unmm$Upper[unmm$Variable %in% "Brazil"] <- round(exp(summ_un_site$`Upper 0.95`[1]), digits =3)
unmm$Upper[unmm$Variable %in% "Chile"] <- round(exp(summ_un_site$`Upper 0.95`[2]), digits =3)
unmm$Upper[unmm$Variable %in% "Honduras"] <- round(exp(summ_un_site$`Upper 0.95`[3]), digits =3)
unmm$Upper[unmm$Variable %in% "Greater then or equal to 50"] <- round(exp(summ_un_age_gt50$`Upper 0.95`[1]), digits =3)

# p values for unadjusted table
unmm$p[unmm$Variable %in% "DTG: post vs. Pre and during DTG (ref = Post DTG)"] <- anova_un_dtg$stats[1]
unmm$p[unmm$Variable %in% "Other AIDS defining illness (ref =No)"] <- anova_un_ade$stats[1]
unmm$p[unmm$Variable %in% "TB (ref = No)"] <- anova_un_tb$stats[1]
unmm$p[unmm$Variable %in% "Site (ref= Haiti)"] <- anova_un_site$stats[1]
unmm$p[unmm$Variable %in% "Age (ref = Less then 50)"] <- anova_un_age_gt50$stats[1]


#Interaction model dummy variables

mod_dummy_interaction <- glm(DTG ~   age_period.factor, data = data_m, family=poisson(link="log"))
mod_dummy_Interaction <- Glm(DTG ~   age_period.factor, data = data_m, family=poisson(link="log"))
exposures <- " "
summ_int <- as.data.frame(my_summary.rms(object=mod_dummy_Interaction, object2 = mod_dummy_interaction))
var <- sandwich(object2)[-1, -1] # Removes Intercept row/column
object2 = mod_dummy_interaction
z1 <- summ_int$Effect[1]/summ_int$S.E.[1]
z2 <- summ_int$Effect[2]/summ_int$S.E.[2]
z3 <- summ_int$Effect[3]/summ_int$S.E.[3]

summ_int$p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3))
summ_int$Effect <- round(exp(summ_int$Effect), digits =3)
summ_int$`Lower 0.95` <- round(exp(summ_int$`Lower 0.95`), digits =3)
summ_int$`Upper 0.95` <- round(exp(summ_int$`Upper 0.95`), digits =3)
variables <- c("older post:younger post","older pre:younger post","younger pre:younger post")
int_m <- as.data.frame(cbind(variables,summ_int$Effect,summ_int$`Lower 0.95`,summ_int$`Upper 0.95`,summ_int$p))
colnames(int_m) <- c("variables","Effect","Lower","Upper","p-value")



## Interaction terms model
mod_interaction <- glm(DTG ~    age_gt50.factor*DTG_pre.factor , data= data_m, family=poisson(link="log"))
summ_int <- summary(mod_interaction)
var <- sandwich(mod_interaction)

#OP
ind_test = 2
test_beta = summ_int$coefficients[ind_test]
OP_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
OP_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
OP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(var[ind_test,ind_test])


#OPre
ind_test = c(2,3,4)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OPre_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPre_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPre_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#YPre
ind_test = 3
test_beta = summ_int$coefficients[ind_test]
YPre_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
YPre_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
YPre_point_test = exp(test_beta)
z3 <- test_beta/sqrt(var[ind_test,ind_test])


#older pre vs Older Post
#OPre
ind_test = c(3,4)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OPOST_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPOST_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPOST_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

## Interaction term
variable <- variables <- c("older post:younger post","older pre:younger post","younger pre:younger post","older pre : older Post")
rr <- c(OP_point_test,OPre_point_test,YPre_point_test,OPOST_point_test)
ll <- c(OP_LB_test,OPre_LB_test,YPre_LB_test,OPOST_LB_test)
ul <- c(OP_UB_test,OPre_UB_test,YPre_UB_test,OPOST_UB_test)
p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3))


int_m_cat <-  data.frame(variable, rr,ll,ul,p)


#Dummy Overall term
mod_dummy_a <- glm(DTG ~    study_site.factor + tb + ade_type + age_period.factor, data= data_m, family=poisson(link="log"))
mod_Dummy_a <- Glm(DTG ~  study_site.factor  + tb + ade_type + age_period.factor , data= data_m, family=poisson(link="log"))
summ_dummy_a <- as.data.frame(my_summary.rms(object = mod_Dummy_a, object2 = mod_dummy_a))
var_dummy_a <- sandwich(mod_dummy_a)

# Table for adjusted model with dummy variables
variable <- c("Brazil","Chile","Honduras","Tb (ref =no)","Other AIDS defining illness (ref =no)","older post:younger post", "older pre:younger post","younger pre:younger post")
rr <- c(round(exp(summ_dummy_a$Effect[1]), digits =3),round(exp(summ_dummy_a$Effect[2]), digits =3),round(exp(summ_dummy_a$Effect[3]), digits =3)
        ,round(exp(summ_dummy_a$Effect[4]), digits =3),round(exp(summ_dummy_a$Effect[5]), digits =3),round(exp(summ_dummy_a$Effect[6]), digits =3)
        ,round(exp(summ_dummy_a$Effect[7]), digits =3),round(exp(summ_dummy_a$Effect[8]), digits =3))

ll <- c(round(exp(summ_dummy_a$`Lower 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[8]), digits =3))

ul <- c(round(exp(summ_dummy_a$`Upper 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[8]), digits =3))                                                                                                                                                                                   

z1 <- summ_dummy_a$Effect[1]/summ_dummy_a$S.E.[1]
z2 <- summ_dummy_a$Effect[2]/summ_dummy_a$S.E.[2]
z3 <- summ_dummy_a$Effect[3]/summ_dummy_a$S.E.[3]
z4 <- summ_dummy_a$Effect[4]/summ_dummy_a$S.E.[4]
z5 <- summ_dummy_a$Effect[5]/summ_dummy_a$S.E.[5]
z6 <- summ_dummy_a$Effect[6]/summ_dummy_a$S.E.[6]
z7 <- summ_dummy_a$Effect[7]/summ_dummy_a$S.E.[7]
z8 <- summ_dummy_a$Effect[8]/summ_dummy_a$S.E.[8]



p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits = 3))
overall_mm_dummy <- data.frame(variable,rr,ll,ul,p)


## Overall terms model
oafm <- glm(DTG ~  age_gt50.factor + tb + ade_type + study_site.factor + DTG_pre.factor +  age_gt50.factor*DTG_pre.factor , data= data_m, family=poisson(link="log"))
summ_overall_a <- summary(oafm)
var_overall_a <- sandwich(oafm)


ind_atest = 3
atest_beta = summ_overall_a$coefficients[ind_atest]
tb_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
tb_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
tb_point_atest = exp(atest_beta)
aztb <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])

ind_atest = 4
atest_beta = summ_overall_a$coefficients[ind_atest]
ade_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
ade_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
ade_point_atest = exp(atest_beta)
azade <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])




ind_atest = 5
atest_beta = summ_overall_a$coefficients[ind_atest]
bra_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
bra_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
bra_point_atest = exp(atest_beta)
azbra <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])

ind_atest = 6
atest_beta = summ_overall_a$coefficients[ind_atest]
chile_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
chile_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
chile_point_atest = exp(atest_beta)
azchile <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])

ind_atest = 7
atest_beta = summ_overall_a$coefficients[ind_atest]
hon_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
hon_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
hon_point_atest = exp(atest_beta)
azhon <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])

#OP
ind_test = 2
test_beta = summ_overall_a$coefficients[ind_test]
OP_UB_test = exp(test_beta + 1.96*sqrt(var_overall_a[ind_test,ind_test]))
OP_LB_test = exp(test_beta - 1.96*sqrt(var_overall_a[ind_test,ind_test]))
OP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(var_overall_a[ind_test,ind_test])


#OPre
ind_test = c(2,8,9)
test_var_overall_a = var_overall_a[ind_test, ind_test]
se_test = sum(diag(test_var_overall_a)) + 2*sum(test_var_overall_a[upper.tri(test_var_overall_a, diag = FALSE)])
test_beta = sum(summ_overall_a$coefficients[ind_test])
OPre_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPre_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPre_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#YPre
ind_test = 8
test_beta = summ_overall_a$coefficients[ind_test]
YPre_UB_test = exp(test_beta + 1.96*sqrt(var_overall_a[ind_test,ind_test]))
YPre_LB_test = exp(test_beta - 1.96*sqrt(var_overall_a[ind_test,ind_test]))
YPre_point_test = exp(test_beta)
z3 <- test_beta/sqrt(var_overall_a[ind_test,ind_test])


#older pre vs Older Post
#OPre
ind_test = c(8,9)
test_var_overall_a = var_overall_a[ind_test, ind_test]
se_test = sum(diag(test_var_overall_a)) + 2*sum(test_var_overall_a[upper.tri(test_var_overall_a, diag = FALSE)])
test_beta = sum(summ_overall_a$coefficients[ind_test])
OPOST_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPOST_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPOST_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)


#Younger Pre DTG vs. Older Post
ind_test = c(8,1)
test_var_overall_a = var_overall_a[ind_test, ind_test]
se_test = sum(diag(test_var_overall_a)) - 2*sum(test_var_overall_a[upper.tri(test_var_overall_a, diag = FALSE)])
test_beta = (summ_overall_a$coefficients[9] - summ_overall_a$coefficients[1])
YPRE_OP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YPRE_OP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YPRE_OP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(se_test)



#Older Pre DTG vs. Younger Pre
ind_test = c(2,9)
test_var_overall_a = var_overall_a[ind_test, ind_test]
se_test = sum(diag(test_var_overall_a)) + 2*sum(test_var_overall_a[upper.tri(test_var_overall_a, diag = FALSE)])
test_beta = sum(summ_overall_a$coefficients[ind_test])
OPRE_YPRE_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPRE_YPRE_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPRE_YPRE_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)



var_overall_aiable <- c("Brazil","Chile","Honduras","Other illness(ref=no)","TB (ref=no)","older post:younger post","older pre:younger post","younger pre:younger post","older pre : older Post", "Younger Pre: Older Post", "Older Pre : Younger Pre")
rr <- c(bra_point_atest,chile_point_atest,hon_point_atest,ade_point_atest,tb_point_atest,OP_point_test,OPre_point_test,YPre_point_test,OPOST_point_test,YPRE_OP_point_test, OPRE_YPRE_point_test)
ll <- c(bra_LB_atest,chile_LB_atest,hon_LB_atest,ade_LB_atest,tb_LB_atest,OP_LB_test,OPre_LB_test,YPre_LB_test,OPOST_LB_test,YPRE_OP_LB_test, OPRE_YPRE_LB_test)
ul <- c(bra_UB_atest,chile_UB_atest,hon_UB_atest,ade_UB_atest,tb_UB_atest,OP_UB_test,OPre_UB_test,YPre_UB_test,OPOST_UB_test,YPRE_OP_UB_test, OPRE_YPRE_UB_test)

p <- c(round((1 - pnorm(abs(azbra))) * 2,digits = 3),
       round((1 - pnorm(abs(azchile))) * 2,digits = 3),round((1 - pnorm(abs(azhon))) * 2,digits = 3),
       round((1 - pnorm(abs(azade))) * 2,digits = 3),round((1 - pnorm(abs(aztb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits =3),round((1 - pnorm(abs(z6))) * 2,digits =3))

overall_mm_cat <-  data.frame(var_overall_aiable, rr,ll,ul,p)





#model
fm <- glm(DTG ~ ns(age, df = 4) + study_site.factor +  DTG_pre.factor +  ade_type + tb,data=data_m,family=poisson(link="log"))


# Spline matrix
ages <- round(min(data_m$age)):round(max(data_m$age))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")


# Summary, coefficients and standard errors from sandwich estimates
fm_summ <- summary(fm, conf.int=TRUE)
covar <- sandwich(fm)
se <- sqrt(covar[row(covar)==col(covar)])
fm_summ <- as.data.frame(cbind(fm_summ$coefficients[,1], se))
fm_summ$vars <- rownames(fm_summ)
colnames(fm_summ) <- c("Coefs", "SE", "Variables")

# Creating datat frame
fm_exp_df <- data.frame(Variable= c(c("Age (ref =35)", 25,45),
                                    c("Site (ref= Haiti)", "Brazil", "Chile", "Honduras"),
                                    c("DTG: post vs. Pre and during DTG (ref = Post DTG)","Pre and during DTG"),
                                    c("Other AIDS defining illness (ref =No)", "Other AIDS defining illness: Yes"),
                                    c("TB (ref = No)", "TB: Yes")),
                        RR= NA, Lower=NA, Upper=NA,  p=NA)



# creating RR and Limits
fm_exp_df$RR[fm_exp_df$Variable %in% "Pre and during DTG"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'ade_typeYes']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "TB: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'tbYes']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "Brazil"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "Chile"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorChile']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "Honduras"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% 25] <- exp(fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)1"]*
                                                  (sp[rownames(sp) %in% 25, 1]- sp[rownames(sp) %in% 35, 1]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)2"]*
                                                  (sp[rownames(sp) %in% 25, 2]- sp[rownames(sp) %in% 35, 2]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)3"]*
                                                  (sp[rownames(sp) %in% 25, 3]- sp[rownames(sp) %in% 35, 3]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)4"]*
                                                  (sp[rownames(sp) %in% 25, 4]-sp[rownames(sp) %in% 35, 4]))

fm_exp_df$RR[fm_exp_df$Variable %in% 45] <- exp(fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)1"]*
                                                  (sp[rownames(sp) %in% 45, 1]-sp[rownames(sp) %in% 35, 1]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)2"]*
                                                  (sp[rownames(sp) %in% 45, 2]- sp[rownames(sp) %in% 35, 2]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)3"]*
                                                  (sp[rownames(sp) %in% 45, 3]- sp[rownames(sp) %in% 35, 3]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)4"]*
                                                  (sp[rownames(sp) %in% 45, 4]- sp[rownames(sp) %in% 35, 4]))

# my function to calculate variance and SE
spline_se_fun <- function(myage){
  
  tmp <- covar
  a <- sp$X1[rownames(sp) %in% myage]-sp$X1[rownames(sp) %in% 35]
  b <- sp$X2[rownames(sp) %in% myage]-sp$X2[rownames(sp) %in% 35]
  c <- sp$X3[rownames(sp) %in% myage]-sp$X3[rownames(sp) %in% 35]
  d <- sp$X4[rownames(sp) %in% myage]-sp$X4[rownames(sp) %in% 35]
  
  Var_x <- tmp["ns(age, df = 4)1", "ns(age, df = 4)1"]
  Var_y <- tmp["ns(age, df = 4)2", "ns(age, df = 4)2"]
  Var_z <- tmp["ns(age, df = 4)3", "ns(age, df = 4)3"]
  Var_w <- tmp["ns(age, df = 4)4", "ns(age, df = 4)4"]
  
  Cov_xy <- tmp["ns(age, df = 4)1", "ns(age, df = 4)2"]
  Cov_xz <- tmp["ns(age, df = 4)1", "ns(age, df = 4)3"]
  Cov_xw <- tmp["ns(age, df = 4)1", "ns(age, df = 4)4"]
  Cov_yz <- tmp["ns(age, df = 4)2", "ns(age, df = 4)3"]
  Cov_yw <- tmp["ns(age, df = 4)2", "ns(age, df = 4)4"]
  Cov_zw <- tmp["ns(age, df = 4)3", "ns(age, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

SE_25 <- spline_se_fun(25)
SE_45 <-spline_se_fun(45)
# SE with sandwich estimates lower limits
fm_exp_df$Lower[fm_exp_df$Variable %in% "Pre and during DTG"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG'] - 
                                                                             1.96*fm_summ$SE[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'ade_typeYes'] - 
                                                                                           1.96*fm_summ$SE[fm_summ$Variables %in% 'ade_typeYes']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% "TB: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'tbYes'] - 
                                                                  1.96*fm_summ$SE[fm_summ$Variables %in% 'tbYes']), digits =3)

fm_exp_df$Lower[fm_exp_df$Variable %in% "Brazil"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorBrazil'] - 
                                                                 1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% "Chile"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorChile'] - 
                                                                1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorChile']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% "Honduras"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorHonduras'] - 
                                                                   1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% 25] <- round(exp(log(fm_exp_df$RR[fm_exp_df$Variable %in% 25]) - 
                                                           1.96*SE_25), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% 45] <- round(exp(log(fm_exp_df$RR[fm_exp_df$Variable %in% 45]) - 
                                                           1.96*SE_45), digits =3)






# SE with sandwich estimates Upper limits
fm_exp_df$Upper[fm_exp_df$Variable %in% "Pre and during DTG"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG'] + 
                                                                             1.96*fm_summ$SE[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'ade_typeYes'] + 
                                                                                           1.96*fm_summ$SE[fm_summ$Variables %in% 'ade_typeYes']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% "TB: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'tbYes'] + 
                                                                  1.96*fm_summ$SE[fm_summ$Variables %in% 'tbYes']), digits =3)

fm_exp_df$Upper[fm_exp_df$Variable %in% "Brazil"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorBrazil'] + 
                                                                 1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% "Chile"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorChile'] + 
                                                                1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorChile']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% "Honduras"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorHonduras'] + 
                                                                   1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% 25] <- round(exp(log(fm_exp_df$RR[fm_exp_df$Variable %in% 25]) + 
                                                           1.96*SE_25), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% 45] <- round(exp(log(fm_exp_df$RR[fm_exp_df$Variable %in% 45]) + 
                                                           1.96*SE_45), digits =3)

# P value computation with sandwich estimates
z1 <- fm_summ$Coefs[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']/fm_summ$SE[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']
z2 <- fm_summ$Coefs[fm_summ$Variables %in% 'ade_typeYes']/fm_summ$SE[fm_summ$Variables %in% 'ade_typeYes']
z3 <- fm_summ$Coefs[fm_summ$Variables %in% 'tbYes']/fm_summ$SE[fm_summ$Variables %in% 'tbYes']
z5 <- fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorBrazil']/fm_summ$SE[fm_summ$Variables %in% 'study_site.factorBrazil']
z6 <- fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorChile']/fm_summ$SE[fm_summ$Variables %in% 'study_site.factorChile']
z7 <- fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorHonduras']/fm_summ$SE[fm_summ$Variables %in% 'study_site.factorHonduras']
z8 <- log(fm_exp_df$RR[fm_exp_df$Variable %in% 25])/SE_25
z9 <- log(fm_exp_df$RR[fm_exp_df$Variable %in% 45])/SE_45



fm_exp_df$p[fm_exp_df$Variable %in% "Pre and during DTG"] <- round((1 - pnorm(abs(z1))) * 2,digits = 3)
fm_exp_df$p[fm_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round((1 - pnorm(abs(z2))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% "TB: Yes"] <- round((1 - pnorm(abs(z3))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% "Brazil"] <- round((1 - pnorm(abs(z5))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% "Chile"] <- round((1 - pnorm(abs(z6))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% "Honduras"] <- round((1 - pnorm(abs(z7))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% 25] <- round((1 - pnorm(abs(z8))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% 45] <- round((1 - pnorm(abs(z9))) * 2,digits =3)


## Interaction term model with 2 categories for DTG warning
data_logit <- data_m %>% select(age , study_site.factor ,  DTG_warning,ade_type , tb, DTG,DTG_pre.factor)

data_logit <- within(data_logit, DTG_pre.factor <- relevel(DTG_pre.factor, ref = "Post DTG"))

m2 <- glm(DTG ~ ns(age, df = 4) + DTG_pre.factor + study_site.factor  +  ade_type + tb  + DTG_pre.factor*(ns(age, df = 4)) ,data= data_logit,family=poisson(link="log"))
summary(m2)

data_logit2 <- within(data_logit, DTG_pre.factor <- relevel(DTG_pre.factor, ref = "Pre and during DTG"))
m4 <- glm(DTG ~ ns(age, df = 4) + DTG_pre.factor + study_site.factor  +  ade_type + tb  + DTG_pre.factor*(ns(age, df = 4)) ,data= data_logit2,family=poisson(link="log"))
summary(m4)


ages <- round(min(data_logit$age)):round(max(data_logit$age))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

m2_summ <- summary(m2, conf.int=TRUE)
covar <- sandwich(m2)
se <- sqrt(covar[row(covar)==col(covar)])
m2_summ <- as.data.frame(cbind(m2_summ$coefficients[,1], se))
m2_summ$vars <- rownames(m2_summ)
colnames(m2_summ) <- c("Coefs", "SE", "Variables")

m4_summ <- summary(m4, conf.int=TRUE)
covar4 <- sandwich(m4)
se4 <- sqrt(covar4[row(covar4)==col(covar4)])
m4_summ <- as.data.frame(cbind(m4_summ$coefficients[,1], se4))
m4_summ$vars <- rownames(m4_summ)
colnames(m4_summ) <- c("Coefs", "SE", "Variables")



# Creating summary data frame from 2 different models
mm2_exp_df <- data.frame(Variable= c(c("Age (ref = 35 age and Post DTG warning) ", 20,25,30,40,45,50,55,60,65),
                                     c("Site (ref= Haiti)", "Brazil", "Chile", "Honduras"),
                                     c("DTG: post vs. Pre and during DTG warning(ref = Post DTG warning) at median age of 35","Pre and during DTG warning"),
                                     c("Other AIDS defining illness (ref =No)", "Other AIDS defining illness: Yes"),
                                     c("TB (ref = No)", "TB: Yes"),
                                     c("Age (ref = 35 age and Pre-During DTG warning) ", "20 (Pre-During DTG warning)","25 (Pre-During DTG warning)",
                                       "30 (Pre-During DTG warning)","40 (Pre-During DTG warning)","45 (Pre-During DTG warning)","50 (Pre-During DTG warning)",
                                       "55 (Pre-During DTG warning)","60 (Pre-During DTG warning)","65 (Pre-During DTG warning)")),
                         RR= NA, Lower=NA, Upper=NA,  p=NA)


# creating RR and Limits
mm2_exp_df$RR[mm2_exp_df$Variable %in% "Pre and during DTG warning"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'DTG_pre.factorPre and during DTG']  +
                                                                                    m2_summ$Coefs[m2_summ$Variables %in% "ns(age, df = 4)1:DTG_pre.factorPre and during DTG"]*(sp[rownames(sp) %in% 35, 1])+ 
                                                                                    m2_summ$Coefs[m2_summ$Variables %in% "ns(age, df = 4)2:DTG_pre.factorPre and during DTG"]*(sp[rownames(sp) %in% 35, 2])+                                                          
                                                                                    m2_summ$Coefs[m2_summ$Variables %in% "ns(age, df = 4)3:DTG_pre.factorPre and during DTG"]*(sp[rownames(sp) %in% 35, 3])+                                                          
                                                                                    m2_summ$Coefs[m2_summ$Variables %in% "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]*(sp[rownames(sp) %in% 35, 4])), digits =3)

mm2_exp_df$RR[mm2_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'ade_typeYes']), digits =3)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "TB: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'tbYes']), digits =3)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "Brazil"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "Chile"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorChile']), digits =3)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "Honduras"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorHonduras']), digits =3)


# function for relative risk


rr <- function(myage,model){
  exp(model$Coefs[model$Variables %in% "ns(age, df = 4)1"]*
        (sp[rownames(sp) %in% myage, 1]-sp[rownames(sp) %in% 35, 1]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)2"]*
        (sp[rownames(sp) %in% myage, 2]- sp[rownames(sp) %in% 35, 2]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)3"]*
        (sp[rownames(sp) %in% myage, 3]- sp[rownames(sp) %in% 35, 3]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)4"]*
        (sp[rownames(sp) %in% myage, 4]- sp[rownames(sp) %in% 35, 4]))
  
}

mm2_exp_df$RR[mm2_exp_df$Variable %in% 20] <- rr(20,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% 25] <- rr(25,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% 30] <- rr(30,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% 40] <- rr(40,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% 45] <- rr(45,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% 50] <- rr(50,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% 55] <- rr(55,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% 60] <- rr(60,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% 65] <- rr(65,m2_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "25 (Pre-During DTG warning)"] <- rr(25,m4_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "20 (Pre-During DTG warning)"] <- rr(20,m4_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "30 (Pre-During DTG warning)"] <- rr(30,m4_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "40 (Pre-During DTG warning)"] <- rr(40,m4_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "45 (Pre-During DTG warning)"] <- rr(45,m4_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "50 (Pre-During DTG warning)"] <- rr(50,m4_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "55 (Pre-During DTG warning)"] <- rr(55,m4_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "60 (Pre-During DTG warning)"] <- rr(60,m4_summ)
mm2_exp_df$RR[mm2_exp_df$Variable %in% "65 (Pre-During DTG warning)"] <- rr(65,m4_summ)





spline_se_fun <- function(myage,variance){
  
  tmp <- variance
  a <- sp$X1[rownames(sp) %in% myage]-sp$X1[rownames(sp) %in% 35]
  b <- sp$X2[rownames(sp) %in% myage]-sp$X2[rownames(sp) %in% 35]
  c <- sp$X3[rownames(sp) %in% myage]-sp$X3[rownames(sp) %in% 35]
  d <- sp$X4[rownames(sp) %in% myage]-sp$X4[rownames(sp) %in% 35]
  
  Var_x <- tmp["ns(age, df = 4)1", "ns(age, df = 4)1"]
  Var_y <- tmp["ns(age, df = 4)2", "ns(age, df = 4)2"]
  Var_z <- tmp["ns(age, df = 4)3", "ns(age, df = 4)3"]
  Var_w <- tmp["ns(age, df = 4)4", "ns(age, df = 4)4"]
  
  
  Cov_xy <- tmp["ns(age, df = 4)1", "ns(age, df = 4)2"]
  Cov_xz <- tmp["ns(age, df = 4)1", "ns(age, df = 4)3"]
  Cov_xw <- tmp["ns(age, df = 4)1", "ns(age, df = 4)4"]
  Cov_yz <- tmp["ns(age, df = 4)2", "ns(age, df = 4)3"]
  Cov_yw <- tmp["ns(age, df = 4)2", "ns(age, df = 4)4"]
  Cov_zw <- tmp["ns(age, df = 4)3", "ns(age, df = 4)4"]
  
  
  
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w 
  2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw 
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

spline_se_exfun <- function(myage,variance){
  
  tmp <- variance
  a <- sp$X1[rownames(sp) %in% myage]
  b <- sp$X2[rownames(sp) %in% myage]
  c <- sp$X3[rownames(sp) %in% myage]
  d <- sp$X4[rownames(sp) %in% myage]
  
  Var_x <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "ns(age, df = 4)1:DTG_pre.factorPre and during DTG"]
  Var_y <- tmp["ns(age, df = 4)2:DTG_pre.factorPre and during DTG", "ns(age, df = 4)2:DTG_pre.factorPre and during DTG"]
  Var_z <- tmp["ns(age, df = 4)3:DTG_pre.factorPre and during DTG", "ns(age, df = 4)3:DTG_pre.factorPre and during DTG"]
  Var_w <- tmp["ns(age, df = 4)4:DTG_pre.factorPre and during DTG", "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]
  var_v <- tmp["DTG_pre.factorPre and during DTG", "DTG_pre.factorPre and during DTG"]
  
  Cov_xy <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "ns(age, df = 4)2:DTG_pre.factorPre and during DTG"]
  Cov_xz <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "ns(age, df = 4)3:DTG_pre.factorPre and during DTG"]
  Cov_xw <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]
  Cov_yz <- tmp["ns(age, df = 4)2:DTG_pre.factorPre and during DTG", "ns(age, df = 4)3:DTG_pre.factorPre and during DTG"]
  Cov_yw <- tmp["ns(age, df = 4)2:DTG_pre.factorPre and during DTG", "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]
  Cov_zw <- tmp["ns(age, df = 4)3:DTG_pre.factorPre and during DTG", "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]
  Cov_wv <- tmp["ns(age, df = 4)4:DTG_pre.factorPre and during DTG","DTG_pre.factorPre and during DTG"]
  Cov_xv <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "DTG_pre.factorPre and during DTG"]
  Cov_yv <- tmp["ns(age, df = 4)2:DTG_pre.factorPre and during DTG", "DTG_pre.factorPre and during DTG"]
  Cov_zv <- tmp["ns(age, df = 4)3:DTG_pre.factorPre and during DTG", "DTG_pre.factorPre and during DTG"]
  
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w + var_v +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw + 2*a*Cov_xv + 2*b*Cov_yv + 2*c*Cov_zv + 2*d*Cov_wv
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  
  print(pooled_se_splines_yr)
  
}

SE_20 <- spline_se_fun(20,covar)
SE_25 <-spline_se_fun(25,covar)
SE_30 <- spline_se_fun(30,covar)
SE_40 <-spline_se_fun(40,covar)
SE_45 <- spline_se_fun(45,covar)
SE_50 <-spline_se_fun(50,covar)
SE_55 <- spline_se_fun(55,covar)
SE_60 <-spline_se_fun(60,covar)
SE_65 <- spline_se_fun(65,covar)
SE_20_m4 <- spline_se_fun(20,covar4)
SE_25_m4 <-spline_se_fun(25,covar4)
SE_30_m4 <- spline_se_fun(30,covar4)
SE_40_m4 <-spline_se_fun(40,covar4)
SE_45_m4 <- spline_se_fun(45,covar4)
SE_50_m4 <-spline_se_fun(50,covar4)
SE_55_m4 <- spline_se_fun(55,covar4)
SE_60_m4 <-spline_se_fun(60,covar4)
SE_65_m4 <- spline_se_fun(65,covar4)
SE_35_DTG <- spline_se_exfun(35,covar)



# SE with sandwich estimates lower limits
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "Pre and during DTG warning"] <- round(exp(log(mm2_exp_df$RR[m2_summ$Variables %in% 'DTG_pre.factorPre and during DTG']) - 
                                                                                       1.96*SE_35_DTG), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'ade_typeYes'] - 
                                                                                             1.96*m2_summ$SE[m2_summ$Variables %in% 'ade_typeYes']), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "TB: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'tbYes'] - 
                                                                    1.96*m2_summ$SE[m2_summ$Variables %in% 'tbYes']), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "Brazil"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorBrazil'] - 
                                                                   1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "Chile"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorChile'] - 
                                                                  1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorChile']), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "Honduras"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorHonduras'] - 
                                                                     1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 20] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 20]) - 
                                                             1.96*SE_20), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 25] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 25]) - 
                                                             1.96*SE_25), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 30] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 30]) - 
                                                             1.96*SE_30), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 40] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 40]) - 
                                                             1.96*SE_40), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 45] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 45]) - 
                                                             1.96*SE_45), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 50] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 50]) - 
                                                             1.96*SE_50), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 55] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 55]) - 
                                                             1.96*SE_55), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 60] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 60]) - 
                                                             1.96*SE_60), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% 65] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 65]) - 
                                                             1.96*SE_65), digits =3)

mm2_exp_df$Lower[mm2_exp_df$Variable %in% "20 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "20 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_20_m4), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "25 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "25 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_25_m4), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "30 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "30 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_30_m4), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "40 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "40 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_40_m4), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "45 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "45 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_45_m4), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "50 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "50 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_50_m4), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "55 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "55 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_55_m4), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "60 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "60 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_60_m4), digits =3)
mm2_exp_df$Lower[mm2_exp_df$Variable %in% "65 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "65 (Pre-During DTG warning)"]) - 
                                                                                        1.96*SE_65_m4), digits =3)





# SE with sandwich estimates Upper limits
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "Pre and during DTG warning"] <- round(exp(log(mm2_exp_df$RR[m2_summ$Variables %in% 'DTG_pre.factorPre and during DTG']) + 
                                                                                       1.96*SE_35_DTG), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'ade_typeYes'] + 
                                                                                             1.96*m2_summ$SE[m2_summ$Variables %in% 'ade_typeYes']), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "TB: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'tbYes'] + 
                                                                    1.96*m2_summ$SE[m2_summ$Variables %in% 'tbYes']), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "Brazil"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorBrazil'] + 
                                                                   1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "Chile"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorChile'] + 
                                                                  1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorChile']), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "Honduras"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorHonduras'] + 
                                                                     1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 20] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 20])  + 
                                                             1.96*SE_20), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 25] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 25])  + 
                                                             1.96*SE_25), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 30] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 30])  + 
                                                             1.96*SE_30), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 40] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 40])  + 
                                                             1.96*SE_40), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 45] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 45])  + 
                                                             1.96*SE_45), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 50] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 50])  + 
                                                             1.96*SE_50), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 55] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 55])  + 
                                                             1.96*SE_55), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 60] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 60])  + 
                                                             1.96*SE_60), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% 65] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 65])  + 
                                                             1.96*SE_65), digits =3)

mm2_exp_df$Upper[mm2_exp_df$Variable %in% "20 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "20 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_20_m4), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "25 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "25 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_25_m4), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "30 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "30 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_30_m4), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "40 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "40 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_40_m4), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "45 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "45 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_45_m4), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "50 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "50 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_50_m4), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "55 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "55 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_55_m4), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "60 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "60 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_60_m4), digits =3)
mm2_exp_df$Upper[mm2_exp_df$Variable %in% "65 (Pre-During DTG warning)"] <- round(exp(log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "65 (Pre-During DTG warning)"])  + 
                                                                                        1.96*SE_65_m4), digits =3)



# P value computation with sandwich estimates
z1 <- log(mm2_exp_df$RR[m2_summ$Variables %in% 'DTG_pre.factorPre and during DTG'])/SE_35_DTG
z2 <- m2_summ$Coefs[m2_summ$Variables %in% 'ade_typeYes']/m2_summ$SE[m2_summ$Variables %in% 'ade_typeYes']
z3 <- m2_summ$Coefs[m2_summ$Variables %in% 'tbYes']/m2_summ$SE[m2_summ$Variables %in% 'tbYes']
z5 <- m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorBrazil']/m2_summ$SE[m2_summ$Variables %in% 'study_site.factorBrazil']
z6 <- m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorChile']/m2_summ$SE[m2_summ$Variables %in% 'study_site.factorChile']
z7 <- m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorHonduras']/m2_summ$SE[m2_summ$Variables %in% 'study_site.factorHonduras']
z8 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 20])/SE_25
z9 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 25])/SE_25
z10 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 30])/SE_30
z11 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 40])/SE_40
z12 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 45])/SE_45
z13 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 50])/SE_50
z14 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 55])/SE_55
z15 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 60])/SE_60
z16 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% 65])/SE_65
z17 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "20 (Pre-During DTG warning)"])/SE_20_m4
z18 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "25 (Pre-During DTG warning)"])/SE_25_m4
z19 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "30 (Pre-During DTG warning)"])/SE_30_m4
z20 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "40 (Pre-During DTG warning)"])/SE_40_m4
z21 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "45 (Pre-During DTG warning)"])/SE_45_m4
z22 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "50 (Pre-During DTG warning)"])/SE_50_m4
z23 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "55 (Pre-During DTG warning)"])/SE_55_m4
z24 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "60 (Pre-During DTG warning)"])/SE_60_m4
z25 <- log(mm2_exp_df$RR[mm2_exp_df$Variable %in% "65 (Pre-During DTG warning)"])/SE_65_m4



mm2_exp_df$p[mm2_exp_df$Variable %in% "Pre and during DTG warning"] <- round((1 - pnorm(abs(z1))) * 2,digits = 3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round((1 - pnorm(abs(z2))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "TB: Yes"] <- round((1 - pnorm(abs(z3))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "Brazil"] <- round((1 - pnorm(abs(z5))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "Chile"] <- round((1 - pnorm(abs(z6))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "Honduras"] <- round((1 - pnorm(abs(z7))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 20] <- round((1 - pnorm(abs(z8))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 25] <- round((1 - pnorm(abs(z9))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 30] <- round((1 - pnorm(abs(z10))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 40] <- round((1 - pnorm(abs(z11))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 45] <- round((1 - pnorm(abs(z12))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 50] <- round((1 - pnorm(abs(z13))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 55] <- round((1 - pnorm(abs(z14))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 60] <- round((1 - pnorm(abs(z15))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% 65] <- round((1 - pnorm(abs(z16))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "20 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z17))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "25 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z18))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "30 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z19))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "40 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z20))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "45 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z21))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "50 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z22))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "55 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z23))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "60 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z24))) * 2,digits =3)
mm2_exp_df$p[mm2_exp_df$Variable %in% "65 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z25))) * 2,digits =3)




save(fm_exp_df,file ="fm_exp_df.Rdata")
save(unmm,file ="unmm.Rdata")
save(int_m,file ="int_m.Rdata")
save(int_m_cat,file ="int_m_cat.Rdata")
save(overall_mm_dummy,file ="overall_mm_dummy.Rdata")
save(overall_mm_cat,file ="overall_mm_cat.Rdata")
save(mm2_exp_df,file ="mm2_exp_df.Rdata")
save(dem_male,file ="dem_male.Rdata")

## Females only analysis 

data_f <- demdata_1 %>% filter(gender.factor == "Female")
dem_female <- demdata_1 %>% filter(gender.factor == "Female")




data_f$age_period.factor <- with(data_f,ifelse(age_gt50.factor == "Greater then or equal to 50" & DTG_pre.factor == "Pre and during DTG", "older pre",
                                               ifelse(age_gt50.factor == "Less then 50"                & DTG_pre.factor == "Pre and during DTG", "younger pre",
                                                      ifelse(age_gt50.factor == "Greater then or equal to 50" & DTG_pre.factor == "Post DTG"           , "older post",
                                                             ifelse(age_gt50.factor == "Less then 50" & DTG_pre.factor == "Post DTG"  ,  "younger post",
                                                                    "99")))))     

#data_f$DTG.factor <- ifelse(data_f$DTG.factor == "No", 0, 1)

dd <- datadist(data_f)
options(datadist="dd")
dd$limits$age_period.factor <- "younger post"
dd$limits$tb <- "No"
dd$limits$ade_type <- "No"

# Setting levels  
data_f <- within(data_f, age_period.factor <- relevel(factor(age_period.factor), ref = "younger post"))
data_f <- within(data_f, tb <- relevel(factor(tb), ref = "No"))
data_f <- within(data_f, ade_type <- relevel(factor(ade_type), ref = "No"))




# Creating data frame
unfm <- data.frame(Variable= c(c("Age (ref =Less then 50)", "Greater then or equal to 50"),
                               c("Site (ref= Haiti)", "Brazil", "Chile", "Honduras"),
                               c("DTG: post vs. Pre and during DTG (ref = Post DTG)","Pre and during DTG"),
                               c("Other AIDS defining illness (ref =No)", "Other AIDS defining illness: Yes"),
                               c("TB (ref = No)", "TB: Yes")),
                   RR= NA, Lower=NA, Upper=NA,  p=NA)

# un-adjusted models
#age_gt50
m1_un_age_gt50 <- glm(DTG ~   age_gt50.factor, data = data_f, family=poisson(link="log"))
m1_un_Age_gt50 <- Glm(DTG ~   age_gt50.factor, data = data_f, family=poisson(link="log"))
summ_un_age_gt50 <- as.data.frame(my_summary.rms(object=m1_un_Age_gt50, object2 = m1_un_age_gt50))
exposures <- " "
anova_un_age_gt50 <- as.data.frame(anova_w_sandwich(object = m1_un_Age_gt50, objectglm = m1_un_age_gt50,coVars_woInt = c("age_gt50")))

#TB
m1_un_tb <- glm(DTG ~   tb, data = data_f, family=poisson(link="log"))
m1_un_Tb <- Glm(DTG ~   tb, data = data_f, family=poisson(link="log"))
summ_un_tb <- as.data.frame(my_summary.rms(object=m1_un_Tb, object2 = m1_un_tb))
exposures <- " "
anova_un_tb <- as.data.frame(anova_w_sandwich(object = m1_un_Tb, objectglm = m1_un_tb,coVars_woInt = c("tb")))

#ADE
m1_un_ade <- glm(DTG ~   ade_type, data = data_f, family=poisson(link="log"))
m1_un_Ade <- Glm(DTG ~   ade_type, data = data_f, family=poisson(link="log"))
summ_un_ade <- as.data.frame(my_summary.rms(object=m1_un_Ade, object2 = m1_un_ade))
exposures <- " "
anova_un_ade <- as.data.frame(anova_w_sandwich(object = m1_un_Ade, objectglm = m1_un_ade,coVars_woInt = c("ade_type")))


#site
m1_un_site <- glm(DTG ~   study_site.factor, data = data_f, family=poisson(link="log"))
m1_un_Site <- Glm(DTG ~   study_site.factor, data = data_f, family=poisson(link="log"))
summ_un_site <- as.data.frame(my_summary.rms(object=m1_un_Site, object2 = m1_un_site))
exposures <- " "
anova_un_site <- as.data.frame(anova_w_sandwich(object = m1_un_Site, objectglm = m1_un_site,coVars_woInt = c("study_site.factor")))


#DTG
m1_un_dtg <- glm(DTG ~   DTG_pre.factor, data = data_f, family=poisson(link="log"))
m1_un_Dtg <- Glm(DTG ~   DTG_pre.factor, data = data_f, family=poisson(link="log"))
summ_un_dtg <- as.data.frame(my_summary.rms(object=m1_un_Dtg, object2 = m1_un_dtg))
exposures <- " "
anova_un_dtg <- as.data.frame(anova_w_sandwich(object = m1_un_Dtg, objectglm = m1_un_dtg,coVars_woInt = c("DTG_pre")))


##Unadjusted Model
unfm$RR[unfm$Variable %in% "Pre and during DTG"] <- round(exp(summ_un_dtg$Effect[1]), digits =3)
unfm$RR[unfm$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(summ_un_ade$Effect[1]), digits =3)
unfm$RR[unfm$Variable %in% "TB: Yes"] <- round(exp(summ_un_tb$Effect[1]), digits =3)
unfm$RR[unfm$Variable %in% "Brazil"] <- round(exp(summ_un_site$Effect[1]), digits =3)
unfm$RR[unfm$Variable %in% "Chile"] <- round(exp(summ_un_site$Effect[2]), digits =3)
unfm$RR[unfm$Variable %in% "Honduras"] <- round(exp(summ_un_site$Effect[3]), digits =3)
unfm$RR[unfm$Variable %in% "Greater then or equal to 50"] <- round(exp(summ_un_age_gt50$Effect[1]), digits =3)
unfm$Lower[unfm$Variable %in% "Pre and during DTG"] <- round(exp(summ_un_dtg$`Lower 0.95`[1]), digits =3)
unfm$Lower[unfm$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(summ_un_ade$`Lower 0.95`[1]), digits =3)
unfm$Lower[unfm$Variable %in% "TB: Yes"] <- round(exp(summ_un_tb$`Lower 0.95`[1]), digits =3)
unfm$Lower[unfm$Variable %in% "Brazil"] <- round(exp(summ_un_site$`Lower 0.95`[1]), digits =3)
unfm$Lower[unfm$Variable %in% "Chile"] <- round(exp(summ_un_site$`Lower 0.95`[2]), digits =3)
unfm$Lower[unfm$Variable %in% "Honduras"] <- round(exp(summ_un_site$`Lower 0.95`[3]), digits =3)
unfm$Lower[unfm$Variable %in% "Greater then or equal to 50"] <- round(exp(summ_un_age_gt50$`Lower 0.95`[1]), digits =3)
unfm$Upper[unfm$Variable %in% "Pre and during DTG"] <- round(exp(summ_un_dtg$`Upper 0.95`[1]), digits =3)
unfm$Upper[unfm$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(summ_un_ade$`Upper 0.95`[1]), digits =3)
unfm$Upper[unfm$Variable %in% "TB: Yes"] <- round(exp(summ_un_tb$`Upper 0.95`[1]), digits =3)
unfm$Upper[unfm$Variable %in% "Brazil"] <- round(exp(summ_un_site$`Upper 0.95`[1]), digits =3)
unfm$Upper[unfm$Variable %in% "Chile"] <- round(exp(summ_un_site$`Upper 0.95`[2]), digits =3)
unfm$Upper[unfm$Variable %in% "Honduras"] <- round(exp(summ_un_site$`Upper 0.95`[3]), digits =3)
unfm$Upper[unfm$Variable %in% "Greater then or equal to 50"] <- round(exp(summ_un_age_gt50$`Upper 0.95`[1]), digits =3)

# p values for unadjusted table
unfm$p[unfm$Variable %in% "DTG: post vs. Pre and during DTG (ref = Post DTG)"] <- anova_un_dtg$stats[1]
unfm$p[unfm$Variable %in% "Other AIDS defining illness (ref =No)"] <- anova_un_ade$stats[1]
unfm$p[unfm$Variable %in% "TB (ref = No)"] <- anova_un_tb$stats[1]
unfm$p[unfm$Variable %in% "Site (ref= Haiti)"] <- anova_un_site$stats[1]
unfm$p[unfm$Variable %in% "Age (ref = Less then 50)"] <- anova_un_age_gt50$stats[1]


#Interaction model dummy variables

mod_dummy_interaction <- glm(DTG ~   age_period.factor, data = data_f, family=poisson(link="log"))
mod_dummy_Interaction <- Glm(DTG ~   age_period.factor, data = data_f, family=poisson(link="log"))
exposures <- " "
summ_int <- as.data.frame(my_summary.rms(object=mod_dummy_Interaction, object2 = mod_dummy_interaction))
var <- sandwich(object2)[-1, -1] # Removes Intercept row/column
object2 = mod_dummy_interaction
z1 <- summ_int$Effect[1]/summ_int$S.E.[1]
z2 <- summ_int$Effect[2]/summ_int$S.E.[2]
z3 <- summ_int$Effect[3]/summ_int$S.E.[3]

summ_int$p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3))
summ_int$Effect <- round(exp(summ_int$Effect), digits =3)
summ_int$`Lower 0.95` <- round(exp(summ_int$`Lower 0.95`), digits =3)
summ_int$`Upper 0.95` <- round(exp(summ_int$`Upper 0.95`), digits =3)
variables <- c("older post:younger post","older pre:younger post","younger pre:younger post")
int_f <- as.data.frame(cbind(variables,summ_int$Effect,summ_int$`Lower 0.95`,summ_int$`Upper 0.95`,summ_int$p))
colnames(int_f) <- c("variables","Effect","Lower","Upper","p-value")



## Interaction terms model
mod_interaction <- glm(DTG ~    age_gt50.factor*DTG_pre.factor , data= data_f, family=poisson(link="log"))
summ_int <- summary(mod_interaction)
var <- sandwich(mod_interaction)

#OP
ind_test = 2
test_beta = summ_int$coefficients[ind_test]
OP_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
OP_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
OP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(var[ind_test,ind_test])


#OPre
ind_test = c(2,3,4)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OPre_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPre_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPre_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#YPre
ind_test = 3
test_beta = summ_int$coefficients[ind_test]
YPre_UB_test = exp(test_beta + 1.96*sqrt(var[ind_test,ind_test]))
YPre_LB_test = exp(test_beta - 1.96*sqrt(var[ind_test,ind_test]))
YPre_point_test = exp(test_beta)
z3 <- test_beta/sqrt(var[ind_test,ind_test])


#older pre vs Older Post
#OPre
ind_test = c(3,4)
test_var = var[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(summ_int$coefficients[ind_test])
OPOST_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPOST_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPOST_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

## Interaction term
variable <- variables <- c("older post:younger post","older pre:younger post","younger pre:younger post","older pre : older Post")
rr <- c(OP_point_test,OPre_point_test,YPre_point_test,OPOST_point_test)
ll <- c(OP_LB_test,OPre_LB_test,YPre_LB_test,OPOST_LB_test)
ul <- c(OP_UB_test,OPre_UB_test,YPre_UB_test,OPOST_UB_test)
p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3))


int_f_cat <-  data.frame(variable, rr,ll,ul,p)


#Dummy Overall term
mod_dummy_a <- glm(DTG ~    study_site.factor + tb + ade_type + age_period.factor, data= data_f, family=poisson(link="log"))
mod_Dummy_a <- Glm(DTG ~  study_site.factor  + tb + ade_type + age_period.factor , data= data_f, family=poisson(link="log"))
summ_dummy_a <- as.data.frame(my_summary.rms(object = mod_Dummy_a, object2 = mod_dummy_a))
var_dummy_a <- sandwich(mod_dummy_a)

# Table for adjusted model with dummy variables
variable <- c("Brazil","Chile","Honduras","Tb (ref =no)","Other AIDS defining illness (ref =no)","older post:younger post", "older pre:younger post","younger pre:younger post")
rr <- c(round(exp(summ_dummy_a$Effect[1]), digits =3),round(exp(summ_dummy_a$Effect[2]), digits =3),round(exp(summ_dummy_a$Effect[3]), digits =3)
        ,round(exp(summ_dummy_a$Effect[4]), digits =3),round(exp(summ_dummy_a$Effect[5]), digits =3),round(exp(summ_dummy_a$Effect[6]), digits =3)
        ,round(exp(summ_dummy_a$Effect[7]), digits =3),round(exp(summ_dummy_a$Effect[8]), digits =3))

ll <- c(round(exp(summ_dummy_a$`Lower 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Lower 0.95`[8]), digits =3))

ul <- c(round(exp(summ_dummy_a$`Upper 0.95`[1]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[2]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[3]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[4]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[5]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[6]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[7]), digits =3),
        round(exp(summ_dummy_a$`Upper 0.95`[8]), digits =3))                                                                                                                                                                                   

z1 <- summ_dummy_a$Effect[1]/summ_dummy_a$S.E.[1]
z2 <- summ_dummy_a$Effect[2]/summ_dummy_a$S.E.[2]
z3 <- summ_dummy_a$Effect[3]/summ_dummy_a$S.E.[3]
z4 <- summ_dummy_a$Effect[4]/summ_dummy_a$S.E.[4]
z5 <- summ_dummy_a$Effect[5]/summ_dummy_a$S.E.[5]
z6 <- summ_dummy_a$Effect[6]/summ_dummy_a$S.E.[6]
z7 <- summ_dummy_a$Effect[7]/summ_dummy_a$S.E.[7]
z8 <- summ_dummy_a$Effect[8]/summ_dummy_a$S.E.[8]



p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),
       round((1 - pnorm(abs(z8))) * 2,digits = 3))
overall_fm_dummy <- data.frame(variable,rr,ll,ul,p)


## Overall terms model
oafm <- glm(DTG ~  age_gt50.factor + tb + ade_type + study_site.factor + DTG_pre.factor +  age_gt50.factor*DTG_pre.factor , data= data_f, family=poisson(link="log"))
summ_overall_a <- summary(oafm)
var_overall_a <- sandwich(oafm)


ind_atest = 3
atest_beta = summ_overall_a$coefficients[ind_atest]
tb_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
tb_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
tb_point_atest = exp(atest_beta)
aztb <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])

ind_atest = 4
atest_beta = summ_overall_a$coefficients[ind_atest]
ade_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
ade_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
ade_point_atest = exp(atest_beta)
azade <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])




ind_atest = 5
atest_beta = summ_overall_a$coefficients[ind_atest]
bra_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
bra_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
bra_point_atest = exp(atest_beta)
azbra <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])

ind_atest = 6
atest_beta = summ_overall_a$coefficients[ind_atest]
chile_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
chile_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
chile_point_atest = exp(atest_beta)
azchile <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])

ind_atest = 7
atest_beta = summ_overall_a$coefficients[ind_atest]
hon_UB_atest = exp(atest_beta + 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
hon_LB_atest = exp(atest_beta - 1.96*sqrt(var_overall_a[ind_atest,ind_atest]))
hon_point_atest = exp(atest_beta)
azhon <- atest_beta/sqrt(var_overall_a[ind_atest,ind_atest])

#OP
ind_test = 2
test_beta = summ_overall_a$coefficients[ind_test]
OP_UB_test = exp(test_beta + 1.96*sqrt(var_overall_a[ind_test,ind_test]))
OP_LB_test = exp(test_beta - 1.96*sqrt(var_overall_a[ind_test,ind_test]))
OP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(var_overall_a[ind_test,ind_test])


#OPre
ind_test = c(2,8,9)
test_var_overall_a = var_overall_a[ind_test, ind_test]
se_test = sum(diag(test_var_overall_a)) + 2*sum(test_var_overall_a[upper.tri(test_var_overall_a, diag = FALSE)])
test_beta = sum(summ_overall_a$coefficients[ind_test])
OPre_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPre_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPre_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#YPre
ind_test = 8
test_beta = summ_overall_a$coefficients[ind_test]
YPre_UB_test = exp(test_beta + 1.96*sqrt(var_overall_a[ind_test,ind_test]))
YPre_LB_test = exp(test_beta - 1.96*sqrt(var_overall_a[ind_test,ind_test]))
YPre_point_test = exp(test_beta)
z3 <- test_beta/sqrt(var_overall_a[ind_test,ind_test])


#older pre vs Older Post
#OPre
ind_test = c(8,9)
test_var_overall_a = var_overall_a[ind_test, ind_test]
se_test = sum(diag(test_var_overall_a)) + 2*sum(test_var_overall_a[upper.tri(test_var_overall_a, diag = FALSE)])
test_beta = sum(summ_overall_a$coefficients[ind_test])
OPOST_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPOST_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPOST_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)


#Younger Pre DTG vs. Older Post
ind_test = c(8,1)
test_var_overall_a = var_overall_a[ind_test, ind_test]
se_test = sum(diag(test_var_overall_a)) - 2*sum(test_var_overall_a[upper.tri(test_var_overall_a, diag = FALSE)])
test_beta = (summ_overall_a$coefficients[9] - summ_overall_a$coefficients[1])
YPRE_OP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YPRE_OP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YPRE_OP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(se_test)



#Older Pre DTG vs. Younger Pre
ind_test = c(2,9)
test_var_overall_a = var_overall_a[ind_test, ind_test]
se_test = sum(diag(test_var_overall_a)) + 2*sum(test_var_overall_a[upper.tri(test_var_overall_a, diag = FALSE)])
test_beta = sum(summ_overall_a$coefficients[ind_test])
OPRE_YPRE_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OPRE_YPRE_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OPRE_YPRE_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)



var_overall_aiable <- c("Brazil","Chile","Honduras","Other illness(ref=no)","TB (ref=no)","older post:younger post","older pre:younger post","younger pre:younger post","older pre : older Post", "Younger Pre: Older Post", "Older Pre : Younger Pre")
rr <- c(bra_point_atest,chile_point_atest,hon_point_atest,ade_point_atest,tb_point_atest,OP_point_test,OPre_point_test,YPre_point_test,OPOST_point_test,YPRE_OP_point_test, OPRE_YPRE_point_test)
ll <- c(bra_LB_atest,chile_LB_atest,hon_LB_atest,ade_LB_atest,tb_LB_atest,OP_LB_test,OPre_LB_test,YPre_LB_test,OPOST_LB_test,YPRE_OP_LB_test, OPRE_YPRE_LB_test)
ul <- c(bra_UB_atest,chile_UB_atest,hon_UB_atest,ade_UB_atest,tb_UB_atest,OP_UB_test,OPre_UB_test,YPre_UB_test,OPOST_UB_test,YPRE_OP_UB_test, OPRE_YPRE_UB_test)

p <- c(round((1 - pnorm(abs(azbra))) * 2,digits = 3),
       round((1 - pnorm(abs(azchile))) * 2,digits = 3),round((1 - pnorm(abs(azhon))) * 2,digits = 3),
       round((1 - pnorm(abs(azade))) * 2,digits = 3),round((1 - pnorm(abs(aztb))) * 2,digits = 3),
       round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits =3),round((1 - pnorm(abs(z6))) * 2,digits =3))

overall_fm_cat <-  data.frame(var_overall_aiable, rr,ll,ul,p)





#model
fm <- glm(DTG ~ ns(age, df = 4) + study_site.factor +  DTG_pre.factor +  ade_type + tb,data=data_f,family=poisson(link="log"))


# Spline matrix
ages <- round(min(data_f$age)):round(max(data_f$age))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")


# Summary, coefficients and standard errors from sandwich estimates
fm_summ <- summary(fm, conf.int=TRUE)
covar <- sandwich(fm)
se <- sqrt(covar[row(covar)==col(covar)])
fm_summ <- as.data.frame(cbind(fm_summ$coefficients[,1], se))
fm_summ$vars <- rownames(fm_summ)
colnames(fm_summ) <- c("Coefs", "SE", "Variables")

# Creating datat frame
fm_exp_df <- data.frame(Variable= c(c("Age (ref =35)", 25,45),
                                    c("Site (ref= Haiti)", "Brazil", "Chile", "Honduras"),
                                    c("DTG: post vs. Pre and during DTG (ref = Post DTG)","Pre and during DTG"),
                                    c("Other AIDS defining illness (ref =No)", "Other AIDS defining illness: Yes"),
                                    c("TB (ref = No)", "TB: Yes")),
                        RR= NA, Lower=NA, Upper=NA,  p=NA)



# creating RR and Limits
fm_exp_df$RR[fm_exp_df$Variable %in% "Pre and during DTG"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'ade_typeYes']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "TB: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'tbYes']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "Brazil"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "Chile"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorChile']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% "Honduras"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
fm_exp_df$RR[fm_exp_df$Variable %in% 25] <- exp(fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)1"]*
                                                  (sp[rownames(sp) %in% 25, 1]- sp[rownames(sp) %in% 35, 1]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)2"]*
                                                  (sp[rownames(sp) %in% 25, 2]- sp[rownames(sp) %in% 35, 2]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)3"]*
                                                  (sp[rownames(sp) %in% 25, 3]- sp[rownames(sp) %in% 35, 3]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)4"]*
                                                  (sp[rownames(sp) %in% 25, 4]-sp[rownames(sp) %in% 35, 4]))

fm_exp_df$RR[fm_exp_df$Variable %in% 45] <- exp(fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)1"]*
                                                  (sp[rownames(sp) %in% 45, 1]-sp[rownames(sp) %in% 35, 1]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)2"]*
                                                  (sp[rownames(sp) %in% 45, 2]- sp[rownames(sp) %in% 35, 2]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)3"]*
                                                  (sp[rownames(sp) %in% 45, 3]- sp[rownames(sp) %in% 35, 3]) +
                                                  fm_summ$Coefs[fm_summ$Variables %in% "ns(age, df = 4)4"]*
                                                  (sp[rownames(sp) %in% 45, 4]- sp[rownames(sp) %in% 35, 4]))

# my function to calculate variance and SE
spline_se_fun <- function(myage){
  
  tmp <- covar
  a <- sp$X1[rownames(sp) %in% myage]-sp$X1[rownames(sp) %in% 35]
  b <- sp$X2[rownames(sp) %in% myage]-sp$X2[rownames(sp) %in% 35]
  c <- sp$X3[rownames(sp) %in% myage]-sp$X3[rownames(sp) %in% 35]
  d <- sp$X4[rownames(sp) %in% myage]-sp$X4[rownames(sp) %in% 35]
  
  Var_x <- tmp["ns(age, df = 4)1", "ns(age, df = 4)1"]
  Var_y <- tmp["ns(age, df = 4)2", "ns(age, df = 4)2"]
  Var_z <- tmp["ns(age, df = 4)3", "ns(age, df = 4)3"]
  Var_w <- tmp["ns(age, df = 4)4", "ns(age, df = 4)4"]
  
  Cov_xy <- tmp["ns(age, df = 4)1", "ns(age, df = 4)2"]
  Cov_xz <- tmp["ns(age, df = 4)1", "ns(age, df = 4)3"]
  Cov_xw <- tmp["ns(age, df = 4)1", "ns(age, df = 4)4"]
  Cov_yz <- tmp["ns(age, df = 4)2", "ns(age, df = 4)3"]
  Cov_yw <- tmp["ns(age, df = 4)2", "ns(age, df = 4)4"]
  Cov_zw <- tmp["ns(age, df = 4)3", "ns(age, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

SE_25 <- spline_se_fun(25)
SE_45 <-spline_se_fun(45)
# SE with sandwich estimates lower limits
fm_exp_df$Lower[fm_exp_df$Variable %in% "Pre and during DTG"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG'] - 
                                                                             1.96*fm_summ$SE[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'ade_typeYes'] - 
                                                                                           1.96*fm_summ$SE[fm_summ$Variables %in% 'ade_typeYes']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% "TB: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'tbYes'] - 
                                                                  1.96*fm_summ$SE[fm_summ$Variables %in% 'tbYes']), digits =3)

fm_exp_df$Lower[fm_exp_df$Variable %in% "Brazil"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorBrazil'] - 
                                                                 1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% "Chile"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorChile'] - 
                                                                1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorChile']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% "Honduras"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorHonduras'] - 
                                                                   1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% 25] <- round(exp(log(fm_exp_df$RR[fm_exp_df$Variable %in% 25]) - 
                                                           1.96*SE_25), digits =3)
fm_exp_df$Lower[fm_exp_df$Variable %in% 45] <- round(exp(log(fm_exp_df$RR[fm_exp_df$Variable %in% 45]) - 
                                                           1.96*SE_45), digits =3)






# SE with sandwich estimates Upper limits
fm_exp_df$Upper[fm_exp_df$Variable %in% "Pre and during DTG"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG'] + 
                                                                             1.96*fm_summ$SE[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'ade_typeYes'] + 
                                                                                           1.96*fm_summ$SE[fm_summ$Variables %in% 'ade_typeYes']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% "TB: Yes"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'tbYes'] + 
                                                                  1.96*fm_summ$SE[fm_summ$Variables %in% 'tbYes']), digits =3)

fm_exp_df$Upper[fm_exp_df$Variable %in% "Brazil"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorBrazil'] + 
                                                                 1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% "Chile"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorChile'] + 
                                                                1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorChile']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% "Honduras"] <- round(exp(fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorHonduras'] + 
                                                                   1.96*fm_summ$SE[fm_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% 25] <- round(exp(log(fm_exp_df$RR[fm_exp_df$Variable %in% 25]) + 
                                                           1.96*SE_25), digits =3)
fm_exp_df$Upper[fm_exp_df$Variable %in% 45] <- round(exp(log(fm_exp_df$RR[fm_exp_df$Variable %in% 45]) + 
                                                           1.96*SE_45), digits =3)

# P value computation with sandwich estimates
z1 <- fm_summ$Coefs[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']/fm_summ$SE[fm_summ$Variables %in% 'DTG_pre.factorPre and during DTG']
z2 <- fm_summ$Coefs[fm_summ$Variables %in% 'ade_typeYes']/fm_summ$SE[fm_summ$Variables %in% 'ade_typeYes']
z3 <- fm_summ$Coefs[fm_summ$Variables %in% 'tbYes']/fm_summ$SE[fm_summ$Variables %in% 'tbYes']
z5 <- fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorBrazil']/fm_summ$SE[fm_summ$Variables %in% 'study_site.factorBrazil']
z6 <- fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorChile']/fm_summ$SE[fm_summ$Variables %in% 'study_site.factorChile']
z7 <- fm_summ$Coefs[fm_summ$Variables %in% 'study_site.factorHonduras']/fm_summ$SE[fm_summ$Variables %in% 'study_site.factorHonduras']
z8 <- log(fm_exp_df$RR[fm_exp_df$Variable %in% 25])/SE_25
z9 <- log(fm_exp_df$RR[fm_exp_df$Variable %in% 45])/SE_45



fm_exp_df$p[fm_exp_df$Variable %in% "Pre and during DTG"] <- round((1 - pnorm(abs(z1))) * 2,digits = 3)
fm_exp_df$p[fm_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round((1 - pnorm(abs(z2))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% "TB: Yes"] <- round((1 - pnorm(abs(z3))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% "Brazil"] <- round((1 - pnorm(abs(z5))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% "Chile"] <- round((1 - pnorm(abs(z6))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% "Honduras"] <- round((1 - pnorm(abs(z7))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% 25] <- round((1 - pnorm(abs(z8))) * 2,digits =3)
fm_exp_df$p[fm_exp_df$Variable %in% 45] <- round((1 - pnorm(abs(z9))) * 2,digits =3)


## Interaction term model with 2 categories for DTG warning
data_logit <- data_f %>% select(age , study_site.factor ,  DTG_warning,ade_type , tb, DTG,DTG_pre.factor)

data_logit <- within(data_logit, DTG_pre.factor <- relevel(DTG_pre.factor, ref = "Post DTG"))

m2 <- glm(DTG ~ ns(age, df = 4) + DTG_pre.factor + study_site.factor  +  ade_type + tb  + DTG_pre.factor*(ns(age, df = 4)) ,data= data_logit,family=poisson(link="log"))
summary(m2)

data_logit2 <- within(data_logit, DTG_pre.factor <- relevel(DTG_pre.factor, ref = "Pre and during DTG"))
m4 <- glm(DTG ~ ns(age, df = 4) + DTG_pre.factor + study_site.factor  +  ade_type + tb  + DTG_pre.factor*(ns(age, df = 4)) ,data= data_logit2,family=poisson(link="log"))
summary(m4)


ages <- round(min(data_logit$age)):round(max(data_logit$age))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

m2_summ <- summary(m2, conf.int=TRUE)
covar <- sandwich(m2)
se <- sqrt(covar[row(covar)==col(covar)])
m2_summ <- as.data.frame(cbind(m2_summ$coefficients[,1], se))
m2_summ$vars <- rownames(m2_summ)
colnames(m2_summ) <- c("Coefs", "SE", "Variables")

m4_summ <- summary(m4, conf.int=TRUE)
covar4 <- sandwich(m4)
se4 <- sqrt(covar4[row(covar4)==col(covar4)])
m4_summ <- as.data.frame(cbind(m4_summ$coefficients[,1], se4))
m4_summ$vars <- rownames(m4_summ)
colnames(m4_summ) <- c("Coefs", "SE", "Variables")



# Creating summary data frame from 2 different models
m2_exp_df <- data.frame(Variable= c(c("Age (ref = 35 age and Post DTG warning) ", 20,25,30,40,45,50,55,60,65),
                                    c("Site (ref= Haiti)", "Brazil", "Chile", "Honduras"),
                                    c("DTG: post vs. Pre and during DTG warning(ref = Post DTG warning) at median age of 35","Pre and during DTG warning"),
                                    c("Other AIDS defining illness (ref =No)", "Other AIDS defining illness: Yes"),
                                    c("TB (ref = No)", "TB: Yes"),
                                    c("Age (ref = 35 age and Pre-During DTG warning) ", "20 (Pre-During DTG warning)","25 (Pre-During DTG warning)",
                                      "30 (Pre-During DTG warning)","40 (Pre-During DTG warning)","45 (Pre-During DTG warning)","50 (Pre-During DTG warning)",
                                      "55 (Pre-During DTG warning)","60 (Pre-During DTG warning)","65 (Pre-During DTG warning)")),
                        RR= NA, Lower=NA, Upper=NA,  p=NA)


# creating RR and Limits
m2_exp_df$RR[m2_exp_df$Variable %in% "Pre and during DTG warning"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'DTG_pre.factorPre and during DTG']  +
                                                                                  m2_summ$Coefs[m2_summ$Variables %in% "ns(age, df = 4)1:DTG_pre.factorPre and during DTG"]*(sp[rownames(sp) %in% 35, 1])+ 
                                                                                  m2_summ$Coefs[m2_summ$Variables %in% "ns(age, df = 4)2:DTG_pre.factorPre and during DTG"]*(sp[rownames(sp) %in% 35, 2])+                                                          
                                                                                  m2_summ$Coefs[m2_summ$Variables %in% "ns(age, df = 4)3:DTG_pre.factorPre and during DTG"]*(sp[rownames(sp) %in% 35, 3])+                                                          
                                                                                  m2_summ$Coefs[m2_summ$Variables %in% "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]*(sp[rownames(sp) %in% 35, 4])), digits =3)

m2_exp_df$RR[m2_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'ade_typeYes']), digits =3)
m2_exp_df$RR[m2_exp_df$Variable %in% "TB: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'tbYes']), digits =3)
m2_exp_df$RR[m2_exp_df$Variable %in% "Brazil"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
m2_exp_df$RR[m2_exp_df$Variable %in% "Chile"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorChile']), digits =3)
m2_exp_df$RR[m2_exp_df$Variable %in% "Honduras"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorHonduras']), digits =3)


# function for relative risk


rr <- function(myage,model){
  exp(model$Coefs[model$Variables %in% "ns(age, df = 4)1"]*
        (sp[rownames(sp) %in% myage, 1]-sp[rownames(sp) %in% 35, 1]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)2"]*
        (sp[rownames(sp) %in% myage, 2]- sp[rownames(sp) %in% 35, 2]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)3"]*
        (sp[rownames(sp) %in% myage, 3]- sp[rownames(sp) %in% 35, 3]) +
        model$Coefs[model$Variables %in% "ns(age, df = 4)4"]*
        (sp[rownames(sp) %in% myage, 4]- sp[rownames(sp) %in% 35, 4]))
  
}

m2_exp_df$RR[m2_exp_df$Variable %in% 20] <- rr(20,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% 25] <- rr(25,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% 30] <- rr(30,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% 40] <- rr(40,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% 45] <- rr(45,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% 50] <- rr(50,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% 55] <- rr(55,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% 60] <- rr(60,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% 65] <- rr(65,m2_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "25 (Pre-During DTG warning)"] <- rr(25,m4_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "20 (Pre-During DTG warning)"] <- rr(20,m4_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "30 (Pre-During DTG warning)"] <- rr(30,m4_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "40 (Pre-During DTG warning)"] <- rr(40,m4_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "45 (Pre-During DTG warning)"] <- rr(45,m4_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "50 (Pre-During DTG warning)"] <- rr(50,m4_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "55 (Pre-During DTG warning)"] <- rr(55,m4_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "60 (Pre-During DTG warning)"] <- rr(60,m4_summ)
m2_exp_df$RR[m2_exp_df$Variable %in% "65 (Pre-During DTG warning)"] <- rr(65,m4_summ)





spline_se_fun <- function(myage,variance){
  
  tmp <- variance
  a <- sp$X1[rownames(sp) %in% myage]-sp$X1[rownames(sp) %in% 35]
  b <- sp$X2[rownames(sp) %in% myage]-sp$X2[rownames(sp) %in% 35]
  c <- sp$X3[rownames(sp) %in% myage]-sp$X3[rownames(sp) %in% 35]
  d <- sp$X4[rownames(sp) %in% myage]-sp$X4[rownames(sp) %in% 35]
  
  Var_x <- tmp["ns(age, df = 4)1", "ns(age, df = 4)1"]
  Var_y <- tmp["ns(age, df = 4)2", "ns(age, df = 4)2"]
  Var_z <- tmp["ns(age, df = 4)3", "ns(age, df = 4)3"]
  Var_w <- tmp["ns(age, df = 4)4", "ns(age, df = 4)4"]
  
  
  Cov_xy <- tmp["ns(age, df = 4)1", "ns(age, df = 4)2"]
  Cov_xz <- tmp["ns(age, df = 4)1", "ns(age, df = 4)3"]
  Cov_xw <- tmp["ns(age, df = 4)1", "ns(age, df = 4)4"]
  Cov_yz <- tmp["ns(age, df = 4)2", "ns(age, df = 4)3"]
  Cov_yw <- tmp["ns(age, df = 4)2", "ns(age, df = 4)4"]
  Cov_zw <- tmp["ns(age, df = 4)3", "ns(age, df = 4)4"]
  
  
  
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w 
  2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw 
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

spline_se_exfun <- function(myage,variance){
  
  tmp <- variance
  a <- sp$X1[rownames(sp) %in% myage]
  b <- sp$X2[rownames(sp) %in% myage]
  c <- sp$X3[rownames(sp) %in% myage]
  d <- sp$X4[rownames(sp) %in% myage]
  
  Var_x <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "ns(age, df = 4)1:DTG_pre.factorPre and during DTG"]
  Var_y <- tmp["ns(age, df = 4)2:DTG_pre.factorPre and during DTG", "ns(age, df = 4)2:DTG_pre.factorPre and during DTG"]
  Var_z <- tmp["ns(age, df = 4)3:DTG_pre.factorPre and during DTG", "ns(age, df = 4)3:DTG_pre.factorPre and during DTG"]
  Var_w <- tmp["ns(age, df = 4)4:DTG_pre.factorPre and during DTG", "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]
  var_v <- tmp["DTG_pre.factorPre and during DTG", "DTG_pre.factorPre and during DTG"]
  
  Cov_xy <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "ns(age, df = 4)2:DTG_pre.factorPre and during DTG"]
  Cov_xz <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "ns(age, df = 4)3:DTG_pre.factorPre and during DTG"]
  Cov_xw <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]
  Cov_yz <- tmp["ns(age, df = 4)2:DTG_pre.factorPre and during DTG", "ns(age, df = 4)3:DTG_pre.factorPre and during DTG"]
  Cov_yw <- tmp["ns(age, df = 4)2:DTG_pre.factorPre and during DTG", "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]
  Cov_zw <- tmp["ns(age, df = 4)3:DTG_pre.factorPre and during DTG", "ns(age, df = 4)4:DTG_pre.factorPre and during DTG"]
  Cov_wv <- tmp["ns(age, df = 4)4:DTG_pre.factorPre and during DTG","DTG_pre.factorPre and during DTG"]
  Cov_xv <- tmp["ns(age, df = 4)1:DTG_pre.factorPre and during DTG", "DTG_pre.factorPre and during DTG"]
  Cov_yv <- tmp["ns(age, df = 4)2:DTG_pre.factorPre and during DTG", "DTG_pre.factorPre and during DTG"]
  Cov_zv <- tmp["ns(age, df = 4)3:DTG_pre.factorPre and during DTG", "DTG_pre.factorPre and during DTG"]
  
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w + var_v +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw + 2*a*Cov_xv + 2*b*Cov_yv + 2*c*Cov_zv + 2*d*Cov_wv
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  
  print(pooled_se_splines_yr)
  
}

SE_20 <- spline_se_fun(20,covar)
SE_25 <-spline_se_fun(25,covar)
SE_30 <- spline_se_fun(30,covar)
SE_40 <-spline_se_fun(40,covar)
SE_45 <- spline_se_fun(45,covar)
SE_50 <-spline_se_fun(50,covar)
SE_55 <- spline_se_fun(55,covar)
SE_60 <-spline_se_fun(60,covar)
SE_65 <- spline_se_fun(65,covar)
SE_20_m4 <- spline_se_fun(20,covar4)
SE_25_m4 <-spline_se_fun(25,covar4)
SE_30_m4 <- spline_se_fun(30,covar4)
SE_40_m4 <-spline_se_fun(40,covar4)
SE_45_m4 <- spline_se_fun(45,covar4)
SE_50_m4 <-spline_se_fun(50,covar4)
SE_55_m4 <- spline_se_fun(55,covar4)
SE_60_m4 <-spline_se_fun(60,covar4)
SE_65_m4 <- spline_se_fun(65,covar4)
SE_35_DTG <- spline_se_exfun(35,covar)



# SE with sandwich estimates lower limits
m2_exp_df$Lower[m2_exp_df$Variable %in% "Pre and during DTG warning"] <- round(exp(log(m2_exp_df$RR[m2_summ$Variables %in% 'DTG_pre.factorPre and during DTG']) - 
                                                                                     1.96*SE_35_DTG), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'ade_typeYes'] - 
                                                                                           1.96*m2_summ$SE[m2_summ$Variables %in% 'ade_typeYes']), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "TB: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'tbYes'] - 
                                                                  1.96*m2_summ$SE[m2_summ$Variables %in% 'tbYes']), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "Brazil"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorBrazil'] - 
                                                                 1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "Chile"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorChile'] - 
                                                                1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorChile']), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "Honduras"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorHonduras'] - 
                                                                   1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 20] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 20]) - 
                                                           1.96*SE_20), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 25] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 25]) - 
                                                           1.96*SE_25), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 30] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 30]) - 
                                                           1.96*SE_30), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 40] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 40]) - 
                                                           1.96*SE_40), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 45] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 45]) - 
                                                           1.96*SE_45), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 50] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 50]) - 
                                                           1.96*SE_50), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 55] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 55]) - 
                                                           1.96*SE_55), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 60] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 60]) - 
                                                           1.96*SE_60), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% 65] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 65]) - 
                                                           1.96*SE_65), digits =3)

m2_exp_df$Lower[m2_exp_df$Variable %in% "20 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "20 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_20_m4), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "25 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "25 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_25_m4), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "30 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "30 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_30_m4), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "40 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "40 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_40_m4), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "45 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "45 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_45_m4), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "50 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "50 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_50_m4), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "55 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "55 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_55_m4), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "60 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "60 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_60_m4), digits =3)
m2_exp_df$Lower[m2_exp_df$Variable %in% "65 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "65 (Pre-During DTG warning)"]) - 
                                                                                      1.96*SE_65_m4), digits =3)





# SE with sandwich estimates Upper limits
m2_exp_df$Upper[m2_exp_df$Variable %in% "Pre and during DTG warning"] <- round(exp(log(m2_exp_df$RR[m2_summ$Variables %in% 'DTG_pre.factorPre and during DTG']) + 
                                                                                     1.96*SE_35_DTG), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'ade_typeYes'] + 
                                                                                           1.96*m2_summ$SE[m2_summ$Variables %in% 'ade_typeYes']), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "TB: Yes"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'tbYes'] + 
                                                                  1.96*m2_summ$SE[m2_summ$Variables %in% 'tbYes']), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "Brazil"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorBrazil'] + 
                                                                 1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorBrazil']), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "Chile"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorChile'] + 
                                                                1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorChile']), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "Honduras"] <- round(exp(m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorHonduras'] + 
                                                                   1.96*m2_summ$SE[m2_summ$Variables %in% 'study_site.factorHonduras']), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 20] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 20])  + 
                                                           1.96*SE_20), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 25] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 25])  + 
                                                           1.96*SE_25), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 30] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 30])  + 
                                                           1.96*SE_30), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 40] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 40])  + 
                                                           1.96*SE_40), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 45] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 45])  + 
                                                           1.96*SE_45), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 50] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 50])  + 
                                                           1.96*SE_50), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 55] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 55])  + 
                                                           1.96*SE_55), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 60] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 60])  + 
                                                           1.96*SE_60), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% 65] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% 65])  + 
                                                           1.96*SE_65), digits =3)

m2_exp_df$Upper[m2_exp_df$Variable %in% "20 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "20 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_20_m4), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "25 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "25 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_25_m4), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "30 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "30 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_30_m4), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "40 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "40 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_40_m4), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "45 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "45 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_45_m4), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "50 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "50 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_50_m4), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "55 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "55 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_55_m4), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "60 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "60 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_60_m4), digits =3)
m2_exp_df$Upper[m2_exp_df$Variable %in% "65 (Pre-During DTG warning)"] <- round(exp(log(m2_exp_df$RR[m2_exp_df$Variable %in% "65 (Pre-During DTG warning)"])  + 
                                                                                      1.96*SE_65_m4), digits =3)



# P value computation with sandwich estimates
z1 <- log(m2_exp_df$RR[m2_summ$Variables %in% 'DTG_pre.factorPre and during DTG'])/SE_35_DTG
z2 <- m2_summ$Coefs[m2_summ$Variables %in% 'ade_typeYes']/m2_summ$SE[m2_summ$Variables %in% 'ade_typeYes']
z3 <- m2_summ$Coefs[m2_summ$Variables %in% 'tbYes']/m2_summ$SE[m2_summ$Variables %in% 'tbYes']
z5 <- m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorBrazil']/m2_summ$SE[m2_summ$Variables %in% 'study_site.factorBrazil']
z6 <- m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorChile']/m2_summ$SE[m2_summ$Variables %in% 'study_site.factorChile']
z7 <- m2_summ$Coefs[m2_summ$Variables %in% 'study_site.factorHonduras']/m2_summ$SE[m2_summ$Variables %in% 'study_site.factorHonduras']
z8 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 20])/SE_25
z9 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 25])/SE_25
z10 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 30])/SE_30
z11 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 40])/SE_40
z12 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 45])/SE_45
z13 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 50])/SE_50
z14 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 55])/SE_55
z15 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 60])/SE_60
z16 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% 65])/SE_65
z17 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "20 (Pre-During DTG warning)"])/SE_20_m4
z18 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "25 (Pre-During DTG warning)"])/SE_25_m4
z19 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "30 (Pre-During DTG warning)"])/SE_30_m4
z20 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "40 (Pre-During DTG warning)"])/SE_40_m4
z21 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "45 (Pre-During DTG warning)"])/SE_45_m4
z22 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "50 (Pre-During DTG warning)"])/SE_50_m4
z23 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "55 (Pre-During DTG warning)"])/SE_55_m4
z24 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "60 (Pre-During DTG warning)"])/SE_60_m4
z25 <- log(m2_exp_df$RR[m2_exp_df$Variable %in% "65 (Pre-During DTG warning)"])/SE_65_m4



m2_exp_df$p[m2_exp_df$Variable %in% "Pre and during DTG warning"] <- round((1 - pnorm(abs(z1))) * 2,digits = 3)
m2_exp_df$p[m2_exp_df$Variable %in% "Other AIDS defining illness: Yes"] <- round((1 - pnorm(abs(z2))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "TB: Yes"] <- round((1 - pnorm(abs(z3))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "Brazil"] <- round((1 - pnorm(abs(z5))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "Chile"] <- round((1 - pnorm(abs(z6))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "Honduras"] <- round((1 - pnorm(abs(z7))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 20] <- round((1 - pnorm(abs(z8))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 25] <- round((1 - pnorm(abs(z9))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 30] <- round((1 - pnorm(abs(z10))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 40] <- round((1 - pnorm(abs(z11))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 45] <- round((1 - pnorm(abs(z12))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 50] <- round((1 - pnorm(abs(z13))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 55] <- round((1 - pnorm(abs(z14))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 60] <- round((1 - pnorm(abs(z15))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% 65] <- round((1 - pnorm(abs(z16))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "20 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z17))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "25 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z18))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "30 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z19))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "40 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z20))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "45 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z21))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "50 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z22))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "55 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z23))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "60 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z24))) * 2,digits =3)
m2_exp_df$p[m2_exp_df$Variable %in% "65 (Pre-During DTG warning)"] <- round((1 - pnorm(abs(z25))) * 2,digits =3)




save(fm_exp_df,file ="workspace/fm_exp_df.Rdata")
save(unfm,file ="workspace/unfm.Rdata")
save(int_f,file ="workspace/int_f.Rdata")
save(int_f_cat,file ="workspace/int_f_cat.Rdata")
save(overall_fm_dummy,file ="workspace/overall_fm_dummy.Rdata")
save(overall_fm_cat,file ="workspace/overall_fm_cat.Rdata")
save(m2_exp_df,file ="workspace/m2_exp_df.Rdata")
save(dem_female,file ="workspace/dem_female.Rdata")


## Sensitivity analysis for viral load only

a.h.vl <- all_haiti %>% filter(viral_failure.n == "No" & num_cat.factor == "One ART regiment before baseline")


## Model ##
m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor*gender.factor*period, data= a.h.vl, cluster = patient_id)

overa.h.vl_a <- data.matrix(m$var)
rownames(overa.h.vl_a) <- NULL
colnames(overa.h.vl_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL


#OMP
ind_test = 1
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))



#OFP
ind_test = c(1,2,4)
test_var = overa.h.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(1,2,3,4,5,6,7)
test_var = overa.h.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(1,3,5)
test_var = overa.h.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 2
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))

#YFPr
ind_test = c(2,3,6)
test_var = overa.h.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))


## young female pre vs young female post
ind_test = c(2,6)
test_var = overa.h.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post", "older female post","older female pre","older male pre","younger female post","younger female pre","younger male pre","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


h.vl_interaction  <-  data.frame(variable, rr,ll,ul,p)


mi  <- coxph(Surv(time0, time1,event) ~ age_gender_period, data= a.h.vl, cluster = patient_id)

overa.h.vl_a <- data.matrix(mi$var)
rownames(overa.h.vl_a) <- NULL
colnames(overa.h.vl_a) <- NULL
m1 <- as.data.frame(mi$coefficients)
rownames(m1) <- NULL
#OMP

ind_test = 3
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))

#OFP
ind_test = 1
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))


#OFPr
ind_test = 2
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))



#OMPr
ind_test = 4
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))

#YFP
ind_test = 5
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))


#YFPr
ind_test = 6
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))

#YMPr
ind_test = 7
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overa.h.vl_a[ind_test,ind_test]))

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3))


h.vl_i  <-  data.frame(variable, rr,ll,ul,p)



## Sensitivity analysis for other sites and viral failure 

a.o.vl <- all_other %>% filter(viral_failure.n == "No" & num_cat.factor == "One ART regiment before baseline")
## Model ##
m  <- coxph(Surv(time0, time1,event) ~ age_gt50.factor*gender.factor*period + strata(study_site.factor), data= a.o.vl, cluster = patient_id)

overa.o.vl_a <- data.matrix(m$var)
rownames(overa.o.vl_a) <- NULL
colnames(overa.o.vl_a) <- NULL
m1 <- as.data.frame(m$coefficients)
rownames(m1) <- NULL


#OMP
ind_test = 1
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))



#OFP
ind_test = c(1,2,4)
test_var = overa.o.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(se_test)

#OFPr
ind_test = c(1,2,3,4,5,6,7)
test_var = overa.o.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(se_test)


#OMPr
ind_test = c(1,3,5)
test_var = overa.o.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
OMPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(se_test)

#YFP
ind_test = 2
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))

#YFPr
ind_test = c(2,3,6)
test_var = overa.o.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(se_test)

#YMPr
ind_test = 3
test_beta = as.numeric(m1$`m$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))


## young female pre vs young female post
ind_test = c(2,6)
test_var = overa.o.vl_a[ind_test, ind_test]
se_test = sum(diag(test_var)) + 2*sum(test_var[upper.tri(test_var, diag = FALSE)])
test_beta = sum(as.numeric(m1$`m$coefficients`[ind_test]))
YFPr_YFP_UB_test = exp(test_beta + 1.96*sqrt(se_test))
YFPr_YFP_LB_test = exp(test_beta - 1.96*sqrt(se_test))
YFPr_YFP_point_test = exp(test_beta)
YFPr_YFP_LB_test
z8 <- test_beta/sqrt(se_test)

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre(ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)","young female pre vs young female post")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test,YFPr_YFP_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test,YFPr_YFP_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test,YFPr_YFP_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3),round((1 - pnorm(abs(z8))) * 2,digits=3))


a.o.vl.interaction  <-  data.frame(variable, rr,ll,ul,p)


mi  <- coxph(Surv(time0, time1,event) ~ age_gender_period + strata(study_site.factor), data= a.o.vl, cluster = patient_id)

overa.o.vl_a <- data.matrix(mi$var)
rownames(overa.o.vl_a) <- NULL
colnames(overa.o.vl_a) <- NULL
m1 <- as.data.frame(mi$coefficients)
rownames(m1) <- NULL
#OMP

ind_test = 3
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OMP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OMP_point_test = exp(test_beta)
z1 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))

#OFP
ind_test = 1
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OFP_point_test = exp(test_beta)
z2 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))


#OFPr
ind_test = 2
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OFPr_point_test = exp(test_beta)
z3 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))



#OMPr
ind_test = 4
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
OMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
OMPr_point_test = exp(test_beta)
z4 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))

#YFP
ind_test = 5
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFP_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YFP_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YFP_point_test = exp(test_beta)
z5 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))


#YFPr
ind_test = 6
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YFPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YFPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YFPr_point_test = exp(test_beta)
z6 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))

#YMPr
ind_test = 7
test_beta = as.numeric(m1$`mi$coefficients`[ind_test])
YMPr_UB_test = exp(test_beta + 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YMPr_LB_test = exp(test_beta - 1.96*sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test])))
YMPr_point_test = exp(test_beta)
z7 <- test_beta/sqrt(as.numeric(overa.o.vl_a[ind_test,ind_test]))

## Interaction term
variable <- c("older male post (ref = younger male post)", "older female post (ref = younger male post)","older female pre (ref = younger male post)","older male pre (ref = younger male post)","younger female post (ref = younger male post)","younger female pre (ref = younger male post)","younger male pre (ref = younger male post)")
rr <- c(OMP_point_test,OFP_point_test,OFPr_point_test,OMPr_point_test,YFP_point_test,YFPr_point_test,YMPr_point_test)
ll <- c(OMP_LB_test,OFP_LB_test,OFPr_LB_test,OMPr_LB_test,YFP_LB_test,YFPr_LB_test,YMPr_LB_test)
ul <- c(OMP_UB_test,OFP_UB_test,OFPr_UB_test,OMPr_UB_test,YFP_UB_test,YFPr_UB_test,YMPr_UB_test)

p <- c(round((1 - pnorm(abs(z1))) * 2,digits = 3),round((1 - pnorm(abs(z2))) * 2,digits =3),round((1 - pnorm(abs(z3))) * 2,digits=3),
       round((1 - pnorm(abs(z4))) * 2,digits =3),round((1 - pnorm(abs(z5))) * 2,digits=3),round((1 - pnorm(abs(z6))) * 2,digits =3),round((1 - pnorm(abs(z7))) * 2,digits=3))


a.o.vl.i  <-  data.frame(variable, rr,ll,ul,p)

save(h.vl_interaction, file="h.vl_interaction.Rdata")
save(h.vl_i, file="h.vl_i.Rdata")
save(a.o.vl.interaction, file="a.o.vl.interaction.Rdata")
save(a.o.vl.i, file="a.o.vl.i.Rdata")


### SUMMARY TABLES ###

VAR1 <- c("OMPost:YMPost", "OFPost:YMPost", "OFPre:YMP", "OMPre:YMPost","YFPost:YMPost", "YFPre:YMPost","YMPre:YMPost", "YFPre: YFPost")
VAR2 <- c("TB (ref = no)",  "Viral Failure (ref = No)" , 
          "YMPost:OMPost", "OFPost:YMPost", "OFPr:YMPost", "OMPre:YMPost","YFPost:YMPost", "YFPre:YMPost","YMPre:YMPost", 
          "YFPre: YFPost", "YFPost:OFPost", "OMPost:OFPost" , 
          "YFPre:YMPre", "OFPre:OMPre", "YFPre:OFPre","Number of ART regiment before baseline (ref = One)","YMPr:OMPr","OFPr:OFP","OMPr:OMP" )
VAR2_o <- c("TB (ref = no)",  "Viral Failure (ref = No)" , 
            "OMPost:YMPost", "OFPost:YMPost", "OFPr:YMPost", "OMPre:YMPost","YFPost:YMPost", "YFPre:YMPost","YMPre:YMPost", "YFPre: YFPost", "OFPost:YFPost", "OMPost:OFPost" , 
            "YFPre:YMPre", "OFPre:OMPre", "YFPre:OFPre","Number of ART regiment before baseline (ref = One)" )

VAR3 <- c("TB (ref = no)", "OMPost:YMPost", "OFPost:YMPost", "OFPr:YMPost", "OMPre:YMPost","YFPost:YMPost", "YFPre:YMPost","YMPre:YMPost", "YFPre: YFPost", 
          "OFPost:YFPost", "OMPost:OFPost" , "YFPre:YMPre", "OFPre:OMPre", "YFPre:OFPre","Number of ART regiment before baseline (ref = One)")

HR1 <- paste0(round(t.up.interaction$rr,digits =3),"(", round(t.up.interaction$ll,digits =3),"-",round(t.up.interaction$ul,digits =3),")")
HR2 <- paste0(round(t.up_a$rr,digits =3),"(", (round(t.up_a$ll,digits =3)),"-",round(t.up_a$ul,digits =3),")")
HR3 <- paste0(round(t.up_imp$rr,digits =3),"(", (round(t.up_imp$ll,digits =3)),"-",round(t.up_imp$ul,digits =3),")")
HR4 <- paste0(round(t.up_vl$rr,digits =3),"(", (round(t.up_vl$ll,digits =3)),"-",round(t.up_vl$ul,digits =3),")")

HR1_H <- paste0(round(t.hup.interaction$rr,digits =3),"(", round(t.hup.interaction$ll,digits =3),"-",round(t.hup.interaction$ul,digits =3),")")
HR2_H <- paste0(round(t.hup_a$rr,digits =3),"(", (round(t.hup_a$ll,digits =3)),"-",round(t.hup_a$ul,digits =3),")")
HR3_H <- paste0(round(t.hup_imp$rr,digits =3),"(", (round(t.hup_imp$ll,digits =3)),"-",round(t.hup_imp$ul,digits =3),")")
HR4_H <- paste0(round(t.hup_vl$rr,digits =3),"(", (round(t.hup_vl$ll,digits =3)),"-",round(t.hup_vl$ul,digits =3),")")

HR1_O <- paste0(round(t.other.up.interaction$rr,digits =3),"(", round(t.other.up.interaction$ll,digits =3),"-",round(t.other.up.interaction$ul,digits =3),")")
HR2_O <- paste0(round(t.other.up_a$rr,digits =3),"(", (round(t.other.up_a$ll,digits =3)),"-",round(t.other.up_a$ul,digits =3),")")
HR3_O <- paste0(round(t.other.up_imp$rr,digits =3),"(", (round(t.other.up_imp$ll,digits =3)),"-",round(t.other.up_imp$ul,digits =3),")")
HR4_O <- paste0(round(t.other.up_vl$rr,digits =3),"(", (round(t.other.up_vl$ll,digits =3)),"-",round(t.other.up_vl$ul,digits =3),")")


MOD1 <- data.frame(VAR1,HR1,t.up.interaction$p)
MOD2 <- data.frame(VAR2_o,HR2,t.up_a$p)
MOD3 <- data.frame(VAR2,HR3,t.up_imp$p)
MOD4 <- data.frame(VAR3,HR4,t.up_vl$p)

MOD1_H <- data.frame(HR1_H,t.hup.interaction$p)
MOD2_H <- data.frame(HR2_H,t.hup_a$p)
MOD3_H <- data.frame(HR3_H,t.hup_imp$p)
MOD4_H <- data.frame(HR4_H,t.hup_vl$p)

MOD1_O <- data.frame(HR1_O,t.other.up.interaction$p)
MOD2_O <- data.frame(HR2_O,t.other.up_a$p)
MOD3_O <- data.frame(HR3_O,t.other.up_imp$p)
MOD4_O <- data.frame(HR4_O,t.other.up_vl$p)

library(data.table) #data.table_1.9.5
TEST1 <- rbindlist(list(MOD1,MOD2,MOD3,MOD4))
TEST2 <- rbindlist(list(MOD1_H,MOD2_H,MOD3_H,MOD4_H))
TEST3 <- rbindlist(list(MOD1_O,MOD2_O,MOD3_O,MOD4_O))


colnames(TEST1) <- NULL
colnames(TEST2) <- NULL
colnames(TEST3) <- NULL



MODELS1 <- cbind(TEST1,TEST2,TEST3)


colnames(MODELS1) <- c("","HR", "p-value","HR", "p-value","HR", "p-value")

save(MODELS1, file="MODELS1.Rdata")
