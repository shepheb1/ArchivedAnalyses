## ----introduction2, include=FALSE-----------------------------------------------------------------------
#Clear workspace
unlink("~/OneDrive/CNICS205_20220912/code/Analysis_DM_2022_12_05_cache", recursive = TRUE)

rm(list = ls()) 
require(openxlsx)
source('~/Dropbox/xxxxfunction12182017.r', echo=TRUE)

#Set options
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(cache = FALSE)
knitr::opts_chunk$set(cache.lazy = FALSE)
knitr::opts_chunk$set(echo = FALSE)
#source("http://biostatdata.app.vumc.org/tgs/misc/r_functions.txt")
#Packages
Packages <- c("readr","purrr","tidyr","lubridate","stringr","docxtractr","kableExtra","knitr","data.table", "ggplot2","openxlsx","xtable","Hmisc","boot", "irr","utils","lme4","TeachingDemos","Hotelling","dplyr", "splines", "table1","tibble", "ggpubr" ,"nlme")
lapply(Packages, require, character.only = TRUE)


## ----include=FALSE--------------------------------------------------------------------------------------
SPRINTF <- function(x){
  ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
}
EST95CI = function(EST, LWR, UPR){
  paste(sprintf("%.2f", EST)
        , " ("
        , sprintf("%.2f", LWR)
        , ", "              
        , sprintf("%.2f", UPR)
        , ")"
        , sep = ""
  )
}


ribbon_gglot2_OR <- function(DATA, YLIMITS, SEL_VAR, XLAB, YLAB, OVERALLP, REF, CAPTION) {
  P =  DATA  %>%
    ggplot( aes(x=CONTVAR, y= point_exp)) +
    geom_line( size =1)+ 
    ylab(str_wrap(YLAB, width = 50)) +
    xlab(XLAB)  +
    # scale_x_continuous(breaks = sort(as.numeric(unique(modeldata$CONTVAR_20160101))), limits = c(2, 80) ) +
    scale_y_continuous(limits = YLIMITS
                       , trans = "log"
                       , breaks = c(0, 0.25,0.5,0.75,1,1.25, 2, 4, 8)) +
    # theme_classic() +
    geom_hline(yintercept = 1, colour = "gray")
  
  P1 = P +
    geom_vline(xintercept = REF, linetype = "dashed", size = 0.8, color = "gray") +
    scale_linetype_manual(values=c("solid","twodash", "dashed"))+
    geom_ribbon(aes(ymin = lwr_exp, ymax = upr_exp), show.legend=FALSE, alpha=0.20, fill = "#56B4E9") +
    geom_vline(xintercept = REF, colour = "gray") +
    #  geom_errorbar(aes(ymin= lwr_exp, ymax= upr_exp), colour="black", width=.2, position=pd) +
    geom_line(position=pd) +
    #  geom_point(position=pd, size=3) +
    #  scale_x_continuous( breaks = as.vector(c(Q1Q2,Q3))) +
    annotate(geom = 'text'
             , label = OVERALLP
             , x = -Inf
             , y = YLIMITS[2]
             , hjust = 0
             , vjust = 1
             , size = 5) +
    theme_bw()+
    theme(
      axis.title.x = element_text( size = 16),
      axis.title.y = element_text( size = 16),
      axis.text.y = element_text( size = 14),
      axis.text.x = element_text( size = 14),
      legend.position = "top"
      , legend.title  = element_text(colour="white", size =20)
      , legend.text   = element_text(size = 16)
    )+ 
    labs(caption = CAPTION)
}



## ----filedirectory, include=FALSE-----------------------------------------------------------------------
#Setting up location to save and read data
   file_path_fitted = "/home/xxxx/OneDrive/CNICS205_20220912/data/modelfit_temp/dm/"
  file_path_rawdata <- "/home/xxxx/OneDrive/CNICS205_20220912/data/data2023_09_13/"
file_path_cleandata <- "/home/xxxx/OneDrive/CNICS205_20220912/data/cleaned_data/"


## -------------------------------------------------------------------------------------------------------
file_path_doc = "/home/xxxx/OneDrive/CNICS205_20220912/documentLibrary"
setwd(file_path_doc)
      doclist <- list.files(".")[grepl(".doc", list.files("."))]


## -------------------------------------------------------------------------------------------------------
cleandatalist = readRDS(paste0(file_path_cleandata, "cleandatalist.rds"))


## -------------------------------------------------------------------------------------------------------
initialvisitdata <- cleandatalist[["Patient.csv"]] %>%
  dplyr::select(studyid, initialvisit)


## -------------------------------------------------------------------------------------------------------
rawdatalist = readRDS(paste0(file_path_cleandata, "/","rawdatalist.rds"))

tmp1 = rawdatalist[["PRO_Drugs_Assist.csv"]] %>%
  dplyr::select(studyid, date, drugever, drugever_p, drugcur, drugcur_p, drug3c, drug3c_p) %>%
  arrange(studyid, date)  %>%
  mutate(studyid = as.character(studyid)) %>%
  mutate(drugcur_p = ifelse(drugcur_p =="", NA, drugcur_p))%>%
  mutate(drugcur_p_date = ifelse(is.na(drugcur_p), NA, date)) %>%
  mutate(drugever_p = ifelse(drugever_p =="", NA, drugever_p)) %>%
  mutate(drugever_p_date = ifelse(is.na(drugever_p), NA, date)) %>%
  mutate(drug3c = ifelse(drug3c =="", NA, drug3c)) %>%
  mutate(drug3c_date = ifelse(is.na(drug3c), NA, date)) %>%
  mutate(drug3c_p = ifelse(drug3c_p =="", NA, drug3c_p)) %>%
  mutate(drug3c_p_date = ifelse(is.na(drug3c_p), NA, date)) %>%
  mutate(drugcur = ifelse(drugcur =="", NA, drugcur)) %>%
  mutate(drugcur_date = ifelse(is.na(drugcur), NA, date)) %>%
  mutate(drugever = ifelse(drugever =="", NA, drugever)) %>%
  mutate(drugever_date = ifelse(is.na(drugever), NA, date)) 

VARNAMES <- c("resultdate",grep("date", names(tmp1), value = TRUE))

cleandatalist[["drug"]] = tmp1 %>%
  dplyr::rename("resultdate" = date)%>%
  # mutate( drug3c_date = dmy(drug3c_date)) %>%
  mutate_if(names(.) %in% VARNAMES, dmy) %>%
  arrange(studyid, resultdate) 


## -------------------------------------------------------------------------------------------------------
tmp = cleandatalist$hba1c %>%
  full_join(cleandatalist$cd4, by = c("studyid", "resultdate")) %>%
  full_join(cleandatalist$hiv1rna, by = c("studyid", "resultdate"))%>%
  full_join(cleandatalist$depsum, by = c("studyid", "resultdate")) %>%
  full_join(cleandatalist[["drug"]], by = c("studyid", "resultdate"))%>%
  full_join(cleandatalist[["depression_diagnosis"]], by = c("studyid", "resultdate")) %>%
  full_join(cleandatalist[["bmi"]], by = c("studyid", "resultdate")) %>%
  full_join(cleandatalist[["alcohol"]], by = c("studyid", "resultdate")) %>%
  left_join(cleandatalist[["hba1c_treatment"]], by = c("studyid", "hba1c_date")) # getting DM and Depression treatment updating

tmp1 = tmp %>%
  right_join(cleandatalist$DM_diagnosis, by = "studyid") %>%
  dplyr::filter(!is.na(dm_diagnosis_date)) %>%
  group_by(studyid) %>%
  arrange(studyid, resultdate) %>%
  fill(cd4_date) %>%
  fill(cd4_num) %>%
  fill(hiv1rna_date) %>%
  fill(hiv1rna_detected) %>%
  fill(depsum_i) %>%
  fill(depsum_i_date) %>%
  fill(depsum_cc) %>%
  fill(depsum_cc_date) %>%
  fill(drugcur_p) %>%
  fill(drugcur_p_date) %>%
  fill(drugcur) %>%
  fill(drugcur_date) %>%
  fill(drugever_p) %>%
  fill(drugever_p_date) %>%
  fill(drugever) %>%
  fill(drugever_date) %>%
  fill(drug3c) %>%
  fill(drug3c_date) %>%
  fill(drug3c_p) %>%
  fill(drug3c_p_date) %>%
  fill(depression) %>%
  fill(bmi) %>%
  fill(bmi_date) %>%
  fill(alcriskhi) %>%
  fill(alcriskhi_date) %>%
  mutate(depression = replace_na(depression, 0)) %>%
  mutate(med_dep = replace_na(med_dep, 0)) %>%
  mutate(med_dm = replace_na(med_dm, 0)) %>%
  mutate(med_hiv = replace_na(med_hiv, 0)) %>%
  mutate(med_htn = replace_na(med_htn, 0)) 


CUTOFF = 456

tmp2 <- tmp1  %>%
  mutate(hiv1rna_detected = ifelse(hba1c_date-hiv1rna_date >= CUTOFF, NA, as.character(hiv1rna_detected))) %>%
  mutate(depsum_i = ifelse(hba1c_date-depsum_i_date >= CUTOFF, NA, depsum_i)) %>%
  mutate(depsum_cc = ifelse(hba1c_date-depsum_cc_date >= CUTOFF, NA, depsum_cc)) %>%
  mutate(cd4_num = ifelse(hba1c_date-cd4_date >= CUTOFF, NA, cd4_num)) %>%
  mutate(drugcur_p = ifelse(hba1c_date-drugcur_p_date >= CUTOFF, NA, drugcur_p)) %>%
  mutate(drugever_p = ifelse(hba1c_date-drugever_p_date >= CUTOFF, NA, drugever_p)) %>%
  mutate(drug3c = ifelse(hba1c_date-drug3c_date >= CUTOFF, NA, drug3c)) %>%
  mutate(drug3c_p = ifelse(hba1c_date-drug3c_p_date >= CUTOFF, NA, drug3c_p)) %>%
  mutate(alcriskhi = ifelse(hba1c_date-alcriskhi_date >= CUTOFF, NA, alcriskhi)) 


## -------------------------------------------------------------------------------------------------------
data4checking = tmp2 %>%
  full_join(initialvisitdata, by = "studyid") %>%
  dplyr::filter(hba1c_date >= dm_diagnosis_date) %>%
  dplyr::filter(hba1c_date >= initialvisit) %>%
  dplyr::filter(!is.na(depsum_i)) %>%
  mutate_at(vars(ends_with("_date")), funs(duration = as.numeric(hba1c_date- . )))


## ----fig.width=16, fig.height=11------------------------------------------------------------------------
data4checking %>%
  ungroup() %>%
  dplyr::select(-hba1c_date_duration, -depression_diagnosis_date_duration, -dm_diagnosis_date_duration) %>%
  dplyr::select(contains("duration")) %>%
  gather() %>%
  mutate(value = as.numeric(as.character(value))) %>%
  dplyr::filter(!is.na(value)) %>%
#  mutate(value = as.numeric(value)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ key, scales = "free") +
  theme_bw()


## ----eval=FALSE, include=FALSE--------------------------------------------------------------------------
# diffdata = tmp2 %>%
#   dplyr::select(contains("cd4"), studyid, contains("hba1c"), contains("date")) %>%
#   mutate(diff = hba1c_date - cd4_date)
# diffdata %>%
#  # dplyr::filter(diff > 16*30) %>%
#   dplyr::select(cd4_num, cd4_date, hba1c_date, diff)


## ----include=FALSE--------------------------------------------------------------------------------------
tmp3 = tmp2 %>%
##  dplyr::filter(hba1c_date > dm_diagnosis_date) %>%
  left_join(cleandatalist$Patient.csv, by  ="studyid") %>%
  mutate(dm_control = ifelse(hba1c <= 7, 0, 1)) %>%
  mutate(age_updating = as.numeric(hba1c_date-dob)/365.25)  %>%
  mutate(sqrt_cd4_num = sqrt(cd4_num)) %>%
  mutate(time_since_dm_diagnosis = as.numeric(hba1c_date-dm_diagnosis_date)/365.25) %>%
  mutate(year_hba1c = year(hba1c_date))%>%
  mutate(year_depsum_i = year(depsum_i_date))%>%
  mutate(year_depsum_i = as.numeric(year_depsum_i))%>%
  mutate(year_depsum_cc = year(depsum_cc_date))%>%
  mutate(year_depsum_cc = as.numeric(year_depsum_cc))%>%
  mutate(hiv1rna_detected= factor(hiv1rna_detected, levels = c("<48", "[48, 400)", ">=400"))) %>%
  mutate(hiv1rna_detected2= ifelse(hiv1rna_detected %in% c("<48", "[48, 400)"), "<400", ">=400")) %>%
  mutate(transgender_model = ifelse(transgender %in% c("ftm transgender", "mtf transgender"), "transgendered", transgender)) %>%
  mutate(time_since_cnics_entry = as.numeric(hba1c_date-initialvisit)/365.25) %>%
  mutate(depsum_i_cat = case_when(
          depsum_i %in% (0:4) ~ "none (0-4)"
         , depsum_i %in% (5:9) ~ "mild (5-9)"
         , depsum_i %in% (10:14) ~ "moderate (10-14)"
         , depsum_i %in% (15:19) ~ "moderately severe (15-19)"
         , depsum_i >= 20 ~ "severe (>= 20)"
         , TRUE ~ ""
         )) %>%
  mutate(depsum_cc_cat = case_when( 
          depsum_cc %in% (0:4) ~ "none (0-4)"
         , depsum_cc %in% (5:9) ~ "mild (5-9)"
         , depsum_cc %in% (10:14) ~ "moderate (10-14)"
         , depsum_cc %in% (15:19) ~ "moderately severe (15-19)"
         , depsum_cc >= 20 ~ "severe (>= 20)"
         , TRUE ~ ""
         )) %>%
  mutate(depsum_cc_cat = factor(depsum_cc_cat, levels = c("none (0-4)"
                                                         , "mild (5-9)"
                                                         , "moderate (10-14)"
                                                         , "moderately severe (15-19)"
                                                         , "severe (>= 20)"))) %>%
  mutate(depsum_i_cat = factor(depsum_i_cat, levels = c("none (0-4)"
                                                         , "mild (5-9)"
                                                         , "moderate (10-14)"
                                                         , "moderately severe (15-19)"
                                                         , "severe (>= 20)"))) 



saveRDS(tmp3, paste0(file_path_cleandata, "/","data4dm_consort.rds")) ## Save data for DM consort diagram
cols_as_factors = c("race_hispanic_model"
                    , "riskfactor_model"
                    , "birthsex"
                    , "race_model"
                    , "drug3c_p"
                    , "drugever_p"
                    , "drugcur_p"
                    , "site"
                    , "alcriskhi"
                    )


## ----include=FALSE--------------------------------------------------------------------------------------
modeldata = tmp3 %>%
  dplyr::filter(hba1c_date > dm_diagnosis_date) %>%
  dplyr::filter(time_since_cnics_entry >= 0) %>%
  dplyr::filter(!is.na(depsum_i) | !is.na(depsum_cc)) %>%
  dplyr::filter(!is.na(hba1c)) %>%
  dplyr::filter(!is.na(cd4_num)) %>%
  dplyr::filter(!is.na(bmi)) %>%
  dplyr::filter(!is.na(alcriskhi)) %>%
  dplyr::filter(!is.na(hiv1rna_detected)) %>%
  dplyr::filter(!is.na(drugcur_p) | !is.na(drugever_p))  %>%
  group_by(studyid) %>%
  arrange(studyid, hba1c_date) %>%
  mutate(time_since_first_hba1c = as.numeric(hba1c_date - min(hba1c_date)))%>% 
  mutate_at(cols_as_factors, factor) %>%
  mutate(med_dep = factor(med_dep, levels=c(0,1), labels=c("No", "Yes")))%>%
  mutate(med_dm = factor(med_dm, levels=c(0,1), labels=c("No", "Yes")))%>%
  mutate(med_hiv = factor(med_hiv, levels=c(0,1), labels=c("No", "Yes")))%>%
  mutate(depression = factor(depression, levels=c(0,1), labels=c("No", "Yes")))%>%
  mutate(dm_control_grp = factor(dm_control, levels=c(0,1), labels=c("<=7%", ">7%")))%>%
  mutate(year_depsum_i_cat = cut(year_depsum_i
                                 , breaks = c(2005,2010,2016,2023)
                                 , include.lowest = TRUE
                                 )) %>%
  mutate(year_depsum_cc_cat = cut(year_depsum_cc
                                 , breaks = c(2005,2010,2016,2023)
                                 , include.lowest = TRUE
  )) %>%
  mutate(year_hba1c_cat = cut(year_hba1c
                                 , breaks = c(2005,2010,2016,2023)
                                 , include.lowest = TRUE
  )) %>%
  mutate(hiv1rna_detected = factor(hiv1rna_detected2))%>%
  mutate(diff_hba1c_phq = as.numeric(hba1c_date - depsum_i_date) ) 


## -------------------------------------------------------------------------------------------------------
baseline_dm <- modeldata %>%
  arrange(studyid, age_updating) %>%
  group_by(studyid) %>%
  mutate(followup_duration = max(time_since_dm_diagnosis) - min(time_since_dm_diagnosis)) %>%
  mutate(number_of_outcome_measures = n()) %>%
  slice(1)


## -------------------------------------------------------------------------------------------------------

dm12data= rawdatalist[["Patient.csv"]] %>%
  dplyr::filter(diabetes == "yes") %>%
  dplyr::select(studyid, dmtype) %>%
  mutate(studyid = as.character(studyid))

describe(modeldata$studyid)
modeldata %>%
  dplyr::distinct(studyid) %>%
  left_join(dm12data) %>%
  ungroup() %>%
  count(dmtype)



## ----include=FALSE--------------------------------------------------------------------------------------
tmp3 %>%
  dplyr::filter(!is.na(hba1c)) %>%
  dplyr::select(studyid, hba1c, hba1c_date, depsum_i, depsum_cc_date, initialvisit) %>%
  dplyr::filter(year(hba1c_date) == 1996) %>%
  dplyr::filter(!is.na(depsum_i) )


## -------------------------------------------------------------------------------------------------------
                Hmisc::label(modeldata$cd4_num) = "CD4 count"
           Hmisc::label(modeldata$sqrt_cd4_num) = "CD4 count (square root transformation)"
Hmisc::label(modeldata$time_since_dm_diagnosis) = "Time since diabetes diagnosis (years)"
 Hmisc::label(modeldata$time_since_cnics_entry) = "Time since CNICS clinic entry (years)"
               Hmisc::label(modeldata$depsum_i) = "Depression Symptoms (depsum_i)"
              Hmisc::label(modeldata$depsum_cc) = "Depression Symptoms (depsum_cc)"
           Hmisc::label(modeldata$depsum_i_cat) = "Depression Symptoms (depsum_i, categorical)"
          Hmisc::label(modeldata$depsum_cc_cat) = "Depression Symptoms (depsum_cc, categorical)"
       Hmisc::label(modeldata$hiv1rna_detected) = "HIV1 virus detectable levels"
         Hmisc::label(modeldata$year_depsum_cc) = "Calendar year of depsum_cc measures"
          Hmisc::label(modeldata$year_depsum_i) = "Calendar year of depsum_i measures"
     Hmisc::label(modeldata$year_depsum_cc_cat) = "Calendar year of depsum_cc measures"
      Hmisc::label(modeldata$year_depsum_i_cat) = "Calendar year of depsum_i measures"
         Hmisc::label(modeldata$year_hba1c_cat) = "Calendar year of hba1c"
           Hmisc::label(modeldata$age_updating) = "Age (years)"
                 Hmisc::label(modeldata$drug3c) = "3 levels: never, past or current use (excluding marijuana)"
               Hmisc::label(modeldata$drug3c_p) = "3 levels: never, past or current use (including marijuana)"
             Hmisc::label(modeldata$drugever_p) = "Binary: ever use drugs, including marijuana" 
               Hmisc::label(modeldata$drugever) = "Binary: ever use drugs, excluding marijuana" 
              Hmisc::label(modeldata$drugcur_p) = "Binary: current use of drugs, including marijuana"
                Hmisc::label(modeldata$drugcur) = "Binary: current use drugs, excluding marijuana" 
                   Hmisc::label(modeldata$race) = "Race (raw)"
                   Hmisc::label(modeldata$site) = "Site"
             Hmisc::label(modeldata$race_model) = "Race (modified)"
               Hmisc::label(modeldata$birthsex) = "Birth Sex"
    Hmisc::label(modeldata$race_hispanic_model) = "Race/ethnicity"
               Hmisc::label(modeldata$hispanic) = "Hispanic"
       Hmisc::label(modeldata$riskfactor_model) = "HIV acquisition risk factors"
             Hmisc::label(modeldata$depression) = "Documented depression diagnosis" 
            Hmisc::label(modeldata$transgender) = "Transgender"
      Hmisc::label(modeldata$transgender_model) = "Transgender (modified)"
             Hmisc::label(modeldata$hbv_status) = "HBV status"
             Hmisc::label(modeldata$hcv_status) = "HCV status"
             Hmisc::label(modeldata$dm_control) = "Hemoglobin A1C: > 7%"
         Hmisc::label(modeldata$dm_control_grp) = "Hemoglobin A1C: > 7%"
                  Hmisc::label(modeldata$hba1c) = "Hemoglobin A1C"
                    Hmisc::label(modeldata$bmi) = "BMI (m2/kg)"
              Hmisc::label(modeldata$alcriskhi) = "Alcohol use (alcriskhi)"
                 Hmisc::label(modeldata$med_dm) = "On diabetes medication (updating)"
                Hmisc::label(modeldata$med_dep) = "On depression medication (updating)"
                Hmisc::label(modeldata$med_hiv) = "On HIV medication (updating)"
         Hmisc::label(modeldata$diff_hba1c_phq) = "Time between hba1c measurement and PHQ collection date (days)"
                    


## -------------------------------------------------------------------------------------------------------
baseline_dm <- modeldata %>%
  arrange(studyid, age_updating) %>%
  group_by(studyid) %>%
  mutate(followup_duration = max(time_since_dm_diagnosis) - min(time_since_dm_diagnosis)) %>%
  mutate(number_of_outcome_measures = n()) %>%
  slice(1)

table1data = baseline_dm %>%
  mutate(year_depsum_cc = factor(as.numeric(year_depsum_cc), levels = 2005:2023)) %>%
  mutate(year_depsum_i = factor(as.numeric(year_depsum_i), levels = 2005:2023)) %>%
  mutate(depression = factor(depression)) %>%
  ungroup()

table1vars <- c( "age_updating" 
                 , "followup_duration"
                 , "number_of_outcome_measures"
                 , "birthsex"
                 , "alcriskhi"
                 , "transgender"
                 , "transgender_model"
                 , "race"
                 , "race_model"
                 , "hispanic"
                 , "race_hispanic_model"
                 , "riskfactor"
                 , "riskfactor_model"
                 , "hbv_status"
                 , "hcv_status"
                 , "time_since_dm_diagnosis"
                 , "cd4_num"
                 , "hiv1rna_detected"
                 , "depression"
             #    , "year_depsum_i"
            #     , "year_depsum_cc"
                 , "time_since_cnics_entry"
                 , "drugcur_p"
                 , "drugcur"
                 , "drugever_p"
                 , "drugever"
                 , "drug3c"
                 , "drug3c_p"
                 , "depsum_i"
                 , "depsum_i_cat"
                 , "depsum_cc"
                 , "depsum_cc_cat"
                 , "bmi"
                 , "hba1c"
                 , "year_hba1c_cat"
                 , "dm_control_grp"
                 , "med_dm"
                 , "med_dep"
                 , "med_hiv"
            , "diff_hba1c_phq"
                 )
GROUPVAR = "site"
    TEST = FALSE
 OVERALL = TRUE
 
       Hmisc::label(table1data$followup_duration) = "Follow up duration (years)"
Hmisc::label(table1data$number_of_outcome_measures) = "Number of measures"



## ----SUMMARYTABLE1TS------------------------------------------------------------------------------------
CAPTION = paste0("Selected characteristics by", Hmisc::label(table1data[,GROUPVAR]), "among Diabetic Patients at Baseline")

table1vars_man <- c("age_updating", "birthsex", "transgender_model", "riskfactor_model" ,"med_hiv","hiv1rna_detected", "time_since_dm_diagnosis", "cd4_num" ,"bmi","hba1c","dm_control_grp","med_dm","year_hba1c_cat","depsum_i_cat","depression","med_dep","hbv_status", "hcv_status","alcriskhi", "time_since_cnics_entry", "followup_duration", "diff_hba1c_phq"   )
tbl1_formula <- as.formula(paste(paste(table1vars, collapse = "+") ,  "~", GROUPVAR))

tbl1 <- summaryM(
  tbl1_formula
  , data = table1data
  , quant = c(0.25,0.5,0.75)
  , continuous = 8
  , overall = OVERALL
  , na.include = TRUE
  , test = TRUE
)

html(tbl1
     , caption = ""
     , vnames = "label"
     , exclude1 = FALSE
     , digits= 5
     , size='normal'
     , prmsd = FALSE
     , middle.bold=TRUE
     , long= TRUE
     , rowsep = TRUE)


## -------------------------------------------------------------------------------------------------------
table1data = table1data %>%
  mutate(dep9 = ifelse(depsum_i> 9, ">9", "<=9"))

  GROUPVAR = "race_hispanic_model"
      TEST = FALSE
   OVERALL = TRUE
table1vars = setdiff(table1vars, c("transgender"
                                   , "race", "race_model"
                                   ,  "hispanic", "riskfactor"
                                  # , "depression"
                                   , "depsum_cc", "depsum_cc_cat"
                                   , "drug3c" , "drug3c_p" 
                                   ) )
table1vars <- c(table1vars, "dep9")



## ----ref.label = "SUMMARYTABLE1TS"----------------------------------------------------------------------


## -------------------------------------------------------------------------------------------------------
##  Poster: Baseline characteristics by `r Hmisc::label(table1data[,GROUPVAR])` [format 2]

table1data = table1data %>%
  mutate(dep9 = ifelse(depsum_i> 9, ">9", "<=9"))

  GROUPVAR = "race_hispanic_model"
      TEST = FALSE
   OVERALL = TRUE
table1vars = setdiff(table1vars, c("transgender"
                                   , "race", "race_model"
                                   ,  "hispanic", "riskfactor"
                                  # , "depression"
                                   , "depsum_cc", "depsum_cc_cat"
                                   , "drug3c" , "drug3c_p" 
                                   ) )
table1vars <- c("age_updating", "birthsex", "cd4_num", "hiv1rna_detected", "time_since_dm_diagnosis"
                , "med_dm", "hba1c", "depression", "med_dep", "dep9" ,"followup_duration" ,"number_of_outcome_measures")



## ----ref.label = "SUMMARYTABLE1TS"----------------------------------------------------------------------


## ----SUMMARYTABLE1TS_TG---------------------------------------------------------------------------------
library(tgsify)
#devtools::install_github("thomasgstewart/tgsify")

tbl1 <- summaryM_to_df(
  tbl1
  , long = TRUE
  , exclude1 = FALSE
#  , prob = c(0.1,0.5,0.9)
  , what = "%"
  , digits = 2
  , prn = TRUE
  , vnames = "labels"
)

for (i in 3: (dim(tbl1)[2]-0)){
  tbl1[,i] <- slice_to_box(tbl1[,i])
  tbl1[,i] <- slice_to_box2(tbl1[,i])
  
}
names(tbl1)[1] <- "Variables"
#kable(tbl1)
kable(tbl1, "pandoc") %>%
  column_spec(1, width = "10em")

write.csv(tbl1, "/home/xxxx/OneDrive/CNICS205_20220912/writing/table1.csv")


## -------------------------------------------------------------------------------------------------------
GROUPVAR = "depsum_i_cat"
    TEST = TRUE
 OVERALL = TRUE


## ----ref.label = "SUMMARYTABLE1TS"----------------------------------------------------------------------


## ----ref.label = "SUMMARYTABLE1TS_TG"-------------------------------------------------------------------


## ----ONLYTABLE1, eval=FALSE-----------------------------------------------------------------------------
# knitr::knit_exit()


## ----fig.width=16, fig.height=11------------------------------------------------------------------------
table1data %>%
  dplyr::select("cd4_num"
                , "depsum_i"
                , "depsum_cc"
                , "age_updating"
                , "time_since_cnics_entry"
                , "time_since_dm_diagnosis"
                , "hba1c"
                , "bmi"
                ) %>%
  gather() %>%
  mutate(value = as.numeric(as.character(value))) %>%
  dplyr::filter(!is.na(value)) %>%
#  mutate(value = as.numeric(value)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ key, scales = "free")


## ----fig.width=16, fig.height=11------------------------------------------------------------------------
table1data %>%
  dplyr::select("cd4_num"
                , "depsum_i"
                , "depsum_cc"
                , "age_updating"
                , "time_since_cnics_entry"
                , "time_since_dm_diagnosis"
                , "hba1c"
                , "bmi"
                , "dm_control"
  ) %>%
  gather("variable", "value", -dm_control) %>%
  mutate(value = as.numeric(as.character(value))) %>%
  dplyr::filter(!is.na(value)) %>%
  mutate(dm_control_grp = ifelse(dm_control==0, "<=7%", ">7%")) %>%
  ggplot( aes(x = dm_control_grp, y = value)) + 
  geom_boxplot(outlier.shape=8, outlier.size = 1)+
  facet_wrap(~ variable, scales = "free")


## ----echo=FALSE-----------------------------------------------------------------------------------------
CATVARS_ALL <- c("race_model"
                 , "race_hispanic_model"
                 , "hiv1rna_detected"
                 , "drug3c_p"
                 , "drugever_p"
                 , "drugcur_p"
                 , "year_depsum_cc_cat"
                 , "year_depsum_i_cat"
                 , "year_hba1c_cat"
                 , "riskfactor_model"
                 , "birthsex"
                 , "site"
                 , "hcv_status"
                 , "hbv_status"
                 , "alcriskhi"
                 , "depsum_i_cat"
                 , "depsum_cc_cat"
                 , "depression"
                 , "med_dm"
              #   , "med_dep"
                 )

CONTVARS_ALL <- c("depsum_i"
                  , "depsum_cc"
                  , "age_updating"
                  , "time_since_dm_diagnosis"
                  , "time_since_cnics_entry"
               #   , "year_depsum_cc"
                #  , "year_depsum_i"
                  , "sqrt_cd4_num"
                  , "bmi"
                  )



CONTVARS_FMLA_ALL = paste0("ns(", CONTVARS_ALL,  ", 3)")
   
continuous_all <- data.frame(
    "CONTVARS_ALL" = CONTVARS_ALL
  , "CONTVARS_FMLA_ALL" = CONTVARS_FMLA_ALL
  , stringsAsFactors = FALSE
)



## ----echo=FALSE-----------------------------------------------------------------------------------------

ADJ_CONVARS <- c("age_updating"
                 , "time_since_dm_diagnosis" 
                 , "sqrt_cd4_num"
                 , "bmi"
                 )

CATVARS <- c("race_hispanic_model"
             , "hiv1rna_detected"
             , "riskfactor_model"
             , "birthsex"
             , "site"
             , "drugcur_p"
             , "year_hba1c_cat"
             , "depression"
             , "alcriskhi"
             , "med_dm"
             , "med_dep"
                 )

EXPOSUREALL <-  c("depsum_i"
                  , "depsum_cc"
                 )  
variables_INTER_TERMS = c("depsum_i* race_hispanic_model"
                 , "depsum_cc* race_hispanic_model"
                 , "depsum_i_cat* race_hispanic_model"
                 , "depsum_cc_cat* race_hispanic_model"
                 , "depsum_i* depression"
                 , "depsum_cc* depression"
                 , "depsum_i_cat* depression"
                 , "depsum_cc_cat* depression"
                 , "depsum_i* med_dm"
                 , "depsum_cc* med_dm"
                 , "depsum_i_cat* med_dm"
                 , "depsum_cc_cat* med_dm"      
                 , "depsum_i* med_dep"
                 , "depsum_cc* med_dep"
                 , "depsum_i_cat* med_dep"
                 , "depsum_cc_cat* med_dep"                         
                 )

INTER_TERMS <- c("ns(depsum_i,3)* race_hispanic_model"
                 , "ns(depsum_cc, 3)* race_hispanic_model"
                 , "depsum_i_cat* race_hispanic_model"
                 , "depsum_cc_cat* race_hispanic_model"
                 , "ns(depsum_i,3)* depression"
                 , "ns(depsum_cc,3)* depression"
                 , "depsum_i_cat* depression"
                 , "depsum_cc_cat* depression"
                 , "ns(depsum_i,3)* med_dm"
                 , "ns(depsum_cc, 3)* med_dm"
                 , "depsum_i_cat* med_dm"
                 , "depsum_cc_cat* med_dm"  
                 , "ns(depsum_i,3)* med_dep"
                 , "ns(depsum_cc, 3)* med_dep"
                 , "depsum_i_cat* med_dep"
                 , "depsum_cc_cat* med_dep"                        
           )


covariates_all = data.frame(
  variables = c(CONTVARS_ALL, CATVARS_ALL, variables_INTER_TERMS)
  , forms  = c(CONTVARS_FMLA_ALL, CATVARS_ALL, INTER_TERMS)
  , stringsAsFactors = FALSE
 
)


## -------------------------------------------------------------------------------------------------------
refdatalist <- NULL
refdatalist[["depsum_i"]] = list("REF" = 0
                                 , "rangetocalculate" = sort(unique(modeldata$depsum_i, na.rm = TRUE))
                               #  , "tabledata" =  setdiff(sort(unique(modeldata$depsum_i, na.rm = TRUE)), c(0))
                                 , "tabledata" =  setdiff(c(4,9,10, 14,15,19,20, 25), c(0))
                                 ,  "KNOTS" = attributes(ns(modeldata$depsum_i, 3))$knots
                                 , "BOUNDARY.KNOTS" = attributes(ns(modeldata$depsum_i, 3))$Boundary.knots
                                   )

refdatalist[["depsum_cc"]] = list("REF" = 0
                                  , "rangetocalculate" = sort(unique(modeldata$depsum_cc, na.rm = TRUE))
                                  , "tabledata" =  setdiff(sort(unique(modeldata$depsum_cc, na.rm = TRUE)), c(0))
                                 ,  "KNOTS" = attributes(ns(modeldata$depsum_cc, 3))$knots
                                 , "BOUNDARY.KNOTS" = attributes(ns(modeldata$depsum_cc, 3))$Boundary.knots
                                  )

refdatalist[["age_updating"]] = list("REF" = 45
                                   , "rangetocalculate" = 20:70
                                   , "tabledata" =  setdiff(c(25,35,55,65), c(45))
                                 ,  "KNOTS" = attributes(ns(modeldata$age_updating, 3))$knots
                                 , "BOUNDARY.KNOTS" = attributes(ns(modeldata$age_updating, 3))$Boundary.knots
                                   )

refdatalist[["bmi"]] = list("REF" = 25
                            , "rangetocalculate" = 200:500/10
                            , "tabledata" =  c(20,30,35,40)
                            , "KNOTS" = attributes(ns(modeldata$bmi, 3))$knots
                            , "BOUNDARY.KNOTS" = attributes(ns(modeldata$bmi, 3))$Boundary.knots
                            )

refdatalist[["time_since_dm_diagnosis"]] = list("REF" = 0
                                                , "rangetocalculate" = 1:180/10
                                                , "tabledata" =  c(0,2,4,8,10)
                                                , "KNOTS" = attributes(ns(modeldata$time_since_dm_diagnosis, 3))$knots
                                                , "BOUNDARY.KNOTS" = attributes(ns(modeldata$time_since_dm_diagnosis, 3))$Boundary.knots
                                   )

refdatalist[["time_since_cnics_entry"]] = list("REF" = 0
                                   , "rangetocalculate" = 1:210/10
                                   , "tabledata" =   c(0,2,4,8,10)
                                   )

refdatalist[["year_depsum_i"]] = list("REF" = min(modeldata$year_depsum_i)
                                   , "rangetocalculate" = sort(unique(modeldata$year_depsum_i, na.rm = TRUE))
                                   , "tabledata" =  setdiff(sort(unique(modeldata$year_depsum_i, na.rm = TRUE)), c(min(modeldata$year_depsum_i)))
                                   )

refdatalist[["year_depsum_cc"]] = list("REF" = min(modeldata$year_depsum_cc)
                                   , "rangetocalculate" = sort(unique(modeldata$year_depsum_cc, na.rm = TRUE))
                                   , "tabledata" =  setdiff(sort(unique(modeldata$year_depsum_cc, na.rm = TRUE)), c(min(modeldata$year_depsum_cc)))
                                   )

refdatalist[["sqrt_cd4_num"]] = list("REF" = 20
                                     , "rangetocalculate" = sqrt(0:1500)
                                     , "tabledata" =  setdiff(sqrt(c(100,225,625,1225)), c(20))
                                     ,  "KNOTS" = attributes(ns(modeldata$sqrt_cd4_num, 3))$knots
                                     , "BOUNDARY.KNOTS" = attributes(ns(modeldata$sqrt_cd4_num, 3))$Boundary.knots
                                   )

saveRDS(refdatalist, paste0(file_path_fitted, "REFDATALIST.rds"))


## ----echo=FALSE-----------------------------------------------------------------------------------------
rawmodeldata = modeldata


## -------------------------------------------------------------------------------------------------------
mod_rawmodeldata <- rawmodeldata %>%
  group_by(studyid) %>%
  dplyr::mutate(base_age_model = min(age_updating))%>%
  dplyr::mutate(base_year_hba1c_model = min(year_hba1c)) %>%
  mutate(base_year_hba1c_model_cat = cut(base_year_hba1c_model
                              , breaks = c(2005,2010,2016,2023)
                              , include.lowest = TRUE
  ))


## -------------------------------------------------------------------------------------------------------
griddata = expand.grid('depsum_i' = sort(unique(mod_rawmodeldata$depsum_i, na.rm = TRUE))
                       ,  'race_hispanic_model'= levels(mod_rawmodeldata$race_hispanic_model)
) %>%
  as.data.frame()


## ----eval=FALSE, include=FALSE--------------------------------------------------------------------------
# fitlist <- NULL
# fitlist[["fit_lme_no_inter"]]  = lme(hba1c ~ ns(depsum_i, 3)
#             + ns(base_age_model, 3)
#             + ns(time_since_dm_diagnosis, 3)
#             + ns(sqrt_cd4_num, 3)
#             + ns(bmi, 3)
#             + race_hispanic_model
#             + hiv1rna_detected
#             + drugcur_p
#             + base_year_hba1c_model_cat
#             + riskfactor_model
#             + birthsex
#             + site
#             + alcriskhi
#             + ns(depsum_i, 3)
#             , data = as.data.frame(mod_rawmodeldata)
#             , random = ~1|studyid
#             , control = lmeControl(opt = "optim")
#             #  , correlation=corCAR1(form = ~time_since_first_hba1c | studyid)
# )
# 
# car::Anova(fitlist[["fit_lme_no_inter"]], type = "III", test = "LR")
# car::Anova(fitlist[["fit_lme_no_inter"]], type = "III")
# 
# fitlist[["fit_lme_inter"]] = lme(hba1c ~ ns(depsum_i, 3)
#                           + ns(base_age_model, 3)
#                           + ns(time_since_dm_diagnosis, 3)
#                           + ns(sqrt_cd4_num, 3)
#                           + ns(bmi, 3)
#                           + race_hispanic_model
#                           + hiv1rna_detected
#                           + drugcur_p
#                           + base_year_hba1c_model_cat
#                           + riskfactor_model
#                           + birthsex
#                           + site
#                           + alcriskhi
#                           + ns(depsum_i, 3) * race_hispanic_model
#                           , data = as.data.frame(mod_rawmodeldata)
#                           , random = ~1|studyid
#                           , control = lmeControl(opt = "optim")
#                           #  , correlation=corCAR1(form = ~time_since_first_hba1c | studyid)
# )
# 
# car::Anova(fitlist[["fit_lme_inter"]], type = "III")
# 
# 
# fitlist[["fit_glmer"]] = glmer(dm_control ~ ns(depsum_i, 3)
#       + ns(base_age_model, 3)
#       + ns(time_since_dm_diagnosis, 3)
#       + ns(sqrt_cd4_num, 3)
#       + ns(bmi, 3)
#       + race_hispanic_model
#       + hiv1rna_detected
#       + drugcur_p
#       + base_year_hba1c_model_cat
#       + riskfactor_model
#       + birthsex
#       + site
#       + alcriskhi
#       + ns(depsum_i, 3)
#       + (1|studyid)
#       , data = as.data.frame(mod_rawmodeldata)
#       , family = binomial
#       , control = glmerControl(optimizer = "bobyqa")
# )
# car::Anova(fitlist[["fit_glmer"]])
# 
# car::Anova(fitlist[["fit_glmer"]], type = "III")
# 
# fitlist[["fit_glmer_inter"]] = glmer(dm_control ~ ns(depsum_i, 3)
#                   + ns(base_age_model, 3)
#                   + ns(time_since_dm_diagnosis, 3)
#                   + ns(sqrt_cd4_num, 3)
#                   + ns(bmi, 3)
#                   + race_hispanic_model
#                   + hiv1rna_detected
#                   + drugcur_p
#                   + base_year_hba1c_model_cat
#                   + riskfactor_model
#                   + birthsex
#                   + site
#                   + alcriskhi
#                   + ns(depsum_i, 3) *race_hispanic_model
#                   + (1|studyid)
#                   , data = as.data.frame(mod_rawmodeldata)
#                   , family = binomial
#                   , control = glmerControl(optimizer = "bobyqa")
# )
# car::Anova(fitlist[["fit_glmer_inter"]])
# 
# car::Anova(fitlist[["fit_glmer_inter"]], type = "III")
# saveRDS(fitlist, "/home/xxxx/OneDrive/CNICS205_20220912/manfitlist.rds")


## -------------------------------------------------------------------------------------------------------
fitlist = readRDS("/home/xxxx/OneDrive/CNICS205_20220912/manfitlist.rds")

griddata = expand.grid('depsum_i' = sort(unique(mod_rawmodeldata$depsum_i, na.rm = TRUE))
) %>%
  as.data.frame()

fit = fitlist[["fit_lme_no_inter"]]
fit_lme_no_inter = fitlist[["fit_lme_no_inter"]]

fakedata = mod_rawmodeldata[rep(1,nrow(griddata)),] %>%
  dplyr::select(-depsum_i) %>%
  mutate(depsum_i = griddata$depsum_i) 

fakedata0 <- fakedata %>%
  mutate(depsum_i = 0)
fakedata_matrix <- model.matrix(terms(fit), data = fakedata)
fakedata_matrix0 <- model.matrix(terms(fit), data = fakedata0)

diff_matrix = fakedata_matrix - fakedata_matrix0

TV = vcov(fit)
Q= fixef(fit)
se_pred_grouping <- sqrt(diag((diff_matrix) %*% TV %*% t(diff_matrix)))
pred_grouping <- diff_matrix %*% Q

## Saving the results
method_result_change <- data.frame( "CONTVAR" = griddata$depsum_i
                                    , "point" = pred_grouping
                                    , "SE_Estimate" = se_pred_grouping
                                    , "lwr" =  pred_grouping - qnorm(0.975)*se_pred_grouping
                                    , "upr" = pred_grouping + qnorm(0.975)*se_pred_grouping
)

method_result_change %>%
  dplyr::filter(CONTVAR %in% c(10,20))

### Abstract fig B
FIGBCOL = "#66CCCC"

FIGBCOL = "black"

CONTVAR_results = method_result_change
YLIMITS = range(c(CONTVAR_results$lwr, CONTVAR_results$upr))
YLIMITS = c(-0.1,0.6)
YLAB = ("Adjusted Mean Difference in HbA1c Levels \n with 95% Confidence Interval")
XLAB = "Depression Symptoms"

p <- CONTVAR_results %>%
  ggplot( aes(x=CONTVAR, y= point)) +
  #geom_line( size =1, col = FIGBCOL)+ 
  ylab(YLAB) +
  xlab(XLAB)  +
  # scale_x_continuous(breaks = sort(as.numeric(unique(modeldata$CONTVAR_20160101))), limits = c(2, 80) ) +
  scale_y_continuous(limits = YLIMITS) +
  # theme_classic() +
  geom_hline(yintercept = 0, colour = FIGBCOL) +
  #  geom_vline(xintercept = 0, linetype = "dashed", size = 0.8, color = FIGBCOL) +
  scale_linetype_manual(values=c("solid","twodash", "dashed"))+
  geom_ribbon(aes(ymin = lwr, ymax = upr), show.legend=FALSE, alpha=0.20) +
#  geom_errorbar(aes(ymin= lwr, ymax= upr), colour=FIGBCOL, width=.2, position=pd) +
  geom_line(size =1) 

p1 = p  +
  annotate(geom = 'text'
           , label = paste0(" overall p = ", round(car::Anova(fit)["ns(depsum_i, 3)", "Pr(>Chisq)"], 3))
           , x = -Inf
           , y = YLIMITS[2]
           , hjust = 0
           , vjust = 1
           , size = 5) +
  theme_bw()+
  theme(
    axis.title.x = element_text( size = 14),
    axis.title.y = element_text( size = 14),
    axis.text.y = element_text( size = 14),
    axis.text.x = element_text( size = 14),
    legend.position="top"
    , legend.title = element_text(colour="white", size =12)
    , legend.text = element_text(size = 12)
  ) 

abstract_fig_b = p1 + ggtitle("B. HbA1c Levels")


## -------------------------------------------------------------------------------------------------------
fit = fitlist[["fit_glmer"]]

fakedata = mod_rawmodeldata[rep(1,nrow(griddata)),] %>%
  dplyr::select(-depsum_i) %>%
  mutate(depsum_i = griddata$depsum_i) 

fakedata0 <- fakedata %>%
  mutate(depsum_i = 0)
fakedata_matrix <- model.matrix(terms(fit), data = fakedata)
fakedata_matrix0 <- model.matrix(terms(fit), data = fakedata0)

diff_matrix = fakedata_matrix - fakedata_matrix0

TV = vcov(fit)
Q= fixef(fit)
se_pred_grouping <- sqrt(diag((diff_matrix) %*% TV %*% t(diff_matrix)))
pred_grouping <- diff_matrix %*% Q

## Saving the results
CONTVAR_results <- data.frame("CONTVAR" = griddata$depsum_i
                                   , "point_exp" = exp(pred_grouping)
                                   , "SE_Estimate" = se_pred_grouping
                                   , "lwr_exp" =  exp(pred_grouping - qnorm(0.975)*se_pred_grouping)
                                   , "upr_exp" = exp(pred_grouping + qnorm(0.975)*se_pred_grouping)
) 

CONTVAR_results %>%
  dplyr::filter(CONTVAR %in% c(10,20))

FIGACOL = "#FF6699"
FIGACOL = "black"

### Abstract fig A
YLIMITS = range(c(CONTVAR_results$lwr_exp, CONTVAR_results$upr_exp))
YLIMITS = c(0.9, 4)
YLAB = paste0("Adjusted Odds Ratio of Uncontrolled Diabetes \n with 95% Confidence Interval")

XLAB = "Depression Symptoms"

p <- CONTVAR_results %>%
  ggplot( aes(x=CONTVAR, y= point_exp)) +
  #geom_line( size =1, col = FIGACOL)+ 
  # ylab(str_wrap(YLAB, width = 40)) +
  ylab(YLAB) +
  xlab(XLAB)  +
  # scale_x_continuous(breaks = sort(as.numeric(unique(modeldata$CONTVAR_20160101))), limits = c(2, 80) ) +
  scale_y_continuous(limits = YLIMITS, trans = "log", breaks = c(1/32, 1/16,0.125,0.25, 0.5, 1, 2, 4,8,16,32)) +
  # theme_classic() +
  geom_hline(yintercept = 1, colour = FIGACOL) +
  # geom_vline(xintercept = REF, linetype = "dashed", size = 0.8, color = FIGACOL) +
  scale_linetype_manual(values=c("solid","twodash", "dashed"))+
  geom_ribbon(aes(ymin = lwr_exp, ymax = upr_exp), show.legend=FALSE, alpha=0.20) +
  #  geom_vline(xintercept = REF, colour = FIGACOL) +
 # geom_errorbar(aes(ymin= lwr_exp, ymax= upr_exp), colour=FIGACOL, width=.2, position=pd) +
  geom_line(size =1) 

p1 = p  +
  #  scale_x_continuous( breaks = as.vector(c(Q1Q2,Q3))) +
   annotate(geom = 'text'
            , label = paste0(" overall p = ", round(car::Anova(fit)["ns(depsum_i, 3)", "Pr(>Chisq)"], 2))
            , x = -Inf
            , y = YLIMITS[2]
            , hjust = 0
            , vjust = 1
            , size = 5) +
  theme_bw()+
  theme(
    axis.title.x = element_text( size = 14),
    axis.title.y = element_text( size = 14),
    axis.text.y = element_text( size = 14),
    axis.text.x = element_text( size = 14),
    legend.position="top"
    , legend.title = element_text(colour="white", size =20)
    , legend.text = element_text(size = 16)
  ) 
#+ labs(title = CAPTION)

abstract_fig_a = p1 + ggtitle("A. HbA1c >7.0% vs. HbA1c \u2264 7.0%")

ggarrange(abstract_fig_a, abstract_fig_b
          , nrow = 1
          , ncol = 2
)

ggsave(paste0("/home/xxxx/OneDrive/CNICS205_20220912/writing/man_figure3_", Sys.Date(), ".jpg"), width = 11, height = 6)


## -------------------------------------------------------------------------------------------------------
griddata = expand.grid('depsum_i' = sort(unique(mod_rawmodeldata$depsum_i, na.rm = TRUE))
                       ,  'race_hispanic_model'= levels(mod_rawmodeldata$race_hispanic_model)
) %>%
  as.data.frame()
fit = fitlist[["fit_lme_inter"]]

fakedata = mod_rawmodeldata[rep(1,nrow(griddata)),] %>%
  dplyr::select(-depsum_i, race_hispanic_model) %>%
  mutate(depsum_i = griddata$depsum_i) %>%
  mutate(race_hispanic_model = griddata$race_hispanic_model)

fakedata0 <- fakedata %>%
  mutate(depsum_i = 0)
fakedata_matrix <- model.matrix(terms(fit), data = fakedata)
fakedata_matrix0 <- model.matrix(terms(fit), data = fakedata0)

diff_matrix = fakedata_matrix - fakedata_matrix0

TV = vcov(fit)
Q= fixef(fit)
se_pred_grouping <- sqrt(diag((diff_matrix) %*% TV %*% t(diff_matrix)))
pred_grouping <- diff_matrix %*% Q

## Saving the results
CONTVAR_results <- data.frame("race_hispanic_model" = griddata$race_hispanic_model
                              , "CONTVAR" = griddata$depsum_i
                              , "point" = pred_grouping
                              , "SE_Estimate" = se_pred_grouping
                              , "lwr" =  pred_grouping - qnorm(0.975)*se_pred_grouping
                              , "upr" = pred_grouping + qnorm(0.975)*se_pred_grouping
)

CONTVAR_results %>%
  dplyr::filter(CONTVAR %in% c(20))%>%
  mutate(est95ci = EST95CI(point, lwr, upr)) 

# Set dodge position
pd <- position_dodge(0.1) 


# Your data
YLIMITS <- range(c(CONTVAR_results$lwr, CONTVAR_results$upr)) + c(0,0.2)
YLAB <- paste0("Adjusted Mean Difference in ", "HbA1c Levels with 95% Confidence Interval Across Race and Ethnicity")
CAPTION <- ""


XLAB <- "Depression Symptoms"

OVERALLP= round(car::Anova(fit)["ns(depsum_i, 3):race_hispanic_model", "Pr(>Chisq)"], 3)

# Plot
p0 <- CONTVAR_results %>%
  ggplot(aes(x = CONTVAR, y = point, group = race_hispanic_model)) +
  geom_line(size = 1) + 
  ylab(str_wrap(YLAB, width = 50)) +
  xlab(XLAB) +
  scale_y_continuous(limits = YLIMITS, breaks = c(-4:2) * 1 ) +
  # geom_vline(xintercept = 0, colour = "gray") +
  geom_hline(yintercept = 0, colour = "black") +
  scale_linetype_manual(values = c("solid", "twodash", "dashed")) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), show.legend=FALSE, alpha=0.20) +
  # geom_errorbar(aes(ymin = lwr, ymax = upr, group = race_hispanic_model),
  #               width = 0.2) +
  geom_line() +
  # geom_point(size = 1) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "none",
    legend.title = element_text(colour = "white", size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "transparent"),
    plot.title = element_text(size = 16)
  ) +
  facet_wrap(~race_hispanic_model) +
  labs(title = CAPTION) +
  theme(strip.background = element_rect(fill = "white", color = "gray90")) +
  theme(strip.text = element_text(colour = 'black', size = 16)) +
  theme(panel.grid.major = element_blank()) +
 annotate("text", x = max(CONTVAR_results$CONTVAR), y = YLIMITS[2] - 0.1, 
          label = paste0("Interaction p = ", OVERALLP), hjust = 1, size = 4)

p0+ ggtitle("HbA1c Levels")

ggsave(paste0("/home/xxxx/OneDrive/CNICS205_20220912/writing/man_figure4_", Sys.Date(), ".jpg"), width = 11, height = 11)



griddata = expand.grid('depsum_i' = sort(unique(mod_rawmodeldata$depsum_i, na.rm = TRUE))
                       ,  'race_hispanic_model'= levels(mod_rawmodeldata$race_hispanic_model)
) %>%
  as.data.frame()


## -------------------------------------------------------------------------------------------------------
fit = fitlist[["fit_glmer_inter"]]
OVERALLP= round(car::Anova(fit)["ns(depsum_i, 3):race_hispanic_model", "Pr(>Chisq)"], 2)
fakedata = mod_rawmodeldata[rep(1,nrow(griddata)),] %>%
  dplyr::select(-depsum_i, race_hispanic_model) %>%
  mutate(depsum_i = griddata$depsum_i) %>%
  mutate(race_hispanic_model = griddata$race_hispanic_model)

fakedata0 <- fakedata %>%
  mutate(depsum_i = 0)
fakedata_matrix <- model.matrix(terms(fit), data = fakedata)
fakedata_matrix0 <- model.matrix(terms(fit), data = fakedata0)

diff_matrix = fakedata_matrix - fakedata_matrix0

TV = vcov(fit)
Q= fixef(fit)
se_pred_grouping <- sqrt(diag((diff_matrix) %*% TV %*% t(diff_matrix)))
pred_grouping <- diff_matrix %*% Q

## Saving the results
CONTVAR_results <- data.frame("race_hispanic_model" = griddata$race_hispanic_model
                              , "CONTVAR" = griddata$depsum_i
                              , "point_exp" = exp(pred_grouping)
                              , "SE_Estimate" = se_pred_grouping
                              , "lwr_exp" =  exp(pred_grouping - qnorm(0.975)*se_pred_grouping)
                              , "upr_exp" = exp(pred_grouping + qnorm(0.975)*se_pred_grouping)
)


  

YLIMITS <- range(c(CONTVAR_results$lwr_exp, CONTVAR_results$upr_exp)) + c(0, 0.2)
YLAB <- "Adjusted Odds Ratio for the Lack of Diabetes Control with 95% Confidence Interval\nAcross Race and Ethnicity"
CAPTION <- ""
REF <- as.numeric(CONTVAR_results$REF[1])
XLAB <- Hmisc::label(rawmodeldata[, "depsum_i"])
XLAB <- "Depression Symptoms"

# Plot
p0 <- CONTVAR_results %>%
  ggplot(aes(x = CONTVAR, y = point_exp, group = race_hispanic_model)) +
  geom_line(size = 1) +
  ylab(YLAB) +
  xlab(XLAB) +
  scale_y_continuous(limits = YLIMITS, breaks = c(0.05, 0.10, 0.50,0.25, 1, 2, 4,8), trans = "log") +
  #  geom_vline(xintercept = 0, colour = "black") +
  geom_ribbon(aes(ymin = lwr_exp, ymax = upr_exp), show.legend=FALSE, alpha=0.20) +
  geom_hline(yintercept = 1, colour = "black") +
  scale_linetype_manual(values = c("solid", "twodash", "dashed")) +
  # geom_errorbar(aes(ymin = lwr_exp, ymax = upr_exp, group = race_hispanic_model),
  #               width = 0.2) +
  geom_line() +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "none",
    legend.title = element_text(colour = "white", size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "transparent"),
    plot.title = element_text(size = 16)
  ) +
  facet_wrap(~race_hispanic_model) +
  theme(strip.background = element_rect(fill = "white", color = "gray90")) +
  theme(strip.text = element_text(colour = 'black', size = 16)) +
  theme(panel.grid.major = element_blank()) +
  annotate("text", x = max(CONTVAR_results$CONTVAR) - 2, y = YLIMITS[2] - 0.1, 
         label = paste0("Interaction p = ", OVERALLP), hjust = 1, size = 5)

# Display plot
pf= p0+ ggtitle("HbA1c >7.0% vs. HbA1c \u2264 7.0%")
ggsave(paste0("/home/xxxx/OneDrive/CNICS205_20220912/writing/man_figure5_", Sys.Date(), ".jpg")
       , p0+ ggtitle("HbA1c >7.0% vs. HbA1c \u2264 7.0%"), width = 11, height = 11)



## -------------------------------------------------------------------------------------------------------
modeldata_sup = tmp3 %>%
  # dplyr::filter(hba1c_date > dm_diagnosis_date) %>%
  # dplyr::filter(time_since_cnics_entry >= 0) %>%
  # dplyr::filter(!is.na(depsum_i) | !is.na(depsum_cc)) %>%
  # dplyr::filter(!is.na(hba1c)) %>%
  # dplyr::filter(!is.na(cd4_num)) %>%
  # dplyr::filter(!is.na(bmi)) %>%
  # dplyr::filter(!is.na(alcriskhi)) %>%
  # dplyr::filter(!is.na(hiv1rna_detected)) %>%
  # dplyr::filter(!is.na(drugcur_p) | !is.na(drugever_p))  %>%
  group_by(studyid) %>%
  arrange(studyid, hba1c_date) %>%
  mutate(time_since_first_hba1c = as.numeric(hba1c_date - min(hba1c_date)))%>% 
  mutate_at(cols_as_factors, factor) %>%
  mutate(med_dep = factor(med_dep, levels=c(0,1), labels=c("No", "Yes")))%>%
  mutate(med_dm = factor(med_dm, levels=c(0,1), labels=c("No", "Yes")))%>%
  mutate(med_hiv = factor(med_hiv, levels=c(0,1), labels=c("No", "Yes")))%>%
  mutate(depression = factor(depression, levels=c(0,1), labels=c("No", "Yes")))%>%
  mutate(dm_control_grp = factor(dm_control, levels=c(0,1), labels=c("<=7%", ">7%")))%>%
  mutate(year_depsum_i_cat = cut(year_depsum_i
                                 , breaks = c(2005,2010,2016,2023)
                                 , include.lowest = TRUE
  )) %>%
  mutate(year_depsum_cc_cat = cut(year_depsum_cc
                                  , breaks = c(2005,2010,2016,2023)
                                  , include.lowest = TRUE
  )) %>%
  mutate(year_hba1c_cat = cut(year_hba1c
                              , breaks = c(2005,2010,2016,2023)
                              , include.lowest = TRUE
  )) %>%
  mutate(hiv1rna_detected = factor(hiv1rna_detected2))%>%
  mutate(diff_hba1c_phq = as.numeric(hba1c_date - depsum_i_date) ) 



baseline_dm_sup <- modeldata_sup %>%
  arrange(studyid, age_updating) %>%
  group_by(studyid) %>%
  mutate(followup_duration = max(time_since_dm_diagnosis) - min(time_since_dm_diagnosis)) %>%
  mutate(number_of_outcome_measures = n()) %>%
  slice(1) %>%
  mutate(included = if_else(studyid %in% baseline_dm$studyid, "Included", "Excluded"))



table1vars_man <- c(#"age_updating"
                     "birthsex"
                     , "race_hispanic_model"
                    , "transgender_model"
                    , "riskfactor_model" 
                   # , "med_hiv"
                  #  , "hiv1rna_detected"
                  #  , "time_since_dm_diagnosis"
                  #  , "cd4_num" 
                  #  , "bmi"
                  #  , "hba1c"
                  #  , "dm_control_grp"
                  #  , "med_dm"
                   # , "year_hba1c_cat"
                  #  , "depsum_i_cat"
                  #  , "depression"
                  #  , "med_dep"
                    , "hbv_status"
                    , "hcv_status"
                  #  , "alcriskhi"
                  #  , "time_since_cnics_entry"
                  #  , "followup_duration"
                  #  , "diff_hba1c_phq"
                    )

tbl1_formula <- as.formula(paste(paste(table1vars_man, collapse = "+") ,  "~", "included"))

tbl1 <- summaryM(
  tbl1_formula
  , data = baseline_dm_sup
  , quant = c(0.25,0.5,0.75)
  , continuous = 8
  , overall = OVERALL
  , na.include = TRUE
  , test = TRUE
)

html(tbl1
     , caption = ""
     , vnames = "label"
     , exclude1 = FALSE
     , digits= 5
     , size='normal'
     , prmsd = FALSE
     , middle.bold=TRUE
     , long= TRUE
     , rowsep = TRUE)



tbl1 <- summaryM_to_df(
  tbl1
  , long = TRUE
  , exclude1 = FALSE
  #  , prob = c(0.1,0.5,0.9)
  , what = "%"
  , digits = 2
  , prn = TRUE
  , vnames = "labels"
)

for (i in 3: (dim(tbl1)[2]-0)){
  tbl1[,i] <- slice_to_box(tbl1[,i])
  tbl1[,i] <- slice_to_box2(tbl1[,i])
  
}
names(tbl1)[1] <- "Variables"
#kable(tbl1)
kable(tbl1, "pandoc") %>%
  column_spec(1, width = "10em")

write.csv(tbl1, "/home/xxxx/OneDrive/CNICS205_20220912/writing/table1_sup.csv")

