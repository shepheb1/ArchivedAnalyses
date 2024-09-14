########## Project: Code for "Understanding Difference-in-Differences methods to evaluate policy effects with staggered adoption: an application to Medicaid and HIV"
########## Author: Julia C. Thome
########## Date: 02/05/24

library(lubridate)
library(Hmisc)
library(svMisc)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(geepack)
library(splines)
library(sjPlot)
library(did)
library(tidyverse)
library(table1)
library(NCmisc)
library(car)
library(latex2exp)
library(rms)

### Read in the data
sim.data <- read.csv("tutorial_data.csv")

# data cleaning
sim.data$Clinical_Retention <- as.numeric(as.character(sim.data$Clinical_Retention))

##### Section 3.1 #####

## subsetting to only 2013 and 2014 data for these models
sim.data.diff <- sim.data %>% subset(Year %in% c(2013,2014))
sim.data.diff$Year <- as.factor(sim.data.diff$Year)
sim.data.diff$Year <- droplevels(sim.data.diff$Year)
sim.data.diff$Year <- relevel(sim.data.diff$Year,ref = "2013")
sim.data.diff$expanded_2014 <- ifelse(sim.data.diff$Expand_Year %in% 2014, 1,0)

## subsetting to people in which clinical retention can be defined 
## details on clinical retention definition:
#Definition: "Retention will be defined as adherence to the Institute of Medicine/National HIV/AIDS Strategy indicator between entry to the study population (first clinic encounter during the study period) and December 31 of the closing calendar year. Adherence to the definition will be assessed in each calendar year: >= 2 HIV primary care visits within the 12-month period, >90 days apart."
#- Defined for patients with encounters after 2011
#- Adherence to the definition was assessed in each calendar year after the year of entry into the study.
#- For patients included in the study with at least one primary care visits:
#   - Clinical retention will be 1 for a year if they had at least 2 primary care visits in that year at least 90 days apart
#   - Clinical retention will be 0 for a year if they had encounters before study period began, are alive, but don't have encounter data for a specific year up to that point in study period
#   - Clinical retention will be 0 for a year if they have had first encounters before year being considered, are alive, but don't have encounter data for that year
#   - Clinical retention will be 0 for a year if they have had first encounters before year being considered, died after year being considered, but don't have encounter data for year being considered
#   - Clinical Retention will not be considered after the last encounter date available during the study period
#   - Clinical Retention will not be considered if death is within 90 days of Jan 1st of death year 

sim.data.diff.sub <- sim.data.diff %>% subset(retention_denom %in% 1) %>% mutate(Year = as.factor(Year))


## estimate of ATT(2014,2014)
att_2014_2014_3.1 <- ((mean(sim.data.diff.sub$Clinical_Retention[(sim.data.diff.sub$Year == 2014) & (sim.data.diff.sub$expanded_2014 == 1)],na.rm = T) 
                       - mean(sim.data.diff.sub$Clinical_Retention[(sim.data.diff.sub$Year == 2013) & (sim.data.diff.sub$expanded_2014 == 1)],na.rm = T) ) 
                      - (mean(sim.data.diff.sub$Clinical_Retention[(sim.data.diff.sub$Year == 2014) & (sim.data.diff.sub$expanded_2014 == 0)],na.rm = T) 
                         - mean(sim.data.diff.sub$Clinical_Retention[(sim.data.diff.sub$Year == 2013) & (sim.data.diff.sub$expanded_2014 == 0)],na.rm = T) ) ) %>% round(3)

## bootstrapping confidence intervals
nboot <- 1000
unique_ids <- unique(sim.data.diff.sub$NAID)
att_2014_2014_3.1_boot_list <- NULL
for(i in 1:nboot){
  # creating bootstrap sample by sampling IDs with replacement and grabbing all data from sim.data.diff.sub for those IDs
  boot.id.table <- table(sort(sample(unique_ids, size = length(unique_ids), replace = TRUE)))
  boot.id.data <- data.frame(NAID = names(boot.id.table), n = as.numeric(boot.id.table))
  merged_data <- merge(boot.id.data, sim.data.diff.sub, by = "NAID")
  dboot <- merged_data %>% uncount(n) 
  
  # estimate of ATT(2014,2014) for this bootstrap iteration
  att_2014_2014_3.1_boot <- ((mean(dboot$Clinical_Retention[(dboot$Year == 2014) & (dboot$expanded_2014 == 1)],na.rm = T) - mean(dboot$Clinical_Retention[(dboot$Year == 2013) & (dboot$expanded_2014 == 1)],na.rm = T) ) - (mean(dboot$Clinical_Retention[(dboot$Year == 2014) & (dboot$expanded_2014 == 0)],na.rm = T) - mean(dboot$Clinical_Retention[(dboot$Year == 2013) & (dboot$expanded_2014 == 0)],na.rm = T) ) ) 
  att_2014_2014_3.1_boot_list <- c(att_2014_2014_3.1_boot_list, att_2014_2014_3.1_boot)
}

## 95% bootstrapped CI for ATT(2014,2014)
att_2014_2014_3.1_lb <- quantile(att_2014_2014_3.1_boot_list,.025) %>% round(3) %>% as.numeric()
att_2014_2014_3.1_ub <- quantile(att_2014_2014_3.1_boot_list,.975) %>% round(3) %>% as.numeric()


##### Section 3.2 #####
## fitting model (1)
mod_att_2014_2014_3.2 <- ols(Clinical_Retention ~ expanded_2014*Year, data = sim.data.diff.sub,x = TRUE, y = TRUE)

## estimate of ATT(2014,2014)
att_2014_2014_3.2 <- mod_att_2014_2014_3.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014) that accounts for clustering
robust_mod_att_2014_2014_3.2 <- robcov(fit = mod_att_2014_2014_3.2, cluster = sim.data.diff.sub$NAID)
robust_se <- vcov(robust_mod_att_2014_2014_3.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014)
att_2014_2014_3.2_lb <- (att_2014_2014_3.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_3.2_ub <- (att_2014_2014_3.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


##### Section 4.1 #####
## data prep for did package
sim.data$Year <- as.numeric(as.character(sim.data$Year))
sim.data.diff.clin.ret <- sim.data %>% subset(Year %in% 2012:2017) %>% subset(retention_denom %in% 1)
sim.data.diff.clin.ret$expanded_2014 <- ifelse(sim.data.diff.clin.ret$Expand_Year == 2014, 1,0)
sim.data.diff.clin.ret$expanded_ever <- as.factor(sim.data.diff.clin.ret$expanded_ever)
sim.data.diff.clin.ret$first_treat <- ifelse(sim.data.diff.clin.ret$Expand_Year %in% 2012:2017, sim.data.diff.clin.ret$Expand_Year, 0)
sim.data.diff.clin.ret2 <- sim.data.diff.clin.ret %>% dplyr::select(-vl_suppress) %>% unique() 
sim.data.diff.clin.ret2$NAID <- as.numeric(sim.data.diff.clin.ret2$NAID)


## first_treat as ITT because of movers (0.62% of patients have itt != actual)
sim.data.diff.clin.ret3 <- sim.data.diff.clin.ret2 %>% group_by(NAID) %>% mutate(n = length(unique(state)), mover = as.numeric(n > 1), first_state = state[which.min(Year)], first_treat_itt = first_treat[which.min(Year)]) %>% ungroup() %>% as.data.frame() 


## DID approach adjusting for race only
mod_att_2014_2014_race_4.1 <- att_gt(yname = "Clinical_Retention", 
                                     gname = "first_treat_itt",
                                     xformla = ~ RACE_black, 
                                     clustervars = "NAID",
                                     tname = "Year", 
                                     data = sim.data.diff.clin.ret3 %>% subset(Year %in% 2013:2014), 
                                     est_method = "reg",
                                     panel = T,
                                     allow_unbalanced_panel = T,
                                     bstrap = T,
                                     control_group="notyettreated",
                                     idname = "NAID")

## estimate and 95% CI for ATT(2014,2014) adjusting for race only
att_2014_2014_race_4.1 <- mod_att_2014_2014_race_4.1$att %>% round(3)
att_2014_2014_race_4.1_lb <- (att_2014_2014_race_4.1 - 1.96*mod_att_2014_2014_race_4.1$se) %>% round(3)
att_2014_2014_race_4.1_ub <- (att_2014_2014_race_4.1 + 1.96*mod_att_2014_2014_race_4.1$se) %>% round(3)

## DID approach adjusting for region only
mod_att_2014_2014_region_4.1 <- att_gt(yname = "Clinical_Retention", 
                                       gname = "first_treat_itt",
                                       xformla = ~ region, 
                                       clustervars = "NAID",
                                       tname = "Year", 
                                       data = sim.data.diff.clin.ret3 %>% subset(Year %in% 2013:2014), 
                                       est_method = "reg",
                                       panel = T,
                                       allow_unbalanced_panel = T,
                                       bstrap = T,
                                       control_group="notyettreated",
                                       idname = "NAID")

## estimate and 95% CI for ATT(2014,2014) adjusting for region only
att_2014_2014_region_4.1 <- mod_att_2014_2014_region_4.1$att %>% round(3)
att_2014_2014_region_4.1_lb <- (att_2014_2014_region_4.1 - 1.96*mod_att_2014_2014_region_4.1$se) %>% round(3)
att_2014_2014_region_4.1_ub <- (att_2014_2014_region_4.1 + 1.96*mod_att_2014_2014_region_4.1$se) %>% round(3)

## DID approach to estimate ATT(2014,2014) adjusting for age, race, sex, and region 
mod_att_2014_2014_all_covariates_4.1 <- att_gt(yname = "Clinical_Retention", 
                                               gname = "first_treat_itt",
                                               xformla = ~ age_first_in_study_cent_10 + RACE_black + SEX + region, 
                                               clustervars = "NAID",
                                               tname = "Year", 
                                               data = sim.data.diff.clin.ret3 %>% subset(Year %in% 2013:2014), 
                                               est_method = "reg",
                                               panel = T,
                                               allow_unbalanced_panel = T,
                                               bstrap = T,
                                               control_group="notyettreated",
                                               idname = "NAID")

## estimate and 95% CI for ATT(2014,2014) adjusting for age, race, sex, and region 
att_2014_2014_all_covariates_4.1 <- mod_att_2014_2014_all_covariates_4.1$att %>% round(3)
att_2014_2014_all_covariates_4.1_lb <- (att_2014_2014_all_covariates_4.1 - 1.96*mod_att_2014_2014_all_covariates_4.1$se) %>% round(3)
att_2014_2014_all_covariates_4.1_ub <- (att_2014_2014_all_covariates_4.1 + 1.96*mod_att_2014_2014_all_covariates_4.1$se) %>% round(3)



##### Section 4.2 #####
### Race
## fitting model (2)
mod_att_2014_2014_race_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year + RACE_black, data = sim.data.diff.sub,x = TRUE, y = TRUE)

## estimate of ATT(2014,2014)
att_2014_2014_race_4.2 <- mod_att_2014_2014_race_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014) that accounts for clustering
robust_mod_att_2014_2014_race_4.2 <- robcov(fit = mod_att_2014_2014_race_4.2, cluster = sim.data.diff.sub$NAID)
robust_se <- vcov(robust_mod_att_2014_2014_race_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014)
att_2014_2014_race_4.2_lb <- (att_2014_2014_race_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_race_4.2_ub <- (att_2014_2014_race_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


### Race-by-time
## fitting model (3)
mod_att_2014_2014_race_by_time_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year + RACE_black*Year, data = sim.data.diff.sub,x = TRUE, y = TRUE)

## estimate of ATT(2014,2014)
att_2014_2014_race_by_time_4.2 <- mod_att_2014_2014_race_by_time_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014) that accounts for clustering
robust_mod_att_2014_2014_race_by_time_4.2 <- robcov(fit = mod_att_2014_2014_race_by_time_4.2, cluster = sim.data.diff.sub$NAID)
robust_se <- vcov(robust_mod_att_2014_2014_race_by_time_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014)
att_2014_2014_race_by_time_4.2_lb <- (att_2014_2014_race_by_time_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_race_by_time_4.2_ub <- (att_2014_2014_race_by_time_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


### One model for each race group

## fitting model (1) among race = not black
mod_att_2014_2014_not_black_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year, data = sim.data.diff.sub %>% subset(RACE_black %in% "Not Black" ),x = TRUE, y = TRUE)

## estimate of ATT(2014,2014|race = not black)
att_2014_2014_not_black_4.2 <- mod_att_2014_2014_not_black_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014|race = not black) that accounts for clustering
robust_mod_att_2014_2014_not_black_4.2 <- robcov(fit = mod_att_2014_2014_not_black_4.2, cluster = sim.data.diff.sub$NAID[sim.data.diff.sub$RACE_black %in% "Not Black"])
robust_se <- vcov(robust_mod_att_2014_2014_not_black_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014|race = not black)
att_2014_2014_not_black_4.2_lb <- (att_2014_2014_not_black_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_not_black_4.2_ub <- (att_2014_2014_not_black_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()

## fitting model (1) among race = black
mod_att_2014_2014_black_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year, data = sim.data.diff.sub %>% subset(RACE_black %in% "Black" ),x = TRUE, y = TRUE)

## estimate of ATT(2014,2014|race = black)
att_2014_2014_black_4.2 <- mod_att_2014_2014_black_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014|race = black) that accounts for clustering
robust_mod_att_2014_2014_black_4.2 <- robcov(fit = mod_att_2014_2014_black_4.2, cluster = sim.data.diff.sub$NAID[sim.data.diff.sub$RACE_black %in% "Black"])
robust_se <- vcov(robust_mod_att_2014_2014_black_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014|race = black)
att_2014_2014_black_4.2_lb <- (att_2014_2014_black_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_black_4.2_ub <- (att_2014_2014_black_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()



### Region
## fitting model adjusting for region only
mod_att_2014_2014_region_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year + region, data = sim.data.diff.sub,x = TRUE, y = TRUE)

## estimate of ATT(2014,2014)
att_2014_2014_region_4.2 <- mod_att_2014_2014_region_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014) that accounts for clustering
robust_mod_att_2014_2014_region_4.2 <- robcov(fit = mod_att_2014_2014_region_4.2, cluster = sim.data.diff.sub$NAID)
robust_se <- vcov(robust_mod_att_2014_2014_region_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014)
att_2014_2014_region_4.2_lb <- (att_2014_2014_region_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_region_4.2_ub <- (att_2014_2014_region_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


### Region-by-Time
## fitting model adjusting for region-by-time only
mod_att_2014_2014_region_by_time_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year + region*Year, data = sim.data.diff.sub,x = TRUE, y = TRUE)

## estimate of ATT(2014,2014)
att_2014_2014_region_by_time_4.2 <- mod_att_2014_2014_region_by_time_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014) that accounts for clustering
robust_mod_att_2014_2014_region_by_time_4.2 <- robcov(fit = mod_att_2014_2014_region_by_time_4.2, cluster = sim.data.diff.sub$NAID)
robust_se <- vcov(robust_mod_att_2014_2014_region_by_time_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014)
att_2014_2014_region_by_time_4.2_lb <- (att_2014_2014_region_by_time_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_region_by_time_4.2_ub <- (att_2014_2014_region_by_time_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()




### One model for each region group
## fitting model (1) among region = SE
mod_att_2014_2014_SE_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year, data = sim.data.diff.sub %>% subset(region %in% "SE" ),x = TRUE, y = TRUE)

## estimate of ATT(2014,2014|region = SE)
att_2014_2014_SE_4.2 <- mod_att_2014_2014_SE_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014|region = SE) that accounts for clustering
robust_mod_att_2014_2014_SE_4.2 <- robcov(fit = mod_att_2014_2014_SE_4.2, cluster = sim.data.diff.sub$NAID[sim.data.diff.sub$region %in% "SE"])
robust_se <- vcov(robust_mod_att_2014_2014_SE_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014|region = SE)
att_2014_2014_SE_4.2_lb <- (att_2014_2014_SE_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_SE_4.2_ub <- (att_2014_2014_SE_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


## fitting model (1) among region = WE
mod_att_2014_2014_WE_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year, data = sim.data.diff.sub %>% subset(region %in% "WE" ),x = TRUE, y = TRUE)

## estimate of ATT(2014,2014|region = WE)
att_2014_2014_WE_4.2 <- mod_att_2014_2014_WE_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014|region = WE) that accounts for clustering
robust_mod_att_2014_2014_WE_4.2 <- robcov(fit = mod_att_2014_2014_WE_4.2, cluster = sim.data.diff.sub$NAID[sim.data.diff.sub$region %in% "WE"])
robust_se <- vcov(robust_mod_att_2014_2014_WE_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014|region = WE)
att_2014_2014_WE_4.2_lb <- (att_2014_2014_WE_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_WE_4.2_ub <- (att_2014_2014_WE_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


## fitting model (1) among region = NE
mod_att_2014_2014_NE_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year, data = sim.data.diff.sub %>% subset(region %in% "NE" ),x = TRUE, y = TRUE)

## estimate of ATT(2014,2014|region = NE)
att_2014_2014_NE_4.2 <- mod_att_2014_2014_NE_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014|region = NE) that accounts for clustering
robust_mod_att_2014_2014_NE_4.2 <- robcov(fit = mod_att_2014_2014_NE_4.2, cluster = sim.data.diff.sub$NAID[sim.data.diff.sub$region %in% "NE"])
robust_se <- vcov(robust_mod_att_2014_2014_NE_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014|region = NE)
att_2014_2014_NE_4.2_lb <- (att_2014_2014_NE_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_NE_4.2_ub <- (att_2014_2014_NE_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


## fitting model (1) among region = MW
mod_att_2014_2014_MW_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year, data = sim.data.diff.sub %>% subset(region %in% "MW" ),x = TRUE, y = TRUE)

## estimate of ATT(2014,2014|region = MW)
att_2014_2014_MW_4.2 <- mod_att_2014_2014_MW_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014|region = MW) that accounts for clustering
robust_mod_att_2014_2014_MW_4.2 <- robcov(fit = mod_att_2014_2014_MW_4.2, cluster = sim.data.diff.sub$NAID[sim.data.diff.sub$region %in% "MW"])
robust_se <- vcov(robust_mod_att_2014_2014_MW_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014|region = MW)
att_2014_2014_MW_4.2_lb <- (att_2014_2014_MW_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_MW_4.2_ub <- (att_2014_2014_MW_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()



## fitting model adjusting for age, race, sex, and region
mod_att_2014_2014_all_covariates_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year + ns(age_first_in_study_cent_10,3) + RACE_black + SEX + region, data = sim.data.diff.sub,x = TRUE, y = TRUE)

## estimate of ATT(2014,2014)
att_2014_2014_all_covariates_4.2 <- mod_att_2014_2014_all_covariates_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014) that accounts for clustering
robust_mod_att_2014_2014_all_covariates_4.2 <- robcov(fit = mod_att_2014_2014_all_covariates_4.2, cluster = sim.data.diff.sub$NAID)
robust_se <- vcov(robust_mod_att_2014_2014_all_covariates_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014)
att_2014_2014_all_covariates_4.2_lb <- (att_2014_2014_all_covariates_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_all_covariates_4.2_ub <- (att_2014_2014_all_covariates_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


## fitting model adjusting for age-by-time, race-by-tim, sex-by-time, and region-by-time
mod_att_2014_2014_all_covariates_by_time_4.2 <- ols(Clinical_Retention ~ expanded_2014*Year + ns(age_first_in_study_cent_10,3)*Year + RACE_black*Year + SEX*Year + region*Year, data = sim.data.diff.sub,x = TRUE, y = TRUE)

## estimate of ATT(2014,2014)
att_2014_2014_all_covariates_by_time_4.2 <- mod_att_2014_2014_all_covariates_by_time_4.2$coefficients['expanded_2014 * Year=2014'] %>% round(3)%>% as.numeric()

## calculating robust standard errors for ATT(2014,2014) that accounts for clustering
robust_mod_att_2014_2014_all_covariates_by_time_4.2 <- robcov(fit = mod_att_2014_2014_all_covariates_by_time_4.2, cluster = sim.data.diff.sub$NAID)
robust_se <- vcov(robust_mod_att_2014_2014_all_covariates_by_time_4.2)["expanded_2014 * Year=2014","expanded_2014 * Year=2014"] %>% sqrt( )

## 95% CI for ATT(2014,2014)
att_2014_2014_all_covariates_by_time_4.2_lb <- (att_2014_2014_all_covariates_by_time_4.2 - 1.96*robust_se) %>% round(3) %>% as.numeric()
att_2014_2014_all_covariates_by_time_4.2_ub <- (att_2014_2014_all_covariates_by_time_4.2 + 1.96*robust_se) %>% round(3) %>% as.numeric()


##### Section 5.1 #####
## Figure 1 values
## DID approach to estimate all possible ATT(g,t)'s adjusting for age, race, sex, and region
did_mod_clin <- att_gt(yname = "Clinical_Retention", 
                       gname = "first_treat_itt",
                       xformla = ~ age_first_in_study_cent_10 + RACE_black + SEX + region, 
                       clustervars = "NAID",tname = "Year", 
                       data = sim.data.diff.clin.ret3, 
                       est_method = "reg",
                       panel = T,
                       allow_unbalanced_panel = T,
                       bstrap = T,
                       control_group="notyettreated",
                       idname = "NAID",cband = FALSE)

## DID approach to aggregate based on years since expansion (w) to calculate all possible ATT_g(w)'s adjusting for age, race, sex, and region
agg.es.clin <- aggte(did_mod_clin, type = "dynamic",cband = FALSE)


## estimate of ATT(2014,2014)
att_2014_2014_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2014 & did_mod_clin$t == 2014] %>% round(3)
## 95% CI for ATT(2014,2014)
att_2014_2014_all_covariates_5.1_lb <-  att_2014_2014_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2014 & did_mod_clin$t == 2014]    
att_2014_2014_all_covariates_5.1_ub <-  att_2014_2014_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2014 & did_mod_clin$t == 2014]   

## estimate of ATT(2014,2015)
att_2014_2015_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2014 & did_mod_clin$t == 2015] %>% round(3)
## 95% CI for ATT(2014,2015)
att_2014_2015_all_covariates_5.1_lb <-  att_2014_2015_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2014 & did_mod_clin$t == 2015]       
att_2014_2015_all_covariates_5.1_ub <-  att_2014_2015_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2014 & did_mod_clin$t == 2015]       

## estimate of ATT(2014,2016)
att_2014_2016_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2014 & did_mod_clin$t == 2016] %>% round(3)
## 95% CI for ATT(2014,2016)
att_2014_2016_all_covariates_5.1_lb <-   att_2014_2016_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2014 & did_mod_clin$t == 2016]         
att_2014_2016_all_covariates_5.1_ub <-   att_2014_2016_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2014 & did_mod_clin$t == 2016]   

## estimate of ATT(2014,2017)
att_2014_2017_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2014 & did_mod_clin$t == 2017] %>% round(3)
## 95% CI for ATT(2014,2017)
att_2014_2017_all_covariates_5.1_lb <-   att_2014_2017_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2014 & did_mod_clin$t == 2017]          
att_2014_2017_all_covariates_5.1_ub <-   att_2014_2017_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2014 & did_mod_clin$t == 2017]  

## estimate of ATT(2015,2015)
att_2015_2015_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2015 & did_mod_clin$t == 2015] %>% round(3)
## 95% CI for ATT(2015,2015)
att_2015_2015_all_covariates_5.1_lb <-  att_2015_2015_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2015 & did_mod_clin$t == 2015]               
att_2015_2015_all_covariates_5.1_ub <-  att_2015_2015_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2015 & did_mod_clin$t == 2015] 

## estimate of ATT(2015,2016)
att_2015_2016_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2015 & did_mod_clin$t == 2016] %>% round(3)
## 95% CI for ATT(2015,2016)
att_2015_2016_all_covariates_5.1_lb <- att_2015_2016_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2015 & did_mod_clin$t == 2016]  
att_2015_2016_all_covariates_5.1_ub <- att_2015_2016_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2015 & did_mod_clin$t == 2016]  

## estimate of ATT(2015,2017)
att_2015_2017_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2015 & did_mod_clin$t == 2017] %>% round(3)
## 95% CI for ATT(2015,2017)
att_2015_2017_all_covariates_5.1_lb <- att_2015_2017_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2015 & did_mod_clin$t == 2017]   
att_2015_2017_all_covariates_5.1_ub <- att_2015_2017_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2015 & did_mod_clin$t == 2017]      

## estimate of ATT(2016,2016)
att_2016_2016_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2016 & did_mod_clin$t == 2016] %>% round(3)
## 95% CI for ATT(2016,2016)
att_2016_2016_all_covariates_5.1_lb <-  att_2016_2016_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2016 & did_mod_clin$t == 2016]        
att_2016_2016_all_covariates_5.1_ub <-   att_2016_2016_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2016 & did_mod_clin$t == 2016]       

## estimate of ATT(2016,2017)
att_2016_2017_all_covariates_5.1 <- did_mod_clin$att[did_mod_clin$group == 2016 & did_mod_clin$t == 2017] %>% round(3)
## 95% CI for ATT(2016,2017)
att_2016_2017_all_covariates_5.1_lb <-  att_2016_2017_all_covariates_5.1 - 1.96*did_mod_clin$se[did_mod_clin$group == 2016 & did_mod_clin$t == 2017]     
att_2016_2017_all_covariates_5.1_ub <-  att_2016_2017_all_covariates_5.1 + 1.96*did_mod_clin$se[did_mod_clin$group == 2016 & did_mod_clin$t == 2017]       

## estimate of ATT_g(0)
att_w0_5.1 <-  agg.es.clin$att.egt[agg.es.clin$egt == 0]%>% round(3)
## 95% CI for ATT_g(0)
att_w0_5.1_lb <- att_w0_5.1 - 1.96*agg.es.clin$se.egt[agg.es.clin$egt == 0]
att_w0_5.1_ub <-  att_w0_5.1 + 1.96*agg.es.clin$se.egt[agg.es.clin$egt == 0]

## estimate of ATT_g(1)
att_w1_5.1 <- agg.es.clin$att.egt[agg.es.clin$egt == 1]%>% round(3)
## 95% CI for ATT_g(1)
att_w1_5.1_lb <-  att_w1_5.1 - 1.96*agg.es.clin$se.egt[agg.es.clin$egt == 1]       
att_w1_5.1_ub <-  att_w1_5.1 + 1.96*agg.es.clin$se.egt[agg.es.clin$egt == 1] 

## estimate of ATT_g(2)
att_w2_5.1 <- agg.es.clin$att.egt[agg.es.clin$egt == 2]%>% round(3)
## 95% CI for ATT_g(2)
att_w2_5.1_lb <-  att_w2_5.1 - 1.96*agg.es.clin$se.egt[agg.es.clin$egt == 2]       
att_w2_5.1_ub <-  att_w2_5.1 + 1.96*agg.es.clin$se.egt[agg.es.clin$egt == 2] 

## estimate of ATT_g(3)
att_w3_5.1 <- agg.es.clin$att.egt[agg.es.clin$egt == 3] %>% round(3)
## 95% CI for ATT_g(3)
att_w3_5.1_lb <-  att_w3_5.1 - 1.96*agg.es.clin$se.egt[agg.es.clin$egt == 3]         
att_w3_5.1_ub <-   att_w3_5.1 + 1.96*agg.es.clin$se.egt[agg.es.clin$egt == 3]        

# estimate of ATT_g,w
att_overall_5.1 <- agg.es.clin$overall.att%>% round(3)
## 95% CI for ATT_g,w
att_overall_5.1_lb <- att_overall_5.1 - 1.96*agg.es.clin$overall.se
att_overall_5.1_ub <- att_overall_5.1 + 1.96*agg.es.clin$overall.se

##### Section 5.2 #####
## subsetting data to include data for years 2013,2014, and 2015
sim.data.diff.sub2 <- sim.data %>% subset(Year %in% 2013:2015) %>% mutate(Year = as.factor(Year))

sim.data.diff.sub2$Expand_Year_factor <- "Not-yet Expanded"
sim.data.diff.sub2$Expand_Year_factor[sim.data.diff.sub2$Expand_Year %in% 2014] <- "2014"
sim.data.diff.sub2$Expand_Year_factor[sim.data.diff.sub2$Expand_Year %in% 2015] <- "2015"

## fitting model (4)
mod_att_overall_5.2 <- ols(Clinical_Retention ~ expanded + Year + Expand_Year_factor, data = sim.data.diff.sub2,x = TRUE, y = TRUE)

## estimate of ATT(2014,2014), ATT(2014,2015), and ATT(2015,2015)
att_overall_5.2 <- coefficients(mod_att_overall_5.2)["expanded"] %>% as.numeric() %>% round(3)

## calculating robust standard errors for ATT(2014,2014), ATT(2014,2015), and ATT(2015,2015) that accounts for clustering
robust_mod_att_overall_5.2 <- robcov(fit = mod_att_overall_5.2, cluster = sim.data.diff.sub2$NAID)
robust_se <- vcov(robust_mod_att_overall_5.2)["expanded","expanded"] %>% sqrt( )

## 95% CI for ATT(2014,2014), ATT(2014,2015), and ATT(2015,2015)
att_overall_5.2_lb <- (att_overall_5.2 - 1.96*robust_se) %>% round(3)
att_overall_5.2_ub <- (att_overall_5.2 + 1.96*robust_se) %>% round(3)


##### Table 1 #####

table1<- cbind(
c("Unadjusted",
"Race",
"Race-by-Time",
"Region",
"Region-by-Time",
"All Covariates",
"All Covaraites-by-Time"),

# unadjusted
rbind(
c(paste0(att_2014_2014_3.2, " (", att_2014_2014_3.2_lb, " , ", att_2014_2014_3.2_ub,")"), 
paste0(att_2014_2014_3.1, " (", att_2014_2014_3.1_lb, " , ", att_2014_2014_3.1_ub,")") ),

# race
c(paste0(att_2014_2014_race_4.2, " (", att_2014_2014_race_4.2_lb, " , ", att_2014_2014_race_4.2_ub,")"), 
  paste0(att_2014_2014_race_4.1, " (", att_2014_2014_race_4.1_lb, " , ", att_2014_2014_race_4.1_ub,")") ),
# race-by-time
c(paste0(att_2014_2014_race_by_time_4.2, " (", att_2014_2014_race_by_time_4.2_lb, " , ", att_2014_2014_race_by_time_4.2_ub,")"), 
  paste0(att_2014_2014_race_4.1, " (", att_2014_2014_race_4.1_lb, " , ", att_2014_2014_race_4.1_ub,")") ),

# region
c(paste0(att_2014_2014_region_4.2, " (", att_2014_2014_region_4.2_lb, " , ", att_2014_2014_region_4.2_ub,")"), 
  paste0(att_2014_2014_region_4.1, " (", att_2014_2014_region_4.1_lb, " , ", att_2014_2014_region_4.1_ub,")") ),
# region-by-time
c(paste0(att_2014_2014_region_by_time_4.2, " (", att_2014_2014_region_by_time_4.2_lb, " , ", att_2014_2014_region_by_time_4.2_ub,")"), 
  paste0(att_2014_2014_region_4.1, " (", att_2014_2014_region_4.1_lb, " , ", att_2014_2014_region_4.1_ub,")") ),
# all covariates
c(paste0(att_2014_2014_all_covariates_4.2, " (", att_2014_2014_all_covariates_4.2_lb, " , ", att_2014_2014_all_covariates_4.2_ub,")"), 
  paste0(att_2014_2014_all_covariates_4.1, " (", att_2014_2014_all_covariates_4.1_lb, " , ", att_2014_2014_all_covariates_4.1_ub,")") ),
# all covariates-by-time
c(paste0(att_2014_2014_all_covariates_by_time_4.2, " (", att_2014_2014_all_covariates_by_time_4.2_lb, " , ", att_2014_2014_all_covariates_by_time_4.2_ub,")"), 
  paste0(att_2014_2014_all_covariates_4.1, " (", att_2014_2014_all_covariates_4.1_lb, " , ", att_2014_2014_all_covariates_4.1_ub,")") ))) %>% as.data.frame()

names(table1) <- c("Covariates in Model", "Linear Regression Model", "Difference-in-differences")



## creating figure with Table 1 data
table1_dat <- rbind(c("Unadjusted",att_2014_2014_3.2, att_2014_2014_3.2_lb %>% round(2), att_2014_2014_3.2_ub %>% round(2), "Linear Regression Model"),
              c("Unadjusted", att_2014_2014_3.1, att_2014_2014_3.1_lb%>% round(2), att_2014_2014_3.1_ub%>% round(2), "Difference-in-differences"), 
              c("Race",att_2014_2014_race_4.2, att_2014_2014_race_4.2_lb %>% round(2), att_2014_2014_race_4.2_ub %>% round(2), "Linear Regression Model"),
              c("Race-by-time",att_2014_2014_race_by_time_4.2, att_2014_2014_race_by_time_4.2_lb %>% round(2), att_2014_2014_race_by_time_4.2_ub %>% round(2), "Linear Regression Model"), 
              c("Race",att_2014_2014_race_4.1, att_2014_2014_race_4.1_lb %>% round(2), att_2014_2014_race_4.1_ub %>% round(2), "Difference-in-differences"), 
              c("Region",att_2014_2014_region_4.2, att_2014_2014_region_4.2_lb %>% round(2), att_2014_2014_region_4.2_ub %>% round(2), "Linear Regression Model"),
              c("Region-by-time",att_2014_2014_region_by_time_4.2, att_2014_2014_region_by_time_4.2_lb %>% round(2), att_2014_2014_region_by_time_4.2_ub %>% round(2), "Linear Regression Model"), 
              c("Region",att_2014_2014_region_4.1, att_2014_2014_region_4.1_lb %>% round(2), att_2014_2014_region_4.1_ub %>% round(2), "Difference-in-differences"),
              c("All Covariates",att_2014_2014_all_covariates_4.2, att_2014_2014_all_covariates_4.2_lb %>% round(2), att_2014_2014_all_covariates_4.2_ub %>% round(2), "Linear Regression Model"),
              c("All Covariates-by-time",att_2014_2014_all_covariates_by_time_4.2, att_2014_2014_all_covariates_by_time_4.2_lb %>% round(2), att_2014_2014_all_covariates_by_time_4.2_ub %>% round(2), "Linear Regression Model"),
              c("All Covariates",att_2014_2014_all_covariates_4.1, att_2014_2014_all_covariates_4.1_lb %>% round(2), att_2014_2014_all_covariates_4.1_ub %>% round(2), "Difference-in-differences")) %>% as.data.frame()

names(table1_dat) <- c("variables","att","lb","ub","method")
table1_dat$att <- as.numeric(table1_dat$att)
table1_dat$lb <- as.numeric(table1_dat$lb)
table1_dat$ub <- as.numeric(table1_dat$ub)
table1_dat$method <- as.factor(table1_dat$method)
table1_dat$method <- relevel(table1_dat$method, ref = "Linear Regression Model")
table1_dat$variables <- as.factor(table1_dat$variables)
table1_dat$variables <- reorder(table1_dat$variables,levels = c("All Covariates-by-time","All Covariates","Region-by-time","Region","Race-by-time","Race","Unadjusted"))  


ggplot(data=table1_dat, aes(y=factor(variables, levels = c("All Covariates-by-time","All Covariates","Region-by-time","Region","Race-by-time","Race","Unadjusted")), x=att,xmin=lb, xmax=ub, label = paste0(att," (",lb,",",ub,")"))) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  #geom_text(vjust = -1,color = "dark grey",position = position_dodge(width = 1))+
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.25) + 
  theme_classic() + 
  facet_grid(cols = vars(method)) +
  labs(x = TeX("$\\widehat{ATT}(2014,2014)$"), y = "Covariates in Model")+
  theme(text = element_text(size=15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
