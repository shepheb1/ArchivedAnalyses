# This is the R file to replicate data pre-process
library(tidyverse)
library(lubridate)

tidydim <- function(x) {
  print(dim(x)[1])
  invisible(x)
}

load("./ccasanet_database_20220830/tbdata.Rda")

# [1] "patient"           "birth_d"           "enrol_d"           "site"
#     "center"            #"country"           "aids_first_d"      "aids_enrol_d"
# [9] "aids_art_d"        "patient_id"        "aids_first_y"      "aids_first_d_a"
#     "aids_art_y"        #"aids_enrol_y"      "aids_enrol_d_a"    "school_lvl"
# [17] "mode"              "mode_oth"          "male_y"            "hivdiagnosis_d_a"
#     "enrol_d_a"        # "recart_y"          "baseline_d_a"      "birth_d_a"
# [25] "clinicaltrial_y"   "pmtct"             "hivdiagnosis_d"    "recart_d"
#     "baseline_d"       # "birthmode"         "male"              "death_d"
# [33] "l_alive_d"         "death_y"           "siteclose"         "sc_minus_lastobs"
#     "ltfu"             # "baseage"           "min_tbdx"          "tb"
# [41] "tbminusbase"       "min_art_sd"        "min_haart_sd"      "artminusbase"
#     "artminusbase2"    # "baselinetb"        "death_before_art"  "futime"
# [49] "log10futime"       "total_follow"      "base_cd4_d"        "cd4_base"
#     "art_cd4_d"        # "cd4_art"           "aidsyear"          "hivyear"
# [57] "baseyear"          "school_lvl.factor" "Mode"              "HAART"
#     "HAART.factor"     # "time_to_HAART"     "tbsite"            "baselinetb_num"
# [65] "ART_event"         "cd4group"          "era"               "tb_before_HAART"
#     "cd4_base.imp"     # "logit.ps"          "fudeath"           "death_num"
# [73] "log10fudeath"      "log10_futime"

school_lvl = na.omit(unique(dt$school_lvl))
school_lvl_dic = list()
school_lvl_dic[school_lvl] = levels(dt$school_lvl.factor)[school_lvl]

### Old data

data_list = list.files("./ccasanet_database_20221005/",
                       pattern = ".csv",
                       full.names = T)

data_list = lapply(data_list, read_csv)

names(data_list) = list.files("./ccasanet_database_20221005/",
                              pattern = ".csv",
                              full.names = F)

#### Fix Haiti
haiti = list.files("./haiti_database_20240627/",
           pattern = ".csv",
           full.names = T) %>% 
  lapply(., read_csv)

names(haiti) = list.files("./haiti_database_20240627/",
                              pattern = ".csv",
                              full.names = F)

haiti_id = haiti$fvisit.csv$patient_id
haiti_id = intersect(haiti_id,haiti$basic.csv %>% 
                       filter(site == 'haiti',
                              data_checked == "Yes"
                              ) %>% 
                       pull(patient_id))

haiti$basic.csv = 
  haiti$basic.csv %>% 
  filter(site != "haiti"|(patient_id %in% haiti_id))

haiti$fvisit.csv =haiti$fvisit.csv %>% 
  group_by(patient_id) %>% 
  mutate(min_date = min(hivdiagnosis_date,baseline_date)) %>% 
  ungroup()

#### Fix haiti basic
haiti$basic.csv[haiti$basic.csv$patient_id %in% haiti_id,
                    c('enrol_d','baseline_d','hivdiagnosis_d')] =
  haiti$fvisit.csv[
    haiti$fvisit.csv$patient_id %in% haiti_id,
    c('min_date',"baseline_date","min_date")
  ]

#### Fix haiti follow
haiti$follow.csv[haiti$follow.csv$patient_id %in% haiti_id,
                     "death_date"] = 
  haiti$fvisit.csv[
    haiti$fvisit.csv$patient_id %in% haiti_id,
    "death_date"
  ]

#### Remove old haiti, add new haiti
data_list = lapply(data_list,function(x) filter(x,site != 'haiti'))
list_name = intersect(names(data_list),
                      names(haiti))
data_list = lapply(list_name,function(i){
  bind_rows(data_list[[i]],haiti[[i]])
})
names(data_list) = list_name


### Fix Peru

Peru_id = data_list$basic.csv %>% filter(site == 'peru') %>% pull(patient_id)
id_to_remove = haiti$Peru_id_to_remove.csv$patient_id
print(intersect(Peru_id,id_to_remove))

data_list$basic.csv = 
  filter(data_list$basic.csv,
         !patient_id %in% haiti$Peru_id_to_remove.csv$patient_id)

#### Fix HIV diangosis Date
data_list$basic.csv = data_list$basic.csv %>% 
  mutate(hivdiagnosis_d = if_else(hivdiagnosis_d == "1900-01-01",NA,
                                  hivdiagnosis_d))


#### The time of TB diagnosis table
first_tb = data_list$ce_tb.csv %>%
  group_by(center, patient_id) %>%
  filter(tbdiagnosis_d == min(tbdiagnosis_d)) %>%
  # keep the first TB records
  ungroup() %>%
  select(colnames(data_list$ce_tb.csv))


#### ART initition table
haart_init = data_list$art.csv %>%
  group_by(site, center, patient_id) %>%
  mutate(
    HAART = 1 * (art_class == "HAART"),
    HARRT.factor = as.factor(HAART),
    # Find Art is HAART
  ) %>%
  filter(art_class == "HAART") %>%
  filter(art_sd == min(art_sd)) %>%
  # keep the first ART treatment
  ungroup() %>%
  select(colnames(data_list$art.csv),
         art_sd,
         starts_with("HAART")) %>% 
  unique()


#### cd4 at enrollment table
baseline_cd4 =
  data_list$lab_cd4.csv %>%
  left_join(data_list$basic.csv) %>%
  mutate(
    date_diff = abs(cd4_d - enrol_d),
    # the days diff between diagnosis and enrollment
    is_baseline = 1 * ((date_diff - 90) <= 0),
    # Baseline CD4 count was defined using the closest CD4 count to
    # enrollment measured within a window of +/- 90 days
    # TODO: a lot of people don't have cd4 baseline
    cd4_base = cd4_v
  ) %>%  #within +-90 day windows
  # replace_na(list(is_baseline = 1)) %>%
  # # if not tb diagnosis, then they are baseline
  group_by(site, center, patient_id) %>%
  filter(is_baseline == 1, date_diff == min(date_diff)) %>% #keep the closest baseline cd4
  ungroup() %>%
  select(colnames(data_list$lab_cd4.csv), cd4_base) %>%
  distinct(patient_id, center, site, .keep_all = T)


#### cd4 at hiv-diagnosis table
hiv_cd4 =
  data_list$lab_cd4.csv %>%
  left_join(data_list$basic.csv) %>%
  mutate(
    date_diff = abs(cd4_d - hivdiagnosis_d),
    # the days diff between diagnosis and enrollment
    is_baseline = 1 * ((date_diff - 90) <= 0),
    # Baseline CD4 count was defined using the closest CD4 count to
    # enrollment measured within a window of +/- 90 days
    # TODO: a lot of people don't have cd4 baseline
    cd4_hiv = cd4_v
  ) %>%  #within +-90 day windows
  # replace_na(list(is_baseline = 1)) %>%
  # # if not tb diagnosis, then they are baseline
  group_by(site, center, patient_id) %>%
  filter(!is.na(date_diff)) %>% 
  filter(is_baseline == 1) %>% 
  filter(date_diff == min(date_diff)) %>% #keep the closest baseline cd4
  ungroup() %>%
  select(colnames(data_list$lab_cd4.csv), cd4_hiv) %>%
  distinct(patient_id, center, site, .keep_all = T)


########## rebuild the data
#### Inclusion
.new_data =
  data_list$basic.csv %>% # 57432 obs
  left_join(haart_init) %>%
  left_join(first_tb) %>%
  left_join(baseline_cd4) %>%
  left_join(hiv_cd4) %>% 
  left_join(data_list$follow.csv) %>%
  select(-ends_with("_d_a")) %>% 
  mutate(
    # indicator of art start 30 days before enrollment
    istartbefore30 = as.numeric(art_sd - enrol_d) > -30,
    istartbefore30 = ifelse(is.na(istartbefore30), T, istartbefore30),
    # indicator of TB happens 30 days - 0 days before enrollment
    istbbefore30 = as.numeric(tbdiagnosis_d - enrol_d) > -30,
    istbbefore30 = ifelse(is.na(istbbefore30), T, istbbefore30),
    
    # indicator of TB happens 180 days - 0 days before enrollment
    istbbefore180 = as.numeric(tbdiagnosis_d - enrol_d) > -180,
    istbbefore180 = ifelse(is.na(istbbefore180), T, istbbefore180)
    # TODO: some patient enroll after event
  ) %>% 
  mutate(close_date = 
           case_when(site == 'haiti' ~ '03-01-2022',
                     site == 'brazil' ~ '01-01-2021',
                     site == 'peru' ~ '02-01-2020',
                     site =='mexico' ~ '05-01-2022',
                     site == 'honduras' ~ '10-01-2021'),
         close_date = mdy(close_date)) %>% 
  group_by(site) %>% 
  mutate(recruit_2_close = as.numeric(enrol_d-close_date)
         )%>% 
  ungroup()

.new_data.copy = .new_data
#### Exclusion
.new_data =
  .new_data.copy%>%
  tidydim() %>% filter(
    site %in% str_to_lower(as.character(dt$site)),
    # the five county
  ) %>% tidydim() %>% 
  filter(
    as.numeric(enrol_d - birth_d) / 365.5>=18,
    # Age greater than 18
  ) %>% tidydim() %>% filter(
    enrol_d >= lubridate::ymd("2006-1-1"),
    # after 2006
  ) %>% tidydim() %>% filter(
    recruit_2_close <0, # New: 790 48407, old:388 50884 
  ) %>% tidydim() %>% filter(
    istartbefore30,  # New: 15910 41522 , old: 4901 54307
  ) %>% tidydim() %>% filter(
    # istbbefore180,
    istbbefore180, # New: 1675 55757 , old: 442 58766
    # we want to filter out patient who have tb way too early
  ) %>% tidydim()

## Old data [1] 34436    89
## New data [1] 21587   102

new_data = .new_data %>% 
  mutate(
    ## TB related variables
    istbbeforeart = as.numeric(art_sd - tbdiagnosis_d) > 0,
    istbbeforeart = ifelse(is.na(istbbeforeart), T, istbbeforeart),

    tb = between(as.numeric(tbdiagnosis_d - enrol_d),-180,30),
    tb = ifelse(is.na(tb), F, tb),

    # who initiated ART before (≥1 day) being diagnosed with TB were considered
    # to not have baseline TB.(this is not required anymore)
    # tb = tb * istbbeforeart,
    
    # if there's a record, you have TB
    baselinetb_num = tb,
    baselinetb = case_when(baselinetb_num > 0 ~ "TB",
                           baselinetb_num == 0 ~ "No TB"),

    tbminusbase = as.numeric(tbdiagnosis_d - enrol_d),
    
    ## baseline vairables
    site = str_to_title(site),
    
    baseyear = year(enrol_d),
    baseage = as.numeric(enrol_d - birth_d) / 365.5,
    aidyear = year(aids_first_d),
    hivyear = year(hivdiagnosis_d),
    
    male_y = ifelse(male_y > 1, NA, male_y),
    male = case_when(male_y > 0 ~ "Male",
                     male_y == 0 ~ "Female"),
    male = fct_relevel(as.factor(male),"Male"),
    
    school_lvl = ifelse(is.na(school_lvl), 9, school_lvl),
    school_lvl.factor = fct_recode(
      as.factor(school_lvl),
      "Unknown" = "9",
      "None" = "0",
      "Primary education" = "1",
      "Lower secondary or end of basic education" = "2",
      "Upper secondary or post-secondary non-tertiary" = "3",
      "University or post-graduate" = "4"
    ),
    
    ## follow up variable
    HAART = ifelse(is.na(HAART), 0, HAART),
    time_to_HAART = ifelse(HAART > 0, as.numeric(art_sd - enrol_d), NA),
    artminusbase = as.numeric(art_sd - enrol_d),
    death_num = death_y,
    death_y = as.factor(case_when(death_y > 0 ~ "Dead",
                                  death_y == 0 ~ "Alive")),
    
    ltfu = as.numeric(close_date-l_alive_d)>365,
    ltfu = ltfu*(1-death_num)>0,
    
    total_follow = as.numeric(l_alive_d - enrol_d),
    futime = as.numeric(art_sd - enrol_d),
    futime = case_when(between(futime, -30, 0) ~ 0,
                       futime >= 0 ~ futime),
    futime = case_when(HAART > 0 ~ futime,
                       HAART == 0 ~ as.numeric(l_alive_d - enrol_d)),
    # Patients who initiated ART within 30 days before enrollment were
    # classified as initiating ART at enrollment.
    fudeath = case_when(
      death_y == "Dead" ~ as.numeric(death_d - enrol_d),
      death_y == "Alive" ~ as.numeric(l_alive_d - enrol_d)
    ),
    fuhiv2art = as.numeric(art_sd- hivdiagnosis_d),
    fuhiv2art = case_when(between(fuhiv2art, -30, 0) ~ 0,
              fuhiv2art >= 0 ~ fuhiv2art),
    fuhiv2art = if_else(fuhiv2art>=0,fuhiv2art,NA),
    fuhiv2death = case_when(
      death_y == "Dead" ~ as.numeric(death_d - hivdiagnosis_d),
      death_y == "Alive" ~ as.numeric(l_alive_d - hivdiagnosis_d)
    ),
    fuhiv2death = if_else(fuhiv2death>=0,fuhiv2art,NA),
    
    log10futime = log10(futime + 1),
    log10fudeath = log10(fudeath + 1),
    log10fuhiv2art = log10(fuhiv2art+1),
    log10fuhiv2death = log10(fuhiv2death+1),
  )

raw_new_data = as.data.frame(new_data)
write_rds(select(raw_new_data,-patient_id, -center,-recart_d) #deidentify
          , "data/rebuilt_data_raw.rds")

new_data = select(new_data,intersect(colnames(new_data),colnames(dt)),
                  contains("hiv"),
                  -patient_id, -center,-recart_d) #deidentify

new_data = as.data.frame(new_data)
write_rds(new_data, "data/rebuilt_data.rds")
