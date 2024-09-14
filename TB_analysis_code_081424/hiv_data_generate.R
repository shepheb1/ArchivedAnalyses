source("./data_generation.R")

.new_data = 
  .new_data.copy %>% 
  select(-ends_with("_d_a"),-contains("_rs")) %>% 
  mutate(
    # indicator of art start 30 days before enrollment
    istartbefore30 = as.numeric(art_sd - hivdiagnosis_d) > -30,
    istartbefore30 = ifelse(is.na(istartbefore30), T, istartbefore30),
    # indicator of TB happens 30 days - 0 days before enrollment
    istbbefore30 = as.numeric(tbdiagnosis_d - hivdiagnosis_d) > -30,
    istbbefore30 = ifelse(is.na(istbbefore30), T, istbbefore30),
    
    # indicator of TB happens 180 days - 0 days before enrollment
    istbbefore180 = as.numeric(tbdiagnosis_d - hivdiagnosis_d) > -180,
    istbbefore180 = ifelse(is.na(istbbefore180), T, istbbefore180)
    # TODO: some patient enroll after event
  ) %>% 
  filter(
    as.numeric(hivdiagnosis_d - birth_d) / 365.5>=18,
    # Age greater than 18
    
    hivdiagnosis_d > lubridate::ymd("2006-1-1"),
    enrol_d > lubridate::ymd("2006-1-1"),
    # after 2006
    
    recruit_2_close <0,
    
    site %in% str_to_lower(as.character(dt$site)),
    # the five county
    
    istartbefore30,
    recart_y < 1,
    
    # istbbefore30,
    istbbefore180,
    # we want to filter out patient who have tb way too early
    
  ) 

new_data = .new_data %>% 
  mutate(
    ## TB related variables
    istbbeforeart = as.numeric(art_sd - tbdiagnosis_d) > 0,
    istbbeforeart = ifelse(is.na(istbbeforeart), T, istbbeforeart),
    
    tb = between(as.numeric(tbdiagnosis_d - hivdiagnosis_d),-180,30),
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
    
    total_follow = as.numeric(l_alive_d - hivdiagnosis_d),
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

follow_tb = 
  data_list$ce_tb.csv %>% 
  group_by(patient_id,site,center) %>% 
  summarise(tb_follow = 1, # this include baseline tb, need to further condition
            tb_follow_d = min(tbdiagnosis_d)
  ) %>% 
  ungroup()

art_tb =
  .new_data %>%
  select(patient_id,site,center,art_sd) %>%
  left_join(data_list$ce_tb.csv) %>%
  mutate(tb_art = as.numeric(art_sd-tbdiagnosis_d)>=-30) %>% # 30 days after
  group_by(patient_id,site,center) %>%
  summarise(tb_art = 1*(sum(tb_art)>0),
            tb_art = replace_na(tb_art,0))

art_drug = 
  data_list$art.csv %>% 
  left_join(follow_tb) %>% 
  group_by(patient_id,site,center) %>% 
  mutate(duration = as.numeric(art_sd - min(art_sd)),
         time_to_tb = as.numeric(tb_follow_d - art_sd), # we want the positive value
         time_to_tb = replace_na(time_to_tb,0)
  ) %>% 
  filter(duration <= 180 & time_to_tb>=0) %>% 
  summarise(insti = 1*(sum(ii1,ii2)>0), #ii if any ii1 or ii2
            nnrti = 1*(sum(nnrti)>0), #nnrti if any nnrti
            pi = 1*(sum(pi)>0)       #pi if any pi
  ) %>% 
  ungroup()

first_art_drug = 
  .new_data %>% #do we want to use only 1 treatment record or all treatment?
  group_by(patient_id,site,center) %>% 
  summarise(insti_first = 1*(sum(ii1,ii2)>0), #ii if any ii1 or ii2
            nnrti_first = 1*(sum(nnrti)>0), #nnrti if any nnrti
            pi_first = 1*(sum(pi)>0),       #pi if any pi
            insti_only = case_when(
              insti_first>0 & sum(nnrti_first,pi_first) ==0 ~ "insti only",
              sum(nnrti_first,pi_first) >0 ~ "other"
            ),
  ) %>% 
  ungroup()

cd4_at_art = 
  data_list$lab_cd4.csv %>% 
  left_join(select(.new_data,patient_id,site,center,art_sd)) %>% 
  # add haart starting date
  mutate(is_art_cd4 = abs(as.numeric(art_sd - cd4_d))) %>% 
  group_by(patient_id,site,center) %>%
  filter(is_art_cd4 == min(is_art_cd4)) %>%  # keep the closest cd4 record to art
  filter(is_art_cd4 <=90) %>% # if the cd4 is taken within a 90 days
  summarise(cd4_haart = mean(cd4_v,na.rm = T)) %>% 
  ungroup()

#   summarise(tb_art = replace_na(1*(sum(tb_art)>0),0))

# making final new_data

raw_new_data = raw_new_data %>% 
  left_join(art_drug %>% mutate(site = str_to_title(site))) %>% # adding art type
  left_join(first_art_drug %>% mutate(site = str_to_title(site))) %>% 
  left_join(follow_tb %>% mutate(site = str_to_title(site))) %>% # adding followup tb
  left_join(cd4_at_art %>% mutate(site = str_to_title(site))) %>% 
  left_join(art_tb %>% mutate(site = str_to_title(site))) %>% 
  mutate(
    tb_follow = replace_na(tb_follow, 0),
    # adding follow up tb info
    art_before_tb = replace_na(art_sd<tb_follow_d,1),
    # some patietn develop tb before receiving haart
    tb_follow_num =
      ifelse(baselinetb_num == 0 & tb_follow == 1, 1, 0),
    #only no baseline and then follow up can be follow up tb
    # tb_follow_num = 
    #   ifelse(tb_follow_num*art_before_tb==1,1,0),
    tb_follow = case_when(tb_follow_num > 0 ~ "TB",
                          tb_follow_num == 0 ~ "No TB"),
    futb = case_when(
      tb_follow_num== 1 ~ as.numeric(tb_follow_d - art_sd),
      tb_follow_num == 0 ~ as.numeric(l_alive_d - art_sd),
      is.na(art_sd) ~ Inf
    ),
    log10futb = log10(futb+1)
  ) 
raw_new_data = raw_new_data %>% 
  select(
    patient_id,
    birth_d,
    enrol_d,
    site,
    center,
    country,
    recart_d,
    aids_first_d,
    aids_enrol_d,
    aids_art_d,
    aids_first_y,
    aids_art_y,
    aids_enrol_y,
    school_lvl,
    mode,
    mode_oth,
    male_y,
    recart_y,
    clinicaltrial_y,
    hivdiagnosis_d,
    baseline_d,
    birthmode,
    male,
    HAART,
    tbsite,
    cd4_base,
    death_d,
    l_alive_d,
    death_y,
    tb,
    baselinetb_num,
    baselinetb,
    tbminusbase,
    baseyear,
    baseage,
    hivyear,
    school_lvl.factor,
    time_to_HAART,
    artminusbase,
    death_num,
    ltfu,
    total_follow,
    futime,
    fudeath,
    log10futime,
    log10fudeath,
    cd4_hiv,
    fuhiv2art,
    fuhiv2death,
    futb,
    log10fuhiv2art,
    log10fuhiv2death,
    log10futb
  )


#deidentify
new_data = select(raw_new_data,-patient_id, -center,-recart_d) 
new_data = as.data.frame(new_data)
write_rds(new_data, "data/hiv_data.rds")