source("data_generation.R")
#tb follow up data

follow_tb = 
  .new_data %>% 
  select(patient_id,enrol_d) %>% 
  left_join(data_list$ce_tb.csv) %>% 
  filter(as.numeric(tbdiagnosis_d-enrol_d)>30) %>%
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

new_data = raw_new_data %>% 
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
    fuenrol2tb = case_when(
      tb_follow_num== 1 ~ as.numeric(tb_follow_d - enrol_d),
      tb_follow_num == 0 ~ as.numeric(l_alive_d - enrol_d),
    ),
    log10futb = log10(futb+1)
  )



#deidentify
new_data = select(new_data,-patient_id, -center,-recart_d) 

new_data = as.data.frame(new_data)

write_rds(new_data,"data/tb_follow_up_data.rds")