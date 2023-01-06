library(here)
library(tidyverse)
library(dplyr)
library(mice)
rm(list = ls())


### ------------  Part 0. Import Dataset  ------------ ###

if(here::here()=="C:/Users/andre/Documents/te_vim"){
  library(boxr)
  box_auth()
  
  df_adsl <- box_read("946034250501")
  df_tte <- box_read("946034258901")
  df_adcm <- box_read("946034255301")
  df_ttse <- box_read("946034261301")
  df_hypo <-box_read("946034248101")
  df_adae <- box_read("946034244501")
  
}else{
  # read in subject level covariate dataset, tte, ttse and cm 
  df_adsl <- haven::read_sas("data/ADaM/adsl.sas7bdat")
  df_tte <- haven::read_sas("data/ADaM/adtte.sas7bdat")
  df_adcm <- haven::read_sas("data/ADaM/adcm.sas7bdat")
  df_ttse <- haven::read_sas("data/ADaM/adttse.sas7bdat")
  df_hypo <- haven::read_sas("data/ADaM/adhypo.sas7bdat")
  df_adae <- haven::read_sas("data/ADaM/adae.sas7bdat")
  # df_adlb <- haven::read_sas("data/ADaM/adlb.sas7bdat")
  # df_adlbx <- haven::read_sas("data/ADaM/adlbx.sas7bdat")
  
}

df_base <- read_csv("data/supp/df_base.csv")

### ------------  Part 2. Process Outcomes  ------------ ###

##------------  Part 2.1 Diabetes progression
# Replace first-time insulin with insulin intensification (code from Edwin, def from Kim)
#Find patients on insulin
adcm_insulin <- df_adcm %>%
  filter(FASFL == "Y" & ANL01FL == "Y") %>% 
  filter(grepl("A10A", CMCLASCD, ignore.case = T)) %>% 
  select(USUBJID, TRTP,CMCLASCD, ASTDY, AENDY,ASTRF,AENRF,ASTDT)

#Separate baseline & after randomization uses of insulin
insulin_baseline <- adcm_insulin %>% 
  filter((ASTRF == 'BEFORE' & AENRF != 'BEFORE') | (is.na(ASTRF) & AENRF != 'BEFORE')) %>%
  select(-ASTRF, -AENRF) %>% 
  mutate(baseline1 = case_when(CMCLASCD %in% c('A10AB', 'A10AD') ~ TRUE),
         baseline2 = TRUE)

insulin_after1 <- adcm_insulin %>% 
  filter((ASTRF == 'AFTER' |ASTRF == 'DURING')) %>%
  select(-ASTRF, -AENRF) %>% 
  filter(CMCLASCD %in% c('A10AB', 'A10AD'))

insulin_after2 <- adcm_insulin %>% 
  filter((ASTRF == 'AFTER' |ASTRF == 'DURING')) %>%
  select(-ASTRF, -AENRF) %>% 
  filter(!CMCLASCD %in% c('A10AB', 'A10AD'))


insulin_after1_1 <- insulin_after1 %>%
  left_join(insulin_baseline  %>% 
              select(-ASTDY, -AENDY, -ASTDT), 
            by = c('USUBJID', 'CMCLASCD', 'TRTP')) %>% 
  filter(is.na(baseline1)) 

insulin_after2_1 <- insulin_after2 %>%
  left_join(insulin_baseline  %>% 
              select(-ASTDY, -AENDY, -ASTDT), 
            by = c('USUBJID', 'CMCLASCD', 'TRTP')) %>% 
  filter(is.na(baseline2)) 

insulin_after <- rbind(insulin_after1_1, insulin_after2_1)

#Identify 'new initiators'
Y_insulin <- insulin_after %>%
  select(-baseline1, -baseline2) %>% #Remove patients on "premix" or "short" at baseline
  group_by(USUBJID, TRTP) %>% 
  arrange(ASTDT, .by_group = TRUE) %>% 
  group_by(USUBJID, TRTP) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  rename('DATE_NEWINS' = 'ASTDT') %>% 
  select(USUBJID, DATE_NEWINS) %>% 
  mutate(DATE_NEWINS = as.Date(DATE_NEWINS))

# #Identify 'new initiators'
# Y_insulin <- insulin_after %>%
#   left_join(insulin_baseline  %>% 
#               select(-ASTDY, -AENDY, -ASTDT), 
#             by = c('USUBJID', 'CMCLASCD', 'TRTP')) %>% #CMCLASCD encodes subcategory of insulin
#   filter(is.na(baseline)) %>% 
#   select(-baseline) %>% #Remove patients on "premix" or "short" at baseline
#   group_by(USUBJID, TRTP) %>% 
#   arrange(ASTDT, .by_group = TRUE) %>% 
#   group_by(USUBJID, TRTP) %>% 
#   filter(row_number()==1) %>% 
#   ungroup() %>% 
#   rename('DATE_NEWINS' = 'ASTDT') %>% 
#   select(USUBJID, DATE_NEWINS) %>% 
#   mutate(DATE_NEWINS = as.Date(DATE_NEWINS))


# First time OAD
df_OADTM <- read_csv("data/supp/FINAL1.csv")
# df_OADINSTM <- read_csv("data/supp/FINAL2.csv")
Y_oad <- df_OADTM %>% 
  dplyr::select("USUBJID", 'PARAMCD', "AVAL", "PARCAT1") %>%
  dplyr::filter(PARAMCD %in% c("OADTM")) %>%
  rename("EVENT" = "PARCAT1") %>% 
  pivot_wider(id_cols = "USUBJID", names_from = "PARAMCD", values_from = c("AVAL", "EVENT")) %>%
  mutate_at(vars(starts_with("EVENT")),
            ~as.numeric(. == "TIME (MONTHS) TO FIRST OAD EVENT"))

# Hyperglycaemia
# aval = AESTDTC - ?
Y_hyper <- df_adae %>% 
  dplyr::filter(AEDECOD %in% c("Hyperglycaemia") & AEONTRFL == "Y") %>% 
  dplyr::select("USUBJID", 'AESTDTC') %>% 
  group_by(USUBJID) %>% filter(row_number() == 1) %>% 
  rename('DATE_HYPER' = 'AESTDTC') %>% 
  mutate(DATE_HYPER = as.Date(DATE_HYPER))

# Hypo
# aval = ASTDT - ?
Y_hypo <- df_hypo %>% 
  dplyr::filter(ACAT01  %in% c("SEVERE") & XHONTRFL == "Y") %>% 
  dplyr::select("USUBJID", 'ASTDT') %>% 
  group_by(USUBJID) %>% filter(row_number() == 1) %>% 
  rename('DATE_HYPO' = 'ASTDT') %>% 
  mutate(DATE_HYPO = as.Date(DATE_HYPO))

# Acid-base disorders
Y_keto <- df_adae %>% 
  dplyr::filter(AEHLGT == "Acid-base disorders" & AEONTRFL == "Y") %>%
  dplyr::select("USUBJID", 'AESTDTC') %>% 
  rename('DATE_ACIDOSIS' = 'AESTDTC') %>% 
  mutate(DATE_ACIDOSIS = as.Date(DATE_ACIDOSIS))



# combine all diabete outcomes
Y_diab <- df_base
list_Y_diab <- list(Y_insulin, Y_oad, Y_hyper, Y_hypo, Y_keto)
for (df in list_Y_diab){
  Y_diab <- left_join(Y_diab, df, by = "USUBJID")
}

Y_diab$AVAL_FRINSLTM <- as.numeric(Y_diab$DATE_NEWINS - Y_diab$RANDDT)/(365.25/12)
Y_diab$EVENT_FRINSLTM <- as.numeric(!is.na(Y_diab$AVAL_FRINSLTM))
Y_diab$AVAL_HYPER <- as.numeric(Y_diab$DATE_HYPER - Y_diab$RANDDT)/(365.25/12)
Y_diab$EVENT_HYPER <- as.numeric(!is.na(Y_diab$AVAL_HYPER))
Y_diab$AVAL_HYPO <- as.numeric(Y_diab$DATE_HYPO - Y_diab$RANDDT)/(365.25/12)
Y_diab$EVENT_HYPO <- as.numeric(!is.na(Y_diab$AVAL_HYPO))
Y_diab$AVAL_ACIDOSIS <- as.numeric(Y_diab$DATE_ACIDOSIS - Y_diab$RANDDT)/(365.25/12)
Y_diab$EVENT_ACIDOSIS <- as.numeric(!is.na(Y_diab$AVAL_ACIDOSIS))

Y_diab <- Y_diab %>% select("USUBJID", "AVAL_FRINSLTM", "AVAL_OADTM",
                            "AVAL_HYPER", "AVAL_HYPO", "AVAL_ACIDOSIS", 
                            "EVENT_FRINSLTM", "EVENT_OADTM", "EVENT_HYPER",
                            "EVENT_HYPO", "EVENT_ACIDOSIS")

Y_diab <- Y_diab %>% 
  mutate_at(vars(starts_with("AVAL")), ~ifelse(is.na(.), 999, .)) %>% 
  mutate_at(vars(starts_with("EVENT")), ~ifelse(is.na(.), 0, .))


### ------------  Part 3. Export Dataset  ------------ ###

df_w <- read_csv(file = "data/supp/df_w.csv")

# 3.1 diabetes progression
df_all_diab2 <- left_join(df_w, Y_diab, by = 'USUBJID')

# write_csv(df_all_diab2, file = "data/supp/df_all_diab2.csv")
