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
  df_cm <- box_read("946034255301")
  df_ttse <- box_read("946034261301")
  df_hypo <-box_read("946034248101")
  df_adae <- box_read("946034244501")
  
}else{
  # read in subject level covariate dataset, tte, ttse and cm 
  df_adsl <- haven::read_sas("data/ADaM/adsl.sas7bdat")
  df_tte <- haven::read_sas("data/ADaM/adtte.sas7bdat")
  df_cm <- haven::read_sas("data/ADaM/adcm.sas7bdat")
  df_ttse <- haven::read_sas("data/ADaM/adttse.sas7bdat")
  df_hypo <- haven::read_sas("data/ADaM/adhypo.sas7bdat")
  df_adae <- haven::read_sas("data/ADaM/adae.sas7bdat")
  # df_adlb <- haven::read_sas("data/ADaM/adlb.sas7bdat")
  # df_adlbx <- haven::read_sas("data/ADaM/adlbx.sas7bdat")
  
}



### ------------  Part 1. Process Covariates  ------------ ###

# remove the subjects who are not in df_tte or df_cm
unique_subid_w = unique(df_adsl$USUBJID)
# 1 subject in subject level dataset, not in the time-to-event dataset
out_subid_tte <- unique_subid_w[which(!(unique_subid_w %in% unique(df_tte$USUBJID)))]
# 2 subjects in subject level dataset, not in the concomitant med dataset
out_subid_cm <- unique_subid_w[which(!(unique_subid_w %in% unique(df_cm$USUBJID)))]
out_subid <- union(out_subid_tte, out_subid_cm)


## -----  Part 1.1. Baseline Covariates  
W <- subset(df_adsl, !USUBJID %in% out_subid)


# REMOVE UNITS, ID VARIABLES, AND POST-BASELINE VARIABLES
W <- dplyr::select(W, -c(STUDYID, SUBJID, DIABDURU, WSTCRBLU,
                         BMIBLU, WGTBLU, WGTTBLU, PULSEBLU, SYSBPBLU,
                         DIABPBLU, HBA1CBLU, HDL1BLU, LDL1BLU,
                         CHOL1BLU, TRIG1BLU, HBC2BLU, HDL2BLU,
                         LDL2BLU, CHOL2BLU, TRIG2BLU, CREATBLU,
                         EGFREPIU, EGFREPBU, HGTTBLU, AGEU, HGTBLU)) %>%
  dplyr::select(-c(DTHDT, EOSDT, EOTDT, LSTCONDT, LSTSVDT,
                   BRTHDT, EOSSTT, TOTREAT, STDURY, TRDURD, TRDURY,
                   TRDU15OD, TRDUROY, TRDU15OY)) # ??? ISHRGRL, ISHRGRN, 

# REMOVE DEGENERATE COLUMNS
W <- select_if(W, ~length(unique(.))>1)

# REMOVE DUPLICATE COLUMNS & LINEARLY DEPENDENT COLUMNS
W <- dplyr::select(W, -c(PREVTFL, HBC2BL, ARMCD, # ??? WGTBL,
                         HDL2BL, LDL2BL, CHOL2BL, TRIG2BL,
                         GERDCMFL, # near dup GERDBLFL, off by 2
                         STROKSFL, # remove the sensitivity analysis stroke flag?,
                         CVHIFL, CVMEDFL, CVRISKN, age_category))

W <- W %>%
  mutate(RETINSEV = case_when(RETINSEV == "" ~ "retinsev NA",
                              T ~ RETINSEV),
         RETINSEV = factor(RETINSEV, ordered = T,
                           levels = c("retinsev NA", "no retinopathy", "non-proliferative", "proliferative")),
         NYHACLAS = factor(case_when(NYHACLAS == "" ~ "NYHA CLAS NA",
                                     T ~ NYHACLAS), ordered = T,
                           levels = c("NYHA CLAS NA", "NYHA CLASS I", "NYHA CLASS II", "NYHA CLASS III")))

# make logical variables logicals
W <- W %>%
  dplyr::mutate_if(~mean(unique(.) %in% c("Y", "N", "y", "n", "", NA)) == 1,
                   ~case_when(. %in% c("Y", "y") ~ TRUE,
                              . %in% c("N", "n") ~ FALSE,
                              T ~ NA)) %>%
  mutate(RACE = as.factor(case_when(RACE == "BLACK OR AFRICAN AMERICAN" ~ "BLACK",
                                    RACE %in% c("AMERICAN INDIAN OR ALASKA NATIVE",
                                                "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER",
                                                "OTHER") ~ "OTHER",
                                    T ~ RACE)),
         SEX = as.factor(SEX),
         ARM = as.numeric(ARM == "Liraglutide")
  ) 


# make factor variables factors, make SMOKER and renal disease cols ordered
W <- W %>% mutate(SMOKER = factor(SMOKER, ordered = T,
                                  levels = c("NEVER SMOKED", "PREVIOUS SMOKER",
                                             "CURRENT SMOKER"))) %>%
  mutate_at(vars(starts_with("RNF"), starts_with("RENF")),
            ~factor(., ordered = T,
                    levels = c("Normal (EGFR>=90)", "Mild (EGFR<90)",
                               "Moderate (EGFR<60)", "Severe (EGFR<30)"))) %>%
  mutate_if(is.character, as_factor)


## -----  Part 1.2. Concomitant medications  
# identify statin use and concomitant flag
df_cm <- df_cm %>%  
  filter(ANL01FL == "Y" & CONBLFL == "Y" & FASFL == "Y")

df_statin <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "C10AA|C10BA02|C10BA05|A10BH51|A10BH52"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    statin_use = TRUE) 


# Low dose aspirin 0
df_ldasp <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "B01AC06|N02BA01"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    low_dose_aspirin = TRUE) 


# Antihypertensives
df_antihypert <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "C09A|C09B|C09C|C09D|C09XA52|C10BX04|C10BX06|
                    C10BX07|C10BX10|C10BX11|C10BX12|C10BX13|C10BX14|C10BX15"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    antihypertensives = TRUE) 

# Beta blockers
df_betab <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "C07"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    betab = TRUE) 


# Mineralocorticoid receptor antagonists
df_minera <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "C03D"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    minera = TRUE) 

# ADP (any ADP, including ASA)
df_adp <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "B01AC"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    adp = TRUE) 

# NOACS
df_noacs <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "B01AF"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    noacs = TRUE) 

# Vitamin K antagonists
df_vkantag <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "B01AA"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    vkantag = TRUE) 

# Ca antagonists 
df_caantag <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "C08")))  %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    caantag = TRUE) 

# Digoxin 0
df_digoxin <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "C01AA05"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    digoxin = TRUE) 

# Thiazide 
df_thiazide <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "C03A"))) %>%
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    thiazide = TRUE) 

# Loop diuretics 
df_loopdiur <- df_cm %>% 
  filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
            any_vars(str_detect(., "C03EB|C03C"))) %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  group_by(USUBJID) %>% 
  summarise(
    loopdiur = TRUE) 


# Cross check numbers with SAS output
# df_adcm5 <- read_csv("data/supp/AD_CM_5.csv")
# df_adcm5 <- df_adcm5 %>% 
#   filter(!is.na(cat)) %>% 
#   group_by(USUBJID, cat, subcat) %>% 
#   filter(row_number() == 1) 
# 
# df_adcm5 %>% group_by(cat, subcat) %>% summarise(n = n())



# add new concomtant medication variables to W 
list_df_cm <- list(df_statin, df_antihypert, df_betab,
                   df_minera, df_adp, df_vkantag, 
                   df_caantag, df_thiazide, df_loopdiur)

for (df_i in list_df_cm){
  W <- left_join(W, df_i, by = 'USUBJID')
}

cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

W <- W %>% mutate_at(cm_names, ~ifelse(is.na(.), FALSE, TRUE))


## -----  Part 1.3. HDL, NLR ratio

W <- W %>% mutate(RC = CHOL1BL - LDL1BL - HDL1BL,
                  RCoverHDL  = RC/HDL1BL)

# df_hdl <- df_adlbx %>% filter(AVISIT %in% c("VISIT 1, V10",
#                                            "VISIT 2, V20",
#                                            "VISIT 3 (DAY 0), V30"))
# write_csv(df_hdl, file = "data/supp/hdl.csv")


# LDL
# HDL
# TRIG
# CHOL
# RC = TC - LDL - HDL
# NonHDL = TC - HDL
# TGoverHDL = TRIG/HDL
# TCoverHDL = TC/HDL
# LDLoverHDL = LDL/HDL
# RC/HDL-C  = (TC - LDL - HDL)/HDL
# nonHDL-C/HDL-C = (TC - HDL)/HDL


# df_hdl <- read_csv("data/supp/hdl.csv")
# 
# df_hdl_v3 <- df_hdl %>% 
#   filter(AVISIT == "VISIT 3 (DAY 0), V30" &
#            PARAMCD %in% c("LDL", "HDL", "TRIG", "CHOL", "LDLE") &
#            !is.na(AVAL))
# 
# df_hdl_v3 <- df_hdl_v3 %>% 
#   group_by(USUBJID, PARAMCD) %>% 
#   filter(row_number() == 1)
# 
# # df_hdl_v1 <- df_hdl %>% filter(AVISIT == "VISIT 1, V10")
# # unique(df_hdl_v3$PARAMCD)
# df_hdl_v3 <- df_hdl_v3 %>% 
#   pivot_wider(id_cols = "USUBJID", 
#               names_from = "PARAMCD", 
#               values_from = c("AVAL")) %>% 
#   mutate(LDLE = coalesce(LDLE, LDL),
#          RC = CHOL - LDLE - HDL,
#          RCoverHDL  = RC/HDL) %>% 
#   select(-"LDL")
# 
# W <- left_join(W, df_hdl_v3, by = 'USUBJID')






## -----  Part 1.4. Final screening, drop/collapse redundant and sparse W
# summary table W
# manually drop some redundant W
df_base <- W %>% select(USUBJID, RANDDT)

W <- W %>% 
  select(c("USUBJID", "ARM", "AGE", "SEX", "RACE", "COUNTRY", "SMOKER", "NYHACLAS",
           "DIABDUR", "INSNVFL", "ANTDBFL", "AHYPERFL", "INCPASSN", 
           "BMIBL", "PULSEBL", "SYSBPBL", "DIABPBL", "HBA1CBL", "HDL1BL", "LDL1BL",
           "CHOL1BL", "RC", "RCoverHDL","TRIG1BL", "CREATBL", "EGFMDRBC", 
           "RETINSEV", "GERDBLFL", "PPIFL", "H2BLFL", 
           "MIFL", "STROKEFL", "REVASFL", "STENFL", "CHDFL", "IHDFL", "CHFFL",
           "KIDFL", "MICFL", "HYPFL", "LVSDFL", "PADFL", "CVRISK", "HBA1CGRN", "DDURGRN", 
           "statin_use", "antihypertensives", "betab", "minera", "adp",
           "vkantag", "caantag", "thiazide", "loopdiur")) %>% 
  rename("A" = "ARM")

# Union subcat of COUNTRY
W <- W %>% mutate(COUNTRY = case_when(str_detect(COUNTRY, "America") ~ "America",
                                      str_detect(COUNTRY, "Asia") ~ "Asia",
                                      str_detect(COUNTRY, "Africa") ~ "Africa",
                                      str_detect(COUNTRY, "Europe") ~ "Europe",
                                      str_detect(COUNTRY, "Pacific") ~ "Pacific"))

# write_csv(W, file = "data/supp/df_w.csv")


## -----  Part 1.5. Table of Baseline Characteristics 
df_w_summary <- W %>% select(-c("USUBJID"))

# cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
#               "vkantag", "caantag", "thiazide", "loopdiur")

df_w_summary <- df_w_summary %>%
  mutate(A = ifelse(A == 1, "Liraglutide", "Placebo"))

df_w_summary <- labelled::remove_labels(df_w_summary)

library(table1)
library(r2rtf)

tbl <- table1(~. | A, data = df_w_summary)
tbl



### ------------  Part 2. Process Outcomes  ------------ ###

##------------  Part 2.1 Diabetes progression
# First time insulin
Y_insulin <- df_ttse %>%
  dplyr::select("USUBJID", 'PARAMCD', "AVAL", "PARCAT1") %>%
  dplyr::filter(PARAMCD %in% c("FRINSLTM")) %>%
  rename("EVENT" = "PARCAT1") %>% 
  pivot_wider(id_cols = "USUBJID", names_from = "PARAMCD", values_from = c("AVAL", "EVENT")) %>%
  mutate_at(vars(starts_with("EVENT")),
            ~as.numeric(. == "FIRST INSULIN EVENT"))

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

# write_csv(Y_diab, file = "data/supp/Y_diab.csv")


##------------  Part 2.1 cv outcomes
Y_cv <- df_tte %>%
  dplyr::select("USUBJID", 'PARAMCD', "AVAL", "PARCAT1") %>%
  dplyr::filter(PARAMCD %in% c("MACEEVTM", "MCECVDTM", "MACEMITM",
                               "MCENFSTM", "NONCVTM"# , "PRMACETM",
  )) %>%
  rename("EVENT" = "PARCAT1") %>%
  mutate(PARAMCD = case_when(PARAMCD == "MACEEVTM" ~ "MACE",
                             PARAMCD == "MACEMITM" ~ "NFMI",
                             PARAMCD == "MCECVDTM" ~ "CVDEATH",
                             PARAMCD == "MCENFSTM" ~ "NFSTROKE",
                             PARAMCD == "NONCVTM" ~ 'NONCVDEATH' #,
                             # PARAMCD == "PRMACETM" ~ 'PRIMACE'
  )) %>%
  pivot_wider(id_cols = "USUBJID", names_from = "PARAMCD", values_from = c("AVAL", "EVENT")) %>%
  mutate_at(vars(starts_with("EVENT")),
            ~as.numeric(. == "TIME TO EVENT"))

# check delta
temp <- Y_cv %>% 
  select("USUBJID", "AVAL", "EVENT") %>% 
  filter(EVENT == "LAST CONTACT DATE") 

temp <- temp[!duplicated(temp), ]

temp <- temp %>% filter(USUBJID %in% W$USUBJID[W$INSNVFL == FALSE])

# hist(temp$AVAL)
# temp_sub <- temp %>% filter(temp$AVAL <= 24)

# df_id <- W %>% select(USUBJID)
# Y_cv <- left_join(df_id, Y_cv, by = "USUBJID")
# write_csv(Y_cv, file = "data/supp/Y_cv.csv")


##------------  Part 2.3 HBA1C
# df_a1c <- df_adlb %>% 
#   filter(PARAMCD == "HBA1C" & str_starts(AVISIT, "VISIT ([3-9]|1[0-3])"))
# write_csv(df_a1c, file = "data/supp/a1c.csv")
df_a1c <- read_csv("data/supp/a1c.csv")

df_a1c <- df_a1c %>%
  group_by(USUBJID, AVISIT) %>%
  filter(row_number() == 1) %>% 
  select(c("USUBJID", "AVISIT", "AVAL"))

df_a1c <- df_a1c %>%
  mutate(AVISIT = str_extract(AVISIT, "VISIT ([3-9]|1[0-3])")) %>% 
  mutate(AVISIT = str_replace(AVISIT, " ", "_")) %>% 
  mutate(AVISIT = paste0("HBA1C_", AVISIT)) %>% 
  pivot_wider(id_cols = "USUBJID",
              names_from = "AVISIT",
              values_from = c("AVAL")) %>% 
  select(-HBA1C_VISIT_4)

# df_id <- W %>% select(USUBJID)
# Y_a1c <- left_join(df_id, df_a1c, by = "USUBJID")
# write_csv(Y_a1c, file = "data/supp/Y_a1c.csv")





### ------------  Part 3. Export Dataset  ------------ ###

df_w <- read_csv(file = "data/supp/df_w.csv")

# 3.1 diabetes progression
Y_diab <- read_csv(file = "data/supp/Y_diab.csv")
df_all_diab <- left_join(df_w, Y_diab, by = 'USUBJID')

# write_csv(df_all_diab, file = "data/supp/df_all_diab.csv")

# 3.2 cv outcomes
Y_cv <- read_csv(file = "data/supp/Y_cv.csv")
df_all_cv <- left_join(df_w, Y_cv, by = 'USUBJID')

# write_csv(df_all_cv, file = "data/supp/df_all_cv.csv")

# 3.2 HBA1C
Y_a1c <- read_csv(file = "data/supp/Y_a1c.csv")
df_all_a1c <- left_join(df_w, Y_a1c, by = 'USUBJID')

# write_csv(df_all_a1c, file = "data/supp/df_all_a1c.csv")
