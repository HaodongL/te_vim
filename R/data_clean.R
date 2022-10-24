library(here)
library(tidyverse)
library(dplyr)
library(mice)


### ------------  Part 0. Import Dataset  ------------ ###

# read in subject level covariate dataset, tte, ttse and cm 
df_w <- haven::read_sas("data/ADaM/adsl.sas7bdat")
df_tte <- haven::read_sas("data/ADaM/adtte.sas7bdat")
df_cm <- haven::read_sas("data/ADaM/adcm.sas7bdat")
df_ttse <- haven::read_sas("data/ADaM/adttse.sas7bdat")
df_hypo <- haven::read_sas("data/ADaM/adhypo.sas7bdat")
df_adae <- haven::read_sas("data/ADaM/adae.sas7bdat")
# df_adlb <- haven::read_sas("data/ADaM/adlb.sas7bdat")


### ------------  Part 1. Process Covariates  ------------ ###

# remove the subjects who are not in df_tte or df_cm
unique_subid_w = unique(df_w$USUBJID)
# 1 subject in subject level dataset, not in the time-to-event dataset
out_subid_tte <- unique_subid_w[which(!(unique_subid_w %in% unique(df_tte$USUBJID)))]
# 2 subjects in subject level dataset, not in the concomitant med dataset
out_subid_cm <- unique_subid_w[which(!(unique_subid_w %in% unique(df_cm$USUBJID)))]
out_subid <- union(out_subid_tte, out_subid_cm)

W <- subset(df_w, !USUBJID %in% out_subid)


# REMOVE UNITS, ID VARIABLES, AND POST-BASELINE VARIABLES
W <- dplyr::select(W, -c(STUDYID, SUBJID, DIABDURU, WSTCRBLU,
                         BMIBLU, WGTBLU, WGTTBLU, PULSEBLU, SYSBPBLU,
                         DIABPBLU, HBA1CBLU, HDL1BLU, LDL1BLU,
                         CHOL1BLU, TRIG1BLU, HBC2BLU, HDL2BLU,
                         LDL2BLU, CHOL2BLU, TRIG2BLU, CREATBLU,
                         EGFREPIU, EGFREPBU, HGTTBLU, AGEU, HGTBLU)) %>%
  dplyr::select(-c(DTHDT, EOSDT, EOTDT, LSTCONDT, LSTSVDT,
                   RANDDT, BRTHDT, EOSSTT, TOTREAT, STDURY, TRDURD, TRDURY,
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
  mutate(RETINSEV = case_when(RETINSEV == "" ~ "NA",
                              T ~ RETINSEV),
         RETINSEV = factor(RETINSEV, ordered = T,
                           levels = c("NA", "no retinopathy", "non-proliferative", "proliferative")),
         NYHACLAS = factor(case_when(NYHACLAS == "" ~ "NA",
                                     T ~ NYHACLAS), ordered = T,
                           levels = c("NA", "NYHA CLASS I", "NYHA CLASS II", "NYHA CLASS III")))

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


# identify statin use and concomitant flag
# temp <- adcm %>% group_by(SUBJID) %>% summarise(count_statin = sum(CMCLASCD == 'C10AA'))
# unique(temp$count_statin)
# temp <- adcm %>% filter(SUBJID == '15506', CMCLASCD == 'C10AA')

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
    minera = TRUE,
    minera_cm = sum(CONBLFL == 'Y') > 0) 

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
df_adcm5 <- read_csv("data/ADaM/AD_CM_5.csv")
df_adcm5 <- df_adcm5 %>% 
  filter(!is.na(cat)) %>% 
  group_by(USUBJID, cat, subcat) %>% 
  filter(row_number() == 1) 

df_adcm5 %>% group_by(cat, subcat) %>% summarise(n = n())



# add new concomtant medication variables to W 

list_df_cm <- list(df_statin, df_antihypert, df_betab,
                   df_minera, df_adp, df_vkantag, 
                   df_caantag, df_thiazide, df_loopdiur)

for (df_i in list_df_cm){
  W <- left_join(W, df_i, by = 'USUBJID')
}



# summary table W
# manually drop some redundant W
W <- W %>% 
  select(c("USUBJID", "ACTARM", "AGE", "SEX", "RACE", "COUNTRY", "SMOKER", 
           "DIABDUR", "INSNVFL", "INSNVFL", "ANTDBFL", "AHYPERFL", "EOTSTT", "INCPASSN", 
           "BMIBL", "PULSEBL", "SYSBPBL", "DIABPBL", "HBA1CBL", "HDL1BL", "LDL1BL",
           "CHOL1BL", "TRIG1BL", "CREATBL", "EGFREPI", "EGFRMDRC", "EGFREPB", "EGFMDRBC", 
           "EGFRMDR", "RETINSEV", "RETISEVN", "GERDBLFL", "PPIFL", "H2BLFL", 
           "MIFL", "STROKEFL", "REVASFL", "STENFL", "CHDFL", "IHDFL", "CHFFL",
           "KIDFL", "MICFL", "HYPFL", "LVSDFL", "PADFL", "CVRISK", "HBA1CGRN", "DDURGRN", 
           "statin_use", "antihypertensives", "betab", "minera", "minera_cm", "adp",
           "vkantag", "caantag", "thiazide", "loopdiur"))

df_w_summary <- W %>% select(-c("USUBJID"))

cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "minera_cm", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

df_w_summary <- df_w_summary %>% 
  mutate_at(cm_names, ~ifelse(is.na(.), FALSE, TRUE))



library(table1)
library(r2rtf)

tbl <- table1(~. | ACTARM, data = df_w_summary)

tbl


# library(papaja)
# # Restructure object
# x <- attr(tbl, "obj")$contents
# names(x) <- lapply(x, function(x){rownames(x)[[1L]]})
# x <- lapply(x, function(x){x[-1L, ]})
# # Use apa_table()
# apa_table(x, caption = "Output from table1 in a pdf document.")


### ------------  Part 2. Process Outcomes  ------------ ###

# diabetes outcome
# First time insulin
Y_insulin <- df_ttse %>%
  dplyr::select("USUBJID", 'PARAMCD', "AVAL", "PARCAT1") %>%
  dplyr::filter(PARAMCD %in% c("FRINSLTM")) %>%
  rename("EVENT" = "PARCAT1") 
# pivot_wider(id_cols = "USUBJID", names_from = "PARAMCD", values_from = c("AVAL", "EVENT")) %>%
# mutate_at(vars(starts_with("EVENT")),
#           ~as.numeric(. == "TIME TO EVENT"))

# First time OAD
df_OADTM <- read_csv("data/ADaM/FINAL1.csv")
# df_OADINSTM <- read_csv("data/ADaM/FINAL2.csv")
Y_oad <- df_OADTM %>% 
  dplyr::select("USUBJID", 'PARAMCD', "AVAL", "PARCAT1") %>%
  dplyr::filter(PARAMCD %in% c("OADTM")) %>%
  rename("EVENT" = "PARCAT1") 

# Hyperglycaemia
# aval = AESTDTC - ?
Y_hyper <- df_adae %>% 
  dplyr::filter(AEDECOD %in% c("Hyperglycaemia") & AEONTRFL == "Y") %>% 
  dplyr::select("USUBJID", 'AEDECOD', 'AESTDTC')

# Hyperglycaemia
# aval = ASTDT - ?
Y_hypo <- df_hypo %>% 
  dplyr::filter(ACAT01  %in% c("SEVERE") & XHONTRFL == "Y") %>% 
  dplyr::select("USUBJID", 'ACAT01', 'ASTDT')

# Acid-base disorders
Y_keto <- df_adae %>% 
  dplyr::filter(AEHLGT == "Acid-base disorders" & AEONTRFL == "Y") %>%
  dplyr::select("USUBJID", 'AEDECOD', 'AESTDTC')



# cv outcomes
Y <- df_tte %>%
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

