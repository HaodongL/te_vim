library(here)
library(tidyverse)
library(dplyr)
library(mice)


# read in subject level covariate dataset, tte, ttse and cm 
df_w <- haven::read_sas("data/ADaM/adsl.sas7bdat")
df_tte <- haven::read_sas("data/ADaM/adtte.sas7bdat")
df_cm <- haven::read_sas("data/ADaM/adcm.sas7bdat")
df_ttse <- haven::read_sas("data/ADaM/adttse.sas7bdat")
df_hypo <- haven::read_sas("data/ADaM/adhypo.sas7bdat")
df_adae <- haven::read_sas("data/ADaM/adae.sas7bdat")
# df_adlb <- haven::read_sas("data/ADaM/adlb.sas7bdat")

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

# Collinearity? - 4 renal failure columns (RENFSEV, RENFCKD, RNFSCSEV, RNFSCCKD)
# dplyr::select(W, RENFSEV, RENFCKD, RNFSCSEV, RNFSCCKD) %>% group_by_all() %>% count() %>% view()

# excessive missingness (ISMHIBFL [99.3%], ISMHGPFL [99.9%], RETINSEV [77.9%], NYHACLAS[82.3%])
# skimr::skim(dplyr::select(W, ISMHIBFL, ISMHGPFL, RETINSEV, NYHACLAS))

# remove vars with >99% missingness
# W <- dplyr::select(W, -c(ISMHIBFL, ISMHGPFL)) # ??? doesn't exist.

# how to deal with missing retinopathy severity & NYHA class?
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
         ARM = as.numeric(ARM == "Liraglutide")# ,
         #ETHNIC = ETHNIC == "HISPANIC OR LATINO"
  ) # %>%
# rename(MALE = SEX,
#        HISPANIC = ETHNIC)

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

CM <- df_cm %>% 
  select(USUBJID, CMCLASCD, CONBLFL) %>% 
  filter(CMCLASCD == 'C10AA') %>% 
  group_by(USUBJID) %>% 
  summarise(
    statin_use = TRUE,
    statin_cm = sum(CONBLFL == 'Y') > 0) 

W <- left_join(W, CM, by = 'USUBJID') 

W <- W %>% mutate(statin_use = ifelse(is.na(statin_use), FALSE, statin_use),
                  statin_cm = ifelse(is.na(statin_cm), FALSE, statin_cm))


# identify disease progression outcome

Y_secondary <- df_ttse %>%
  dplyr::select("USUBJID", 'PARAMCD', "AVAL", "PARCAT1") %>%
  dplyr::filter(PARAMCD %in% c("FRINSLTM")) %>%
  rename("EVENT" = "PARCAT1") 
  # pivot_wider(id_cols = "USUBJID", names_from = "PARAMCD", values_from = c("AVAL", "EVENT")) %>%
  # mutate_at(vars(starts_with("EVENT")),
  #           ~as.numeric(. == "TIME TO EVENT"))

# First time OAD


df_OADTM <- read_csv("data/ADaM/FINAL1.csv")
df_OADINSTM <- read_csv("data/ADaM/FINAL2.csv")

# tmp <- df_cm %>% 
#   filter(ASTRF %in% c("DURING","AFTER"), FASFL == "Y")
# 
# drop_sujid_oad <- tmp %>% 
#   filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
#             any_vars(str_detect(., "A10BX"))) %>% 
#   filter(CMDECOD == 'LIRAGLUTIDE') 
# 
# df_oad <- tmp %>% 
#   filter_at(vars(matches("CMCLASCD|DCL[0-9]*C")),
#             any_vars(str_detect(., "A10BA|A10BB|A10BC|A10BD|A10BF|
#                                 A10BG|A10BH|A10BX09|A10BX11|A10BX12|A10BX"))) 

  


# 'OADTM' %in% unique(df_ttse$PARAMCD)

# outcomes
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

