# upstream: data_clean.R
# Prepare data for formal analysis

library(dplyr)
library(readr)

get_data <- function(outcome, t = 24, rm_baseIns = FALSE, drop_censor = TRUE){
  
  assertthat::assert_that(outcome %in% c("diab", "diab2", "cv", 'a1c'))
  
  # read data based on outcome
  filepath <- paste0(repo_path,"/data/supp/df_all_", outcome, ".csv")
  df <- read_csv(filepath)
  
  # drop those censored before t
  df_C <- read_csv(paste0(repo_path,"/data/supp/df_censor.csv"))
  if (drop_censor){
    n1 = nrow(df)
    df <- df %>% left_join(df_C, by='USUBJID') %>% filter(AVAL > t) %>% select(-AVAL)
    n2 = nrow(df)
    cat(paste0(n1-n2, " people censored before ", t, " months are dropped."))
  }
  
  # convert character/logical to factor
  df <- df %>% 
    mutate(across(where(is.character), ~ as.factor(.))) %>% 
    mutate(across(where(is.logical), ~ as.factor(.))) 
  
  # remove INSNVFL == TRUE (confirm wtih novo)
  if (rm_baseIns){
    df <- df %>% filter(INSNVFL == FALSE) %>% select(-INSNVFL)
  }
  
  # define Y
  if (outcome %in% c('diab', 'diab2')){
    df <- df %>% mutate(Y = ifelse((AVAL_FRINSLTM <= t & EVENT_FRINSLTM == 1 )|
                                   (AVAL_OADTM <= t & EVENT_OADTM == 1 )|
                                   (AVAL_HYPER <= t & EVENT_HYPER == 1 )|
                                   (AVAL_HYPO <= t & EVENT_HYPO == 1 )|
                                   (AVAL_ACIDOSIS <= t & EVENT_ACIDOSIS == 1),
                                   1, 0)) %>%
      select(-c("USUBJID", "AVAL_FRINSLTM",
                "AVAL_OADTM", "AVAL_HYPER",
                "AVAL_HYPO", "AVAL_ACIDOSIS",
                "EVENT_FRINSLTM", "EVENT_OADTM",
                "EVENT_HYPER", "EVENT_HYPO", "EVENT_ACIDOSIS"))
  }else if (outcome == "cv"){
    df <- df %>% mutate(Y = ifelse((AVAL_MACE <= t & EVENT_MACE == 1 )|
                                   (AVAL_NFMI <= t & EVENT_NFMI == 1 )|
                                   (AVAL_CVDEATH <= t & EVENT_CVDEATH == 1 )|
                                   (AVAL_NFSTROKE <= t & EVENT_NFSTROKE == 1 )|
                                   (AVAL_NONCVDEATH <= t & EVENT_NONCVDEATH == 1),
                                   1, 0)) %>%
                 select(-c("USUBJID", "AVAL_MACE",
                           "AVAL_NFMI", "AVAL_CVDEATH",
                           "AVAL_NFSTROKE", "AVAL_NONCVDEATH",
                           "EVENT_MACE", "EVENT_NFMI",
                           "EVENT_CVDEATH", "EVENT_NFSTROKE", "EVENT_NONCVDEATH"))
  }else {
    if (t == 24){
      df <- df %>% mutate(Y = HBA1C_VISIT_9 -  HBA1C_VISIT_3)
    }else if (t == 12){
      df <- df %>% mutate(Y = HBA1C_VISIT_7 -  HBA1C_VISIT_3)
    }else if (t == 36){
      df <- df %>% mutate(Y = HBA1C_VISIT_11 -  HBA1C_VISIT_3)
    }else if (t == 48){
      df <- df %>% mutate(Y = HBA1C_VISIT_13 -  HBA1C_VISIT_3)
    }
    df <- df %>% select(-c("USUBJID", "HBA1C_VISIT_3", "HBA1C_VISIT_5",  
                      "HBA1C_VISIT_6", "HBA1C_VISIT_7", "HBA1C_VISIT_8",     
                      "HBA1C_VISIT_9", "HBA1C_VISIT_10", "HBA1C_VISIT_11",    
                      "HBA1C_VISIT_12", "HBA1C_VISIT_13"))
  }
  
  return(df)
}



get_data_tte <- function(outcome){
  
  assertthat::assert_that(outcome %in% c("diab", "cv"))
  
  # read data based on outcome
  filepath <- paste0(repo_path,"/data/supp/df_all_", outcome, ".csv")
  df <- read_csv(filepath)
  
  # convert character/logical to factor
  df <- df %>% 
    mutate(across(where(is.character), ~ as.factor(.))) %>% 
    mutate(across(where(is.logical), ~ as.factor(.))) 
  
  # remove INSNVFL == TRUE (confirm wtih novo)
  df <- df %>% filter(INSNVFL == FALSE) %>% select(-INSNVFL)
  
  # remove highly correlated cols (already removed from data_clean)
  # var_corr <- c("EGFRMDRC", "EGFREPB", "EGFMDRBC", "EGFRMDR")
  # df <- df %>% select(-all_of(var_corr))
  
  # define Y
  if (outcome == "diab"){
    df <- df %>% mutate(Y = pmin(AVAL_FRINSLTM, AVAL_OADTM, AVAL_HYPER,
                                AVAL_HYPO, AVAL_ACIDOSIS),
                        delta = as.numeric(psum(EVENT_FRINSLTM, EVENT_OADTM,
                                               EVENT_HYPER, EVENT_HYPO,
                                               EVENT_ACIDOSIS) > 0)) %>%
      select(-c("USUBJID", "AVAL_FRINSLTM",
                "AVAL_OADTM", "AVAL_HYPER",
                "AVAL_HYPO", "AVAL_ACIDOSIS",
                "EVENT_FRINSLTM", "EVENT_OADTM",
                "EVENT_HYPER", "EVENT_HYPO", "EVENT_ACIDOSIS"))
  }else if (outcome == "cv"){
    df <- df %>% mutate(Y = pmin(AVAL_MACE, AVAL_NFMI, AVAL_CVDEATH,
                                AVAL_NFSTROKE, AVAL_NONCVDEATH),
                        delta = as.numeric(psum(EVENT_MACE, EVENT_NFMI,
                                               EVENT_CVDEATH, EVENT_NFSTROKE,
                                               EVENT_NONCVDEATH) > 0)) %>%
      select(-c("USUBJID", "AVAL_MACE",
                "AVAL_NFMI", "AVAL_CVDEATH",
                "AVAL_NFSTROKE", "AVAL_NONCVDEATH",
                "EVENT_MACE", "EVENT_NFMI",
                "EVENT_CVDEATH", "EVENT_NFSTROKE", "EVENT_NONCVDEATH"))
    }
  return(df)
}

psum <- function(...,na.rm=FALSE) { 
  rowSums(do.call(cbind,list(...)),na.rm=na.rm) 
  } 

