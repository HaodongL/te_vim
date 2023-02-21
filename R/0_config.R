
library(here)
library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(ggplot2)
library(tidyverse)
library(table1)
library(ggpubr)
library(kableExtra)
library(readr)
library(R6)

if(here::here()=="C:/Users/andre/Documents/jici/te_vim"){
  
  repo_path = "C:/Users/andre/Documents/jici/te_vim/"
  
}else{
  repo_path = "~/Repo/te_vim/"
}


#Variable labels

# $W
# [1] "AGE"               "SEX"               "RACE"              "COUNTRY"           "SMOKER"            "NYHACLAS"          "DIABDUR"           "ANTDBFL"           "AHYPERFL"         
# [10] "INCPASSN"          "BMIBL"             "PULSEBL"           "SYSBPBL"           "DIABPBL"           "HBA1CBL"           "HDL1BL"            "LDL1BL"            "CHOL1BL"          
# [19] "RC"                "RCoverHDL"         "TRIG1BL"           "CREATBL"           "EGFMDRBC"          "RETINSEV"          "GERDBLFL"          "PPIFL"             "H2BLFL"           
# [28] "MIFL"              "STROKEFL"          "REVASFL"           "STENFL"            "CHDFL"             "IHDFL"             "CHFFL"             "KIDFL"             "MICFL"            
# [37] "HYPFL"             "LVSDFL"            "PADFL"             "CVRISK"            "HBA1CGRN"          "DDURGRN"           "statin_use"        "antihypertensives" "betab"            
# [46] "minera"            "adp"               "vkantag"           "caantag"           "thiazide"          "loopdiur"         

# covariate_levels <- c(age = "AGE", sex = "SEX",
#                       race= = "RACE"
#                       country,= "COUNTRY", 
#                       smoker=   "SMOKER", 
#                       `NYHACLAS`=    "NYHACLAS", 
#                       `Diabetes duration`=  "DIABDUR",
#                       `Antidiabetic therapy at BL`=   "ANTDBFL", 
#                       =   "AHYPERFL", 
#                       = "INCPASSN",
#                       =  "BMIBL", 
#                       =     "PULSEBL", 
#                       =   "SYSBPBL", 
#                       =   "DIABPBL", 
#                       =   "HBA1CBL", 
#                       =   "HDL1BL", 
#                       =    "LDL1BL", 
#                       =    "CHOL1BL",
#                       =   "RC",
#                       =   "RCoverHDL", 
#                       = "TRIG1BL",
#                       =   "CREATBL", 
#                       =   "EGFMDRBC", 
#                       =  "RETINSEV", 
#                       =  "GERDBLFL", 
#                       =  "PPIFL", 
#                       =     "H2BLFL", 
#                       =     "MIFL",
#                       = "STROKEFL",
#                       =  "REVASFL", 
#                       =   "STENFL", 
#                       =    "CHDFL", 
#                       =     "IHDFL",
#                       =     "CHFFL", 
#                       =     "KIDFL", 
#                       =     "MICFL", 
#                       =     "HYPFL", 
#                       =     "LVSDFL", 
#                       =    "PADFL", 
#                       =     "CVRISK", 
#                       =    "HBA1CGRN", 
#                       =  "DDURGRN", 
#   `Statin use` = "statin_use", 
# Antihypertensives = "antihypertensives",
# `Beta blockers` = "betab",
# `Mineralocorticoid receptor	antagonists` = "minera", 
# ADP = "adp",
# `Vitamin K antagonists` = "vkantag",
# `Ca antagonists` = "caantag", 
# Thiazide = "thiazide",
# `Loop diuretic` = "loopdiur")

# [1] "statin_use"   "" "betab"            
# [6] "betab"             ""            "minera"            "adp"               "adp"              
# [11] "vkantag"           "vkantag"           "caantag"           "caantag"           "thiazide"         
# [16] "thiazide"          "loopdiur"          "loopdiur"   

subgroup_levels <- c(`Statin use` = "statin_use", 
                     Antihypertensives = "antihypertensives",
                     `Beta blockers` = "betab",
                     `Mineralocorticoid receptor	antagonists` = "minera", 
                     ADP = "adp",
                     `Vitamin K antagonists` = "vkantag",
                     `Ca antagonists` = "caantag", 
                     Thiazide = "thiazide",
                    `Loop diuretic` = "loopdiur")
