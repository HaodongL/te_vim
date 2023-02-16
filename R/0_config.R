
library(here)
library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(ggplot2)
library(table1)
library(ggpubr)


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


# [1] "statin_use"        "statin_use"        "antihypertensives" "antihypertensives" "betab"            
# [6] "betab"             "minera"            "minera"            "adp"               "adp"              
# [11] "vkantag"           "vkantag"           "caantag"           "caantag"           "thiazide"         
# [16] "thiazide"          "loopdiur"          "loopdiur"   