
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
subgroup_levels <- c(`Statin use` = "statin_use", 
                     Antihypertensives = "antihypertensives",
                     `Beta blockers` = "betab",
                     `Mineralocorticoid receptor	antagonists` = "minera", 
                     ADP = "adp",
                     `Vitamin K antagonists` = "vkantag",
                     `Ca antagonists` = "caantag", 
                     Thiazide = "thiazide",
                    `Loop diuretic` = "loopdiur")


# #fix labels
# library(rvest)
# library(xml2)
# file<-read_html(paste0(repo_path, "data/Documents/define.html"))
# tables<-html_nodes(file, "table")
# var_tab <- NULL
# for(i in 2:length(tables)){
#   temp_tab <- html_table(tables[i], fill = TRUE)[[1]]
#   temp_tab$`Length or Display Format` <- as.numeric(temp_tab$`Length or Display Format`)
#   var_tab <- bind_rows(var_tab, temp_tab)
# }
# 
# 
# ws = c("AGE", "SEX", "RACE", "COUNTRY", "SMOKER", "NYHACLAS",
#        "DIABDUR", "ANTDBFL", "AHYPERFL", "INCPASSN", 
#        "BMIBL", "PULSEBL", "SYSBPBL", "DIABPBL", "HBA1CBL", "HDL1BL", "LDL1BL",
#        "CHOL1BL", "RC", "RCoverHDL","TRIG1BL", "CREATBL", "EGFMDRBC", 
#        "RETINSEV", "GERDBLFL", "PPIFL", "H2BLFL", 
#        "MIFL", "STROKEFL", "REVASFL", "STENFL", "CHDFL", "IHDFL", "CHFFL",
#        "KIDFL", "MICFL", "HYPFL", "LVSDFL", "PADFL", "CVRISK", "HBA1CGRN", "DDURGRN", 
#        "statin_use", "antihypertensives", "betab", "minera", "adp",
#        "vkantag", "caantag", "thiazide", "loopdiur")
# 
# 
# 
# var_tab <- var_tab[var_tab$Variable %in% ws,c(1,2)] %>% distinct()
# var_tab
# 
# var_labs <- var_tab$Variable
# names(var_labs) <- var_tab$`Label / Description`
# dput(var_labs)


covar_levels = c(Country = "COUNTRY", Age = "AGE", Sex = "SEX", Race = "RACE", 
  `Smoker status` = "SMOKER", `Diabetes Duration` = "DIABDUR", 
  `Antidiabetic therapy at Baseline` = "ANTDBFL", `NYHA CLASS (I - IV)` = "NYHACLAS", 
  `Serum Creatinine at Baseline` = "CREATBL", `eGFR (Baseline) using MDRD (Calc)` = "EGFMDRBC", 
  `Diabetic Retinopathy Severity` = "RETINSEV", `Myocardial Infarction Flag` = "MIFL", 
  `Revascularization Flag` = "REVASFL", `Carotid >50% stenosis on angiography` = "STENFL", 
  `Coronary Heart Disease Flag` = "CHDFL", `Ischaemic Heart Disease Flag` = "IHDFL", 
  `Chronic Heart Failure NYHA II-III Flag` = "CHFFL", `Chronic Kidney Failure Flag` = "KIDFL", 
  `Microalbuminuria or Proteinuria Flag` = "MICFL", `Hypertension and LVH Flag` = "HYPFL", 
  `Left Ventricular Systolic and Diastolic Dysfunction Flag` = "LVSDFL", 
  `Ankle/Brachial index <0.9 Flag` = "PADFL", `CV Risk Category` = "CVRISK", 
  `HbA1c Group (N)` = "HBA1CGRN", `Diabetes Duration Group (N)` = "DDURGRN", 
  `Antihypertensive therapy` = "AHYPERFL", `Number of Inclusion Criterion Passed` = "INCPASSN", 
  `Body Mass Index at Baseline` = "BMIBL", `Pulse at Baseline` = "PULSEBL", 
  `Systolic BP at Baseline` = "SYSBPBL", `Diastolic BP at Baseline` = "DIABPBL", 
  `HbA1c at Baseline (SI)` = "HBA1CBL", `HDL Cholesterol at Baseline (SI)` = "HDL1BL", 
  `Calc. LDL Cholesterol at Baseline (SI)` = "LDL1BL", `Total Cholesterol at Baseline (SI)` = "CHOL1BL", 
  `Triglycerides at Baseline (SI)` = "TRIG1BL", `eGFR (Baseline) using MDRD (Calculated)` = "EGFMDRBC", 
  `Peptic Ulcer and GERD Meds Baseline Flag` = "GERDBLFL", `Proton Pump Inhibitors Flag` = "PPIFL", 
  `H2 Blockers Flag` = "H2BLFL", `Stroke Flag` = "STROKEFL",
  `Statin use` = "statin_use", 
  Antihypertensives = "antihypertensives",
  `Beta blockers` = "betab",
  `Mineralocorticoid receptor	antagonists` = "minera", 
  ADP = "adp",
  `Vitamin K antagonists` = "vkantag",
  `Ca antagonists` = "caantag", 
  Thiazide = "thiazide",
  `Loop diuretic` = "loopdiur")
