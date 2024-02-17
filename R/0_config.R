
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
subgroup_levels <- c(`Overall` = "overall", 
                     `Statin use` = "statin_use", 
                     Antihypertensives = "antihypertensives",
                     `Beta blockers` = "betab",
                     `Mineralocorticoid receptor antagonists` = "minera", 
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



drug_levels = c(`Statin use` = "statin", 
                 `Anti-hypertensives` = "antihy",
                 `Beta blockers` = "betab",
                 `Mineralocorticoid receptor antagonists` = "minera", 
                 `ADP receptor inhibitors`= "adp",
                 `Vitamin K antagonists` = "vkanta",
                 `Ca antagonists` = "caanta", 
                 Thiazide = "thiazi",
                 `Loop diuretic` = "loopdi")


covar_levels = c(`Country` = "COUNTRY", `Age` = "AGE", `Sex` = "SEX", `Race` = "RACE", 
                 `Smoker status` = "SMOKER", `Diabetes duration` = "DIABDUR", 
                 `Antidiabetic therapy at baseline` = "ANTDBFL", `NYHA class (I - IV)` = "NYHACLAS", 
                 `Serum creatinine at baseline` = "CREATBL", `eGFR (Baseline) using MDRD (Calc)` = "EGFMDRBC", 
                 `Diabetic retinopathy severity` = "RETINSEV", `Myocardial infarction flag` = "MIFL", 
                 `Revascularization flag` = "REVASFL", `Carotid >50% stenosis on angiography` = "STENFL", 
                 `Coronary heart disease flag` = "CHDFL", `Ischaemic heart disease flag` = "IHDFL", 
                 `Chronic heart failure NYHA II-III flag` = "CHFFL", `Chronic kidney failure flag` = "KIDFL", 
                 `Microalbuminuria or proteinuria flag` = "MICFL", `Hypertension and LVH flag` = "HYPFL", 
                 `Left ventricular systolic and diastolic dysfunction flag` = "LVSDFL", 
                 `Ankle/Brachial index <0.9 flag` = "PADFL", `CV risk category` = "CVRISK", 
                 `HbA1c group (N)` = "HBA1CGRN", `Diabetes duration group (N)` = "DDURGRN", 
                 `Antihypertensive therapy` = "AHYPERFL", `Number of inclusion criterion passed` = "INCPASSN", 
                 `BMI at baseline` = "BMIBL", `Pulse at baseline` = "PULSEBL", 
                 `Systolic BP at baseline` = "SYSBPBL", `Diastolic BP at baseline` = "DIABPBL", 
                 `HbA1c at baseline (SI)` = "HBA1CBL", `HDL cholesterol at baseline (SI)` = "HDL1BL", 
                 `Calc. LDL cholesterol at baseline (SI)` = "LDL1BL", `Total cholesterol at baseline (SI)` = "CHOL1BL", 
                 `Triglycerides at baseline (SI)` = "TRIG1BL", `eGFR (Baseline) using MDRD (Calculated)` = "EGFMDRBC", 
                 `Peptic Ulcer and GERD Meds baseline flag` = "GERDBLFL", `Proton pump inhibitors flag` = "PPIFL", 
                 `H2 blockers flag` = "H2BLFL", `Stroke flag` = "STROKEFL",
                 `Statin use` = "statin_use", 
                 `Antihypertensives` = "antihypertensives",
                 `Beta blockers` = "betab",
                 `Remant cholesterol` ="RC",
                 `Mineralocorticoid receptor antagonists` = "minera", 
                 `ADP receptor inhibitors`= "adp",
                 `Vitamin K antagonists` = "vkantag",
                 `Ca antagonists` = "caantag", 
                 `Thiazide` = "thiazide",
                 `Loop diuretic` = "loopdiur")
