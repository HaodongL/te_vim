
rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))

### ------------  Part 1. import data  ------------ ###
# outcome = 'diab2'; t = 24
# df <- get_data(outcome, t, rm_baseIns=F)

outcome = 'diab'; t = 24
df <- get_data(outcome, t, rm_baseIns=T)

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data


### ------------  Part 2. Estimation ------------ ###
# vim loop over all covariates
set.seed(3)
list_ws <- list(
  "Country" = c("COUNTRY"),
  "Age" = c("AGE"), 
  "Sex" = c("SEX"), 
  "Race" = c("RACE"),
  "Smoker status" = c("SMOKER"),
  "Diabetes duration" = c("DIABDUR"),
  "Antidiabetic therapy at baseline" = c("ANTDBFL"),
  "NYHA class (I - IV)" = c("NYHACLAS"),
  "Kidney function" = c("CREATBL", "EGFMDRBC"),
  "Diabetic retinopathy severity" = c("RETINSEV"),
  "Previous CVD" = c("MIFL", "REVASFL", "STENFL", "CHDFL", "IHDFL", "HYPFL", "PADFL", "STROKEFL"),
  "Heart failure" = c("CHFFL", "LVSDFL", "CVRISK"),
  "Kidney disease" = c("KIDFL", "MICFL"),
  "BP medication" = c("AHYPERFL", "antihypertensives", "betab", "minera", "caantag", "thiazide", "loopdiur"),
  "Number of inclusion criterion passed" = c("INCPASSN"),
  "BMI at baseline" = c("BMIBL"),
  "Pulse at baseline" = c("PULSEBL"),
  "Blood pressure" = c("SYSBPBL", "DIABPBL"),
  "Lipids" = c("HDL1BL", "LDL1BL", "CHOL1BL", "TRIG1BL", "statin_use", "RC", "RCoverHDL"),
  "Peptic Ulcer and GERD Meds" = c("GERDBLFL", "PPIFL", "H2BLFL"),
  "HbA1c at Baseline (SI)" = c("HBA1CBL"),
  "ADP" = c("adp"),
  "Vitamin K antagonists" = c("vkantag"),
  "Insulin naive" = c("INSNVFL")
)

# fit para
df_fit <- fit_cvpara(df, 
                     sl_Q, 
                     sl_g,
                     sl_x,
                     ws = list_ws)

save_path <- paste0("~/Repo/te_vim/analy_res/df_fit_", outcome,".RDS")
saveRDS(df_fit, file = save_path)
