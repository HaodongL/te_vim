
rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))

### ------------  Part 1. import data  ------------ ###
outcome = 'diab2'; t = 24
df <- get_data(outcome, t, rm_baseIns=F)

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
df_fit <- fit_para(df, 
                   sl_Q, 
                   sl_g,
                   sl_x,
                   ws = list_ws,
                   dr = T)

# vte
res_vte_tmle <- TMLE_VTE(data = df_fit[[1]], max_it = 1e4, lr = 1e-4)
vte_n <- res_vte_tmle$coef
ic_vte <- res_vte_tmle$ic

# vim and scaled vim
df_res <- foreach(i = 1:length(list_ws), .combine = 'rbind') %do% {
  # vim
  res_vim_tmle <- TMLE_VIM(data = df_fit[[i]], max_it = 1e4, lr = 1e-4)
  vim_n <- res_vim_tmle$coef
  ic_vim <- res_vte_tmle$ic
  
  # scaled vim
  svim_n <- vim_n/vte_n
  ic_svim <- (ic_vim - svim_n*ic_vte)/vte_n
  ci_l_svim <- svim_n - 1.96*sd(ic_svim)/sqrt(nrow(df))
  ci_u_svim <- svim_n + 1.96*sd(ic_svim)/sqrt(nrow(df))
  
  
  df_res_i <- data.frame("varname" = rep(names(list_ws)[i], 2), 
                         "psiname" = c("VIM", "sVIM"),
                         "coef" = c(vim_n, svim_n), 
                         "ci_l" = c(res_vim_tmle$ci_l, ci_l_svim), 
                         "ci_u" = c(res_vim_tmle$ci_u, ci_u_svim))
  
  return(df_res_i)
}




