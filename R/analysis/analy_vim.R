library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(foreach)
library(doParallel)
# library(ggplot2)
# library(table1)
# library(ggpubr)

rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/fit_hal.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/vim.R"))

### ------------  Part 1. import data  ------------ ###
# outcome = 'diab'; t = 24
# df <- get_data(outcome, t, rm_baseIns=T)

outcome = 'diab2'; t = 24
df <- get_data(outcome, t, rm_baseIns=F)

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data


### ------------  Part 2. Estimation ------------ ###
# vim loop over all covariates
set.seed(1)
ws = c("AGE", "SEX", "RACE", "COUNTRY", "SMOKER", "NYHACLAS",
       "DIABDUR", "ANTDBFL", "AHYPERFL", "INCPASSN", 
       "BMIBL", "PULSEBL", "SYSBPBL", "DIABPBL", "HBA1CBL", "HDL1BL", "LDL1BL",
       "CHOL1BL", "RC", "RCoverHDL","TRIG1BL", "CREATBL", "EGFMDRBC", 
       "RETINSEV", "GERDBLFL", "PPIFL", "H2BLFL", 
       "MIFL", "STROKEFL", "REVASFL", "STENFL", "CHDFL", "IHDFL", "CHFFL",
       "KIDFL", "MICFL", "HYPFL", "LVSDFL", "PADFL", "CVRISK", "HBA1CGRN", "DDURGRN", 
       "statin_use", "antihypertensives", "betab", "minera", "adp",
       "vkantag", "caantag", "thiazide", "loopdiur")
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

# ws = c("AGE", "SEX")
n_ws <- length(ws)

tic()
registerDoParallel(9)
df_vim <- foreach(i = 1:n_ws, .combine = 'rbind') %do% {
  print(paste0("Covar name: ", ws[i]))
  res <- run_VIM(df = df,
                 sl_Q = sl_Q, 
                 sl_g = sl_g,
                 sl_x = sl_x,
                 ws = ws[i], 
                 cv = T,
                 dr = F,
                 max_it = 1e4)
  res_ee <- res$resEE
  res_tmle <- res$resTMLE
  res_ss <- res$resSS
  
  df_vim_i <- data.frame("varname" = c(ws[i], ws[i], ws[i]), 
                         "importance" = c(res_ee$coef, res_tmle$coef, res_ss$coef), 
                         "ci_l" = c(res_ee$ci_l, res_tmle$ci_l, res_ss$ci_l), 
                         "ci_u" = c(res_ee$ci_u, res_tmle$ci_u, res_ss$ci_u), 
                         "method" = c('EE', 'TMLE', 'SS'))
  
  # if (ws[i] %in% cm_names){
  #   df_fit_hal <- fit_para_hal(df = df, 
  #                              sl_Q = sl_Q_hal, 
  #                              sl_g = sl_g_hal,
  #                              sl_ws = sl_ws_hal,
  #                              ws = ws[i])
  #   res_hal <- HAL_VIM(df_fit_hal)
  #   df_vim_i2 <- data.frame("varname" = ws[i], 
  #                          "importance" = res_hal$coef, 
  #                          "ci_l" = res_hal$ci_l, 
  #                          "ci_u" = res_hal$ci_u, 
  #                          "method" = 'SS-HAL')
  #   df_vim_i <- rbind(df_vim_i, df_vim_i2)
  # }
  return(df_vim_i)
}
toc()

saveRDS(df_vim, file = "~/Repo/te_vim/data/df_vim_dr_diab2.RDS")