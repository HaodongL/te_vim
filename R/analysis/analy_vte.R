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
outcome = 'diab'; t = 24
df <- get_data(outcome, t, rm_baseIns=T)

# outcome = 'cv'; t = 24
# df <- get_data(outcome, t)

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data


### ------------  Part 2. import data  ------------ ###
res <- run_VTE(df = df, 
               sl_Q = sl_Q, 
               sl_g = sl_g,
               sl_x = sl_x,
               ws = c('statin_use'), 
               cv = T,
               dr = F,
               max_it = 1e4,
               lr = 1e-4)
res_ee <- res$resEE
res_tmle <- res$resTMLE
res_ss <- res$resSS

res_vte_t <- list('res_ee' = res_ee,
                  'res_tmle' = res_tmle,
                  'res_ss' = res_ss)

saveRDS(res_vte_t, file = "~/Repo/te_vim/data/res_vte_t_diab.RDS")

# df_fit <- fit_para(df = df,
#                    sl_Q = sl_Q, 
#                    sl_g = sl_g,
#                    sl_x = sl_x,
#                    ws = c('statin_use'), 
#                    dr = F)
# res_ee <- EE_VTE(df_fit)
# res_tmle <- TMLE_VTE(df_fit, max_it = 1e4, lr = 1e-4)


