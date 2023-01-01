rm(list = ls())
require(tidyverse)
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/simu/simu_dgd.R")) #Used for the current examples
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))
source(paste0(repo_path, "R/est_function/fit_hal.R"))
library(tictoc)

set.seed(123)
N <- 5e2 #size of generated data
# df <- generate_data_v2(N, print_truth = TRUE) #VIM_Theta_s: 1.929, 0.11
df <- generate_data_simple(N, print_truth = TRUE) 

# ee, tmle
tic()
res <- run_VIM_Theta(df = df,
                     sl_Q = sl_Q, 
                     sl_g = sl_g,
                     sl_x = sl_x,
                     ws = c('X2'), 
                     cv = F,
                     dr = F,
                     lfm_linear = F, 
                     tmle_dr_update = F, 
                     max.it = 1e4, 
                     Q_bounds = c(0.001, 0.999), 
                     g_bounds = c(0.025, 0.975),
                     tau_bounds = c(-1+1e-3, 1-1e-3),
                     tau_s_bounds = c(-1+1e-3, 1-1e-3),
                     gamma_s_bounds = c(1e-6, 1-1e-6)
)
toc()
res_ee <- res$resEE
res_tmle <- res$resTMLE


# test HAL
tic()
df_fit_hal <- fit_para_hal(df,
                           sl_Q_hal,
                           sl_g_hal,
                           sl_ws_hal,
                           ws = c('X2'))

res_hal <- HAL_VIM(df_fit_hal)
toc()