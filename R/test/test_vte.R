rm(list = ls())
require(tidyverse)
library(tictoc)
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/simu/simu_dgd.R")) #Used for the current examples
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))


set.seed(123)
N <- 5e2 #size of generated data
df <- generate_data_simple(N, print_truth = TRUE)


# wrapper
tic()
res <- run_VTE(df = df,
                 sl_Q = sl_Q, 
                 sl_g = sl_g,
                 sl_x = sl_x,
                 ws = c('X2'), 
                 cv = F,
                 dr = T,
                 max_it = 1e4)
toc()
res_ee <- res$resEE
res_tmle <- res$resTMLE
res_ss <- res$resSS


# wrapper
tic()
res <- run_VIM2(df = df,
               sl_Q = sl_Q, 
               sl_g = sl_g,
               sl_x = sl_x,
               ws = c('X2'), 
               cv = F,
               dr = T,
               max_it = 1e4)
toc()
res_ee <- res$resEE
res_tmle <- res$resTMLE
res_ss <- res$resSS
