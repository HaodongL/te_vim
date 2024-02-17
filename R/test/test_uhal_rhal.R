rm(list = ls())
require(tidyverse)
library(hal9001)
library(glmnet)
library(R6)
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/est_function/relax_hal.R"))
source(paste0(repo_path, "R/est_function/under_hal.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) #Used for the current examples
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))
library(tictoc)

set.seed(123)
N <- 5e2 #size of generated data
# df <- generate_data_v2(N, print_truth = TRUE) #VIM_Theta_s: 1.929, 0.11
df <- generate_data_simple(N, print_truth = TRUE) 


lrnr_uhal <- Lrnr_uhal9001$new()
lrnr_rhal <- Lrnr_rhal9001$new()

sl_Q <- lrnr_uhal
sl_g <- lrnr_earth
sl_x <- lrnr_earth


# ee, tmle
tic()
res <- run_VIM(df = df,
               sl_Q = sl_Q, 
               sl_g = sl_g,
               sl_x = sl_x,
               ws = c('X2'), 
               cv = F,
               dr = F
)
toc()
res_ss <- res$resSS
res_ee <- res$resEE
res_tmle <- res$resTMLE