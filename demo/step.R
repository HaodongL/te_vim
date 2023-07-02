library(speff2trial)
library(dplyr)
library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)


rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/fit_hal.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/vim.R"))

### ------------  Part 1. Data ------------ ###
data(ACTG175)
df <- ACTG175
# df$cd420 <- df$cd420 - mean(df$cd420)
# df$cd420 <- df$cd420/sd(df$cd420)

df <- df %>% 
  filter(arms == 1 | arms == 3) %>% 
  select(c("arms", "cd420", "age", "wtkg", "karnof", "cd40", "cd80",
           "homo", "gender", "race", "symptom", "drugs", "hemo", "str2")) %>% 
  rename(A = "arms",
         Y = "cd420") %>% 
  mutate(A = as.numeric(A == 3))

### ------------  Part 2. Estimation ------------ ###
# vim loop over all covariates
set.seed(1)
ws = c("age", "wtkg", "karnof", "cd40", "cd80", "homo", 
       "gender", "race", "symptom", "drugs", "hemo", "str2")

n_ws <- length(ws)

df_fit <- fit_cvpara(df = df,
                   sl_Q = sl_Q, 
                   sl_g = sl_g,
                   sl_x = sl_x,
                   ws = ws[4],
                   dr = F,
                   Q_bounds = NULL,
                   g_bounds = c(0.025, 0.975),
                   tau_bounds = NULL,
                   tau_s_bounds = NULL,
                   gamma_s_bounds = NULL)

saveRDS(df_fit, file = "~/Repo/te_vim/data/step_df_fit.RDS")

resTMLE <- TMLE_VIM(df_fit, max_it=1e4, lr=0.0001)
resEE <- EE_VIM(df_fit)

# # Q
# task_Q <- sl3::make_sl3_Task(
#   data = df,
#   covariates = setdiff(names(df), c('Y')),
#   outcome = 'Y',
# )
# Q_fit <- sl_Q$train(task_Q)
# 
# # tau
# 
# # tau_s

df_fit_new <- resTMLE$df_fit_star

gamma_s = df_fit_new$gamma_s
tau_s = df_fit_new$tau_s
tau = df_fit_new$tau
QA = df_fit_new$mua_hat
Q1 = df_fit_new$mu1_hat
Q0 = df_fit_new$mu0_hat
po = df_fit_new$po
gn = df_fit_new$pi_hat

gamma_s_star = df_fit_new$gamma_s_star
tau_s_star = df_fit_new$tau_s_star
tau_star = df_fit_new$tau_star
QA_star = df_fit_new$mua_hat_star
A = df_fit_new$A
Y = df_fit_new$Y

N = nrow(df_fit_new)

# risk Q
mean((Y - QA)^2)
mean((Y - QA_star)^2)


# risk tau


# risk tau_s
mean((tau - tau_s)^2)
mean((tau - tau_s_star)^2)


### ------------  Part 3. Synthetic Data ------------ ###










