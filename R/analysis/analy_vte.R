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
df <- get_data(outcome, t)

# outcome = 'cv'; t = 24
# df <- get_data(outcome, t)

# outcome = 'a1c'; t = 24
# df <- get_data(outcome, t)


nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data

sl_Q <- lrnr_earth
sl_g <- lrnr_earth
sl_x <- lrnr_earth

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

# saveRDS(res_vte_t, file = "~/Repo/te_vim/data/res_vte_t_cv_diab.RDS")

# df_fit <- fit_para(df = df,
#                    sl_Q = sl_Q, 
#                    sl_g = sl_g,
#                    sl_x = sl_x,
#                    ws = c('statin_use'), 
#                    dr = F)
# res_ee <- EE_VTE(df_fit)
# res_tmle <- TMLE_VTE(df_fit, max_it = 1e4, lr = 1e-4)


#  oh
rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source("~/Repo/TE-Heterogeneity/R/vte.R")
source("~/Repo/TE-Heterogeneity/R/aipw.R")
source("~/Repo/TE-Heterogeneity/R/tmle.R")
source("~/Repo/TE-Heterogeneity/R/example_helpers.R")



### ------------  Part 1. import data  ------------ ###
outcome = 'diab'; t = 24
df <- get_data(outcome, t, rm_baseIns=T)

# outcome = 'cv'; t = 24
# df <- get_data(outcome, t)

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data

sl_Q <- lrnr_xgb
sl_g <- lrnr_xgb
sl_x <- lrnr_xgb

df_fit <- fit_para(df, 
                     sl_Q, 
                     sl_g,
                     sl_x,
                     ws = c("statin"))

aipw_cv   <- VTE(df_fit, method="AIPW")
tmle_cv   <- VTE(df_fit, method="TMLE")


# test 
sl_Q <- lrnr_xgb
sl_g <- lrnr_xgb
sl_x <- lrnr_xgb

df_fit <- fit_cvpara(df, 
                     sl_Q, 
                     sl_g,
                     sl_x,
                     ws = c("SMOKER"))

res_vim <- TMLE_VIM(data = df_fit, max_it = 1e4, lr = 1e-4)
res_vim$coef

res_vte <- TMLE_VTE(data = df_fit, max_it = 1e4, lr = 1e-4)
res_vte$coef

res_vim_a <- TMLE_VIM_a(data = df_fit, y_l = 0, y_u = 1, max_it = 1e4, lr = 1e-4)
res_vim_a$coef


EE_VTE(data = df_fit)$coef
EE_VIM(data = df_fit)$coef

SS_VTE(data = df_fit)$coef
SS_VIM(data = df_fit)$coef

TEVIM_VTE(data = df_fit, marginal_cates = df_fit$tau_s)$coef

res_both <- TMLE_VIMnVTE(data = df_fit, max_it = 1e4, lr = 1e-4)




var(df_fit$mu1_hat - df_fit$mu0_hat)

ate = mean(df_fit$mu1_hat - df_fit$mu0_hat)
mean((df_fit$tau - ate)^2)


# Estimator: TMLE
data = df_fit
max_it = 1e4
lr = 1e-4

N <- nrow(data)
A <- data$A
Y <- data$Y

# Q, g, tau, tau_s, gamma_s
QA_0 <- data$mua_hat
Q1_0 <- data$mu1_hat
Q0_0 <- data$mu0_hat
gn <- data$pi_hat

tau_0 <- data$tau
ate_0 <- rep(mean(tau_0), N)

# sig1
sig1 <- sd(2*(tau_0 - ate_0)*(A/gn - (1-A)/(1-gn))*(Y - QA_0))
# eps1
eps1 <- lr

# i
i <- 1
QA_i <- QA_0
Q1_i <- Q1_0
Q0_i <- Q0_0
tau_i <- tau_0
ate_i <- ate_0

# logging <- matrix(NA, max_it, 3)

while(i <= max_it){
  # step 1. update Q and tau
  if (i == 1){
    HQ_1i <- 2*(tau_i - ate_i)
    H1_1i <- 2*(tau_i - ate_i)*(1/gn)
    H0_1i <- 2*(tau_i - ate_i)*(-1/(1-gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i)
    D_1i <- HA_1i*(Y - QA_i)
    PnD_1i <- mean(D_1i)
  }
  
  Q1_i <- Q1_i + eps1*(HQ_1i)*sign(PnD_1i)
  Q0_i <- Q0_i + eps1*(-HQ_1i)*sign(PnD_1i)
  QA_i <- ifelse(A, Q1_i, Q0_i)
  
  tau_i <- Q1_i - Q0_i
  ate_i <- rep(mean(tau_i), N)
  
  # update D1 D2, check criteria
  HQ_1i <- 2*(tau_i - ate_i)
  H1_1i <- 2*(tau_i - ate_i)*(1/gn)
  H0_1i <- 2*(tau_i - ate_i)*(-1/(1-gn))
  HA_1i <- ifelse(A, H1_1i, H0_1i)
  D_1i <- HA_1i*(Y - QA_i)
  PnD_1i <- mean(D_1i)
  
  c1 <- abs(PnD_1i) <= sig1/(sqrt(N)*log(N))
  
  if (c1){
    break
  }
  i <- i + 1
}

QA_star <- QA_i
Q1_star <- Q1_i
Q0_star <- Q0_i
tau_star <- tau_i
ate_star <- rep(mean(tau_i), N)
gamma_s_star <- rep(mean(tau_i^2), N)

if(i>=max_it) warning("Max iterations reached in TMLE")

# calculate Theta_s, scale back
# theta_s_star <- mean(gamma_s_star - ate_star^2)
theta_s_star <- mean((tau_star - ate_star)^2)

mean(tau_star^2) -  mean(tau_star)^2 

# ic <- 2*(tau_star - ate_star)*(A/gn - (1-A)/(1-gn))*(Y-QA_star) + (tau_star - ate_star)^2 - theta_s_star
# use initial est for ic
theta_s_0 <- mean((tau_0 - ate_0)^2)
ic <- 2*(tau_0 - ate_0)*(A/gn - (1-A)/(1-gn))*(Y-QA_0) + (tau_0 - ate_0)^2 - theta_s_0
se <- sqrt(var(ic)/N)
  

par(mfrow = c(1,2))
hist(QA_star)
hist(QA_0)


par(mfrow = c(1,2))
hist(tau_star)
hist(tau_0)

ate_star[1]
ate_0

mean(tau_star)
mean(tau_0)

var(tau_star)
var(tau_0)

mean(Q1_star)
mean(Q1_0)


mean(Q0_star)
mean(Q0_0)

mean(QA_star)
mean(QA_0)
