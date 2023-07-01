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

data(ACTG175)
df <- ACTG175
# df$cd420 <- df$cd420 - mean(df$cd420)
# df$cd420 <- df$cd420/sd(df$cd420)

df_sub <- df %>% 
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

# ws = c("gender")

# ws = c("AGE", "SEX")
n_ws <- length(ws)



# One
# res <- run_VIM(df = df_sub,
#                sl_Q = sl_Q, 
#                sl_g = sl_g,
#                sl_x = sl_x,
#                ws = ws[2], 
#                cv = F,
#                dr = F,
#                max_it = 1e4)
# res_ee <- res$resEE
# res_tmle <- res$resTMLE


df_fit <- fit_para(df = df_sub,
                   sl_Q = sl_Q, 
                   sl_g = sl_g,
                   sl_x = sl_x,
                   ws = ws[2],
                   dr = F,
                   Q_bounds = NULL,
                   g_bounds = c(0.025, 0.975),
                   tau_bounds = NULL,
                   tau_s_bounds = NULL,
                   gamma_s_bounds = NULL)

resTMLE <- TMLE_VIM(df_fit, max_it=1e4, lr=0.0001)
resEE <- EE_VIM(df_fit)
# resSS <- SS_VIM(df_fit)



df_fit_new <- resTMLE$df_fit_star

plot_compare <- function(x,y){
  par(mfrow = c(1,2))
  hist(x)
  hist(y)
}

# Q
plot_compare(df_fit_new$mua_hat, df_fit_new$mua_hat_star)

# tau
plot_compare(df_fit_new$tau, df_fit_new$tau_star)

# tau_s
plot_compare(df_fit_new$tau_s, df_fit_new$tau_s_star)

# gamma_s
plot_compare(df_fit_new$gamma_s, 
             df_fit_new$gamma_s_star)

# df_fit_new$gamma_s[which(df_fit_new$gamma_s>5000)] = 5000
# df_fit_new$gamma_s_star[which(df_fit_new$gamma_s_star>5000)] = 5000

plot_compare(df_fit_new$tau, df_fit_new$tau_s)

plot_compare(df_fit_new$tau_star, df_fit_new$tau_s_star)


# tau fit
# DR
po = df_fit_new$po
tau_fit <- fit_x(df = df_sub, sl_x = sl_x, po = po, outcome = 'po', para = 'tau')
tau <- tau_fit$predict()

# T
tau <- df_fit_new$mu1_hat - df_fit_new$mu0_hat 


# tau_s fit 
tau_s_fit <- fit_x(df = df_sub, sl_x = sl_x, tau = tau, #
                   outcome = 'tau', para = 'tau_s', ws = ws[2])
tau_s <- tau_s_fit$predict()


plot_compare(tau, tau_s)
df_fit_new$tau = tau
df_fit_new$tau_s = tau_s

resEE <- EE_VIM(df_fit_new)




gamma_s = df_fit_new$gamma_s
tau_s = df_fit_new$tau_s
tau = df_fit_new$tau
QA = df_fit_new$mua_hat
po = df_fit_new$po
gn = df_fit_new$pi_hat
  
gamma_s_star = df_fit_new$gamma_s_star
tau_s_star = df_fit_new$tau_s_star
tau_star = df_fit_new$tau_star
QA_star = df_fit_new$mua_hat_star
A = df_fit_new$A
Y = df_fit_new$Y

N = nrow(df_fit_new)

# ee
theta_s <- mean((po - tau_s)^2 - (po - tau)^2)
ic_ee <- 2*(tau - tau_s)*(A/gn - (1-A)/(1-gn))*(Y-QA) + (tau - tau_s)^2 - theta_s
ic_ee2 <- (po - tau_s)^2 - (po - tau)^2 - theta_s
sqrt(var(ic_ee)/N)
sqrt(var(ic_ee2)/N)

# tmle
theta_s_star <- mean(gamma_s_star - tau_s_star^2)
ic <- 2*(tau_star - tau_s_star)*(A/gn - (1-A)/(1-gn))*(Y-QA_star) + (tau_star - tau_s_star)^2 - theta_s_star
sqrt(var(ic)/N)

plot_compare(ic_ee, ic)

hist(gn)


ic_ee <- 2*(tau - tau_s)*(A/gn - (1-A)/(1-gn))*(Y-QA)
ic <- 2*(tau_star - tau_s_star)*(A/gn - (1-A)/(1-gn))*(Y-QA_star)


ic_ee <- (tau - tau_s)^2
ic <- (tau_star - tau_s_star)^2


ic_ee <- 2*(tau - tau_s)*(A/gn - (1-A)/(1-gn))
ic <- 2*(tau_star - tau_s_star)*(A/gn - (1-A)/(1-gn))

sqrt(var(ic_ee)/N)
sqrt(var(ic)/N)


ic_ee[which(ic_ee >=6000)] = 6000
ic_ee[which(ic_ee <=-6000)] = -6000


# mse
mean((Y - QA_star)^2)
mean((Y - QA)^2)

gamma_s[1:20]
gamma_s_star[1:20]


# bootstrap

res_b = rep(NA, 500)

for (i in c(1:500)){
  cat(paste0("begin ", i, "th ", "job."))
  df_sub_new = sample_n(df_sub, 1083, replace = T)
  df_fit <- fit_para(df = df_sub_new,
                     sl_Q = sl_Q, 
                     sl_g = sl_g,
                     sl_x = sl_x,
                     ws = ws[2],
                     dr = F,
                     Q_bounds = NULL,
                     g_bounds = c(0.025, 0.975),
                     tau_bounds = NULL,
                     tau_s_bounds = NULL,
                     gamma_s_bounds = NULL)
  res_b[i] <- EE_VIM(df_fit)$coef
}

saveRDS(res_b, file = "~/Repo/te_vim/data/diag_res_b.RDS")

sd(res_b)


res_b <- readRDS("~/Repo/te_vim/data/diag_res_b.RDS")


hist(res_b)


df_sub_new = sample_n(df_sub, 1083, replace = F)
df_fit <- fit_para(df = df_sub_new,
                   sl_Q = sl_Q, 
                   sl_g = sl_g,
                   sl_x = sl_x,
                   ws = ws[2],
                   dr = F,
                   Q_bounds = NULL,
                   g_bounds = c(0.025, 0.975),
                   tau_bounds = NULL,
                   tau_s_bounds = NULL,
                   gamma_s_bounds = NULL)
EE_VIM(df_fit)$coef
hist(df_fit$pi_hat)



# check Q predict
df = df_sub
folds <- origami::make_folds(strata_ids = df$A)
all_covar = setdiff(names(df), 'Y')

# setup sl3
task_Q <- sl3::make_sl3_Task(
  data = df,
  covariates = setdiff(names(df), c('Y')),
  outcome = 'Y',
  folds = folds
)

# fit Q and g
Q_fit <- sl_Q$train(task_Q)

Q0_task <- make_sl3_Task(data = df %>% mutate(A=0), covariates = all_covar)
Q1_task <- make_sl3_Task(data = df %>% mutate(A=1), covariates = all_covar)

Qbar0W <- Q_fit$predict(Q0_task)
Qbar1W <- Q_fit$predict(Q1_task)

hist(Qbar1W - Qbar0W)

hist(Qbar1W)

Q_fit$print()
