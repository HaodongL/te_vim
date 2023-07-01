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

df_sub <- df %>% 
  filter(arms == 1 | arms == 3) %>% 
  select(c("arms", "cd420", "age", "wtkg", "karnof", "cd40", "cd80",
           "homo", "gender", "race", "symptom", "drugs", "hemo", "oprior")) %>% 
  rename(A = "arms",
         Y = "cd420") %>% 
  mutate(A = as.numeric(A == 3))

### ------------  Part 2. Estimation ------------ ###
# vim loop over all covariates
set.seed(1)
ws = c("age", "wtkg", "karnof", "cd40", "cd80", "homo", 
       "gender", "race", "symptom", "drugs", "hemo", "oprior")

# ws = c("gender")

# ws = c("AGE", "SEX")
n_ws <- length(ws)

tic()
registerDoParallel(6)
df_vim <- foreach(i = 1:n_ws, .combine = 'rbind') %do% {
  print(paste0("Covar name: ", ws[i]))
  res <- run_VIM(df = df_sub,
                 sl_Q = sl_Q, 
                 sl_g = sl_g,
                 sl_x = sl_x,
                 ws = ws[i], 
                 cv = F,
                 dr = T,
                 max_it = 1e4)
  res_ee <- res$resEE
  res_tmle <- res$resTMLE
  res_ss <- res$resSS
  
  df_vim_i <- data.frame("varname" = c(ws[i], ws[i], ws[i]), 
                         "importance" = c(res_ee$coef, res_tmle$coef, res_ss$coef), 
                         "ci_l" = c(res_ee$ci_l, res_tmle$ci_l, res_ss$ci_l), 
                         "ci_u" = c(res_ee$ci_u, res_tmle$ci_u, res_ss$ci_u), 
                         "method" = c('EE', 'TMLE', 'SS'))
  return(df_vim_i)
}
toc()


# One
res <- run_VIM(df = df_sub,
               sl_Q = sl_Q, 
               sl_g = sl_g,
               sl_x = sl_x,
               ws = ws[2], 
               cv = F,
               dr = F,
               max_it = 1e4)
res_ee <- res$resEE
res_tmle <- res$resTMLE


df_fit <- fit_para(df = df_sub,
                   sl_Q = sl_Q, 
                   sl_g = sl_g,
                   sl_x = sl_x,
                   ws = ws[2],
                   dr = T,
                   Q_bounds = NULL,
                   g_bounds = c(0.001, 0.999),
                   tau_bounds = NULL,
                   tau_s_bounds = NULL,
                   gamma_s_bounds = NULL)

resTMLE <- TMLE_VIM(df_fit, max_it=1e4, lr=0.0001)
resEE <- EE_VIM(df_fit)
resSS <- SS_VIM(df_fit)


saveRDS(df_vim, file = "~/Repo/te_vim/data/df_vim_t_demo_sl.RDS")


df_theta <-  readRDS("~/Repo/te_vim/data/df_vim_t_demo_sl.RDS")




df_theta <- df_vim

p_theta_ee <- plot_theta(df_theta, estimator = "EE")
p_theta_tmle <- plot_theta(df_theta, estimator = "TMLE")

p_theta_all <-
  ggarrange(p_theta_ee + ggtitle("EE VIM estimates (T-learner)"),
            p_theta_tmle + ggtitle("TMLE VIM estimates (T-learner)"))

p_theta_all


p_theta_ss <- plot_theta(df_theta, estimator = "SS")




