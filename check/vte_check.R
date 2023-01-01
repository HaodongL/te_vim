library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(ggplot2)
library(table1)
library(ggpubr)

rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))
source(paste0(repo_path, "R/est_function/vte.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/fit_paraloop.R"))


set.seed(1234)
N <- 500 #size of generated data
df <- generate_data_simple(N)

df_fit <- fitMods(df)
print(df_fit)


foldIDs <- getFoldIDs(N,Nfolds=5) #5 fold cross validation in this example
df_xfit <- crossFit(df,foldIDs=foldIDs)
print(df_xfit)

aipw_nocv <- VTE(df_fit,method="AIPW")
tmle_nocv <- VTE(df_fit,method="TMLE")
aipw_cv   <- VTE(df_xfit,method="AIPW")
tmle_cv   <- VTE(df_xfit,method="TMLE")



# tmle3
node_list <- list(
  W = c("X1","X2"),
  A = "A",
  Y = "Y"
)


ate_spec <- tmle_ATE(
  treatment_level = 1,
  control_level = 0
)

sl_Q <- Lrnr_sl$new(
  learners = lrnr_stack_Q,
  metalearner = ls_metalearner
)

sl_g <- Lrnr_sl$new(
  learners = lrnr_stack_g,
  metalearner = lb_metalearner,
  outcome_type = 'binomial'
)

learner_list <- list(A = sl_g, Y = sl_Q)

tmle_fit <- tmle3(ate_spec, df, node_list, learner_list)

print(tmle_fit)
test_stat <- tmle_fit$summary$tmle_est/tmle_fit$summary$se
wald <- (test_stat)^2


2*pnorm(-abs(test_stat))












