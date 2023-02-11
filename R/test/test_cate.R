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
res <- run_VIM(df = df,
               sl_Q = sl_Q, 
               sl_g = sl_g,
               sl_x = sl_x,
               ws = c('X2'), 
               cv = F,
               dr = F,
               max_it = 1e4)
toc()
res_ee <- res$resEE
res_tmle <- res$resTMLE
res_ss <- res$resSS





# test
# lrnr_gam <- Lrnr_gam$new()
# lrnr_gam_Q <- Lrnr_gam$new('Y ~ s(X1) + s(X2) + ti(X1,X2) + s(X1,by=A) + s(X2,by=A) + ti(X1,X2,by=A)')
# lrnr_gam_A <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# 
# 
# folds <- origami::make_folds(strata_ids = df$A)
# 
# task_Q <- sl3::make_sl3_Task(
#   data = df,
#   covariates = setdiff(names(df), c('Y')),
#   outcome = 'Y',
#   folds = folds
# )
# 
# # gam_fit <- lrnr_gam_Q$train(task_Q)
# # gam_fit$predict()
# 
# lrnr_lm_inter<- Lrnr_glm$new(formula = "~.^2")
# lm_inter_fit <- lrnr_lm_inter$train(task_Q)
# 
# 
# task_g <- sl3::make_sl3_Task(
#   data = df,
#   covariates = setdiff(names(df), c('Y', 'A')),
#   outcome = 'A',
#   folds = folds
# )
# 
# sl_g <- Lrnr_sl$new(
#   learners = lrnr_stack_g,
#   metalearner = lb_metalearner,
#   outcome_type = 'binomial'
# )
# 
# gam_fit <- lrnr_gam_A$train(task_g)
# gam_fit$predict()







