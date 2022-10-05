require(tidyverse)
repo_path = "/Users/haodongli/Repo/te_vim/"
source(paste0(repo_path, "R/fit_sl3.R"))
source(paste0(repo_path, "R/example_helpers.R")) #Used for the current examples
source(paste0(repo_path, "R/fit_cate.R"))
source(paste0(repo_path, "R/vim.R"))

set.seed(1234)
N <- 1e3 #size of generated data
df <- generate_data_simple(N, print_truth = TRUE)

# y_l <- min(df$Y)
# y_u <- max(df$Y)
# df$Y <- scale01(df$Y, y_l, y_u)
covar = c('X2')
max.it = 1e3

# df_fit_sl <- fitSL(df)
# df_fit_sl <- fit_tau(df, df_fit_sl, option = "T-Learner")
# df_fit_sl <- fit_tau_s(df, df_fit_sl, covar = covar)
# df_fit_sl <- fit_gamma_s(df, df_fit_sl, covar = covar)

# res_ee <- run_EE_VIM(df, covar)
# res_tmle <- run_TMLE_VIM(df, covar, max.it)

# wrapper
res <- run_ALL_VIM(df = df, 
                   covar = covar, 
                   cv = TRUE,
                   max.it = max.it, 
                   Q_bounds = c(0.001, 0.999), 
                   g_bounds = c(0.025, 0.975),
                   tau_bounds = c(-1+1e-3, 1-1e-3),
                   tau_s_bounds = c(-1+1e-3, 1-1e-3),
                   gamma_s_bounds = c(1e-6, 1-1e-6),
                   cate_option = "DR-Learner")
res_ee <- res$resEE
res_tmle <- res$resTMLE


# test
lrnr_gam <- Lrnr_gam$new()
lrnr_gam_Q <- Lrnr_gam$new('Y ~ s(X1) + s(X2) + ti(X1,X2) + s(X1,by=A) + s(X2,by=A) + ti(X1,X2,by=A)')
lrnr_gam_A <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')


folds <- origami::make_folds(strata_ids = df$A)

task_Q <- sl3::make_sl3_Task(
  data = df,
  covariates = setdiff(names(df), c('Y')),
  outcome = 'Y',
  folds = folds
)

# gam_fit <- lrnr_gam_Q$train(task_Q)
# gam_fit$predict()

lrnr_lm_inter<- Lrnr_glm$new(formula = "~.^2")
lm_inter_fit <- lrnr_lm_inter$train(task_Q)


task_g <- sl3::make_sl3_Task(
  data = df,
  covariates = setdiff(names(df), c('Y', 'A')),
  outcome = 'A',
  folds = folds
)

sl_g <- Lrnr_sl$new(
  learners = lrnr_stack_g,
  metalearner = lb_metalearner,
  outcome_type = 'binomial'
)

gam_fit <- lrnr_gam_A$train(task_g)
gam_fit$predict()







