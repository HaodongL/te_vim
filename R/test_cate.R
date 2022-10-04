require(tidyverse)
repo_path = "/Users/haodongli/Repo/te_vim/"
source(paste0(repo_path, "R/fit_sl3.R"))
source(paste0(repo_path, "R/example_helpers.R")) #Used for the current examples
source(paste0(repo_path, "R/fit_cate.R"))
source(paste0(repo_path, "R/vim.R"))

set.seed(1234)
N <- 1e6 #size of generated data
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

res_ee <- run_EE_VIM(df, covar)
res_tmle <- run_TMLE_VIM(df, covar, max.it)



