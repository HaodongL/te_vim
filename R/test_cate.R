require(tidyverse)
repo_path = "/Users/haodongli/Repo/te_vim/"
source(paste0(repo_path, "R/fit_sl3.R"))
source(paste0(repo_path, "R/example_helpers.R")) #Used for the current examples
source(paste0(repo_path, "R/fit_cate.R"))

set.seed(1234)
N <- 500 #size of generated data
df <- generate_data_simple(N)

df_fit_sl <- fitSL(df)
df_fit_sl <- fit_tau(df, df_fit_sl, option = "T-Learner")
df_fit_sl <- fit_tau_s(df, df_fit_sl, covar = c('X2'))

