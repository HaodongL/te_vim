require(tidyverse)
repo_path = "/Users/haodongli/Repo/te_vim/"
source(paste0(repo_path, "R/fit_sl3.R"))
source(paste0(repo_path, "R/example_helpers.R")) #Used for the current examples


set.seed(1234)
N <- 500 #size of generated data
df <- generate_data_simple(N)

df_fit_sl <- fitSL(df)