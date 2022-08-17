require(tidyverse)
repo_path = "/Users/haodongli/Repo/TE-Heterogeneity/"
source(paste0(repo_path, "R/aipw.R"))
source(paste0(repo_path, "R/tmle.R"))
source(paste0(repo_path, "R/vte.R"))
source(paste0(repo_path, "R/example_helpers.R")) #Used for the current examples

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
