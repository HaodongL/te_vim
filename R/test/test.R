require(tidyverse)
repo_path = "/Users/haodongli/Repo/te_vim/"
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


# SL
source(paste0(repo_path, "R/fit_sl3.R"))

df_fit_sl <- fitSL(df)
print(df_fit_sl)
tmle_nocv_sl <- VTE(df_fit_sl,method="TMLE")

foldIDs <- getFoldIDs(N,Nfolds=5) #5 fold cross validation in this example
df_xfit_sl <- crossFitSL(df,foldIDs=foldIDs)
print(df_xfit_sl)

aipw_nocv_sl <- VTE(df_fit_sl,method="AIPW")
tmle_nocv_sl <- VTE(df_fit_sl,method="TMLE")
aipw_cv_sl   <- VTE(df_xfit_sl,method="AIPW")
tmle_cv_sl   <- VTE(df_xfit_sl,method="TMLE")











# par(mfrow=c(1,2))
# hist(df_fit_sl$mu1_hat, xlim = c(-1,6))
# hist(df_fit$mu1_hat, xlim = c(-1,6))
# 
# 
# par(mfrow=c(1,2))
# hist(df_fit_sl$mu0_hat, xlim = c(-2,6))
# hist(df_fit$mu0_hat, xlim = c(-2,6))
# 
# 
# mean(df_fit_sl$mu1_hat - df_fit_sl$mu0_hat)
# mean(df_fit$mu1_hat - df_fit$mu0_hat)
# 
# 
# 
# mod.m <- gam(Y~s(X1) + s(X2) + ti(X1,X2)+s(X1,by=A) + s(X2,by=A) + ti(X1,X2,by=A),
#              family = gaussian(),data=df)
# 
# 
# mean((predict(mod.m,df,type="response") - df$Y)^2)
# 
# mod.m <- ranger::ranger(Y ~ .,data = df, num.trees =  200)
# mu1_hat = predict(mod.m,mutate(df,A=1),type="response")
# 
# 
# # setup sl3
# task_Q <- sl3::make_sl3_Task(
#   data = df,
#   covariates = setdiff(names(df), c('Y')),
#   outcome = 'Y',
#   folds = origami::make_folds(strata_ids = df$A)
# )
# 
# sl_Q <- Lrnr_sl$new(
#   learners = lrnr_stack,
#   metalearner = ls_metalearner,
#   outcome_type = 'continuous'
# )
# 
# # fit Q and g
# Q_fit <- sl_Q$train(task_Q)
# 
# mean((Q_fit$predict() - df$Y)^2)




