library(speff2trial)
require(parallel)
require(tidyverse)
require(SuperLearner)


rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "oh/TevimFitting/fitCATE.R"))
source(paste0(repo_path, "oh/TevimFitting/fitmods_sl.R"))
source(paste0(repo_path, "oh/TevimFitting/tevim_fit_noCV.R"))

# repo_path = "~/Repo/te_vim/"
# source(paste0(repo_path, "R/est_function/sl3_config.R"))
# source(paste0(repo_path, "R/est_function/fit_para.R"))
# source(paste0(repo_path, "R/est_function/fit_hal.R"))
# source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/vim.R"))




data(ACTG175,package="speff2trial") #install.packages("speff2trial")

##Subset data and define Treatment and outcome variables
df <- as_tibble(ACTG175) %>% filter(arms %in% c(1,3)) %>%
  mutate(A = as.numeric(arms==1),
         Y = sqrt(cd420))

#Create covariate matrices
x <- df %>% select(age,wtkg,karnof,cd40,cd80,
                   gender,homo,race,symptom,drugs,str2,hemo)

#For the sub CATEs
cov_list <- list("age","wtkg","karnof","cd40","cd80",
                 "gender","homo","race","symptom","drugs","str2","hemo")
names(cov_list) <- cov_list


RANGERlearners = create.Learner("SL.ranger",params=list(num.trees=2000),
                                tune=list(mtry=c(3,4)),
                                name_prefix = "RANGER")

GAMlearners = create.Learner("SL.gam",tune = list(deg.gam=c(2,3,4)),
                             name_prefix = "GAM")

XGBlearners = create.Learner("SL.xgboost",params=list(minobspernode=10,ntrees = 2000,shrinkage=0.01),
                             tune = list(max_depth =c(2,3)),
                             name_prefix = "XGB")

SL.glmnet.interaction <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, 
                                   nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + .^2, X)
    newX <- model.matrix(~-1 + .^2, newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet.interaction"
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.glmnet.interaction <- function (object, newdata, remove_extra_cols = T, add_missing_cols = T,...){
  SuperLearner:::.SL.require("glmnet")
  if (!is.matrix(newdata)) {
    newdata <- model.matrix(~-1 + .^2, newdata)
  }
  original_cols = rownames(object$object$glmnet.fit$beta)
  if (remove_extra_cols) {
    extra_cols = setdiff(colnames(newdata), original_cols)
    if (length(extra_cols) > 0) {
      warning(paste("Removing extra columns in prediction data:", 
                    paste(extra_cols, collapse = ", ")))
      newdata = newdata[, !colnames(newdata) %in% extra_cols, 
                        drop = FALSE]
    }
  }
  if (add_missing_cols) {
    missing_cols = setdiff(original_cols, colnames(newdata))
    if (length(missing_cols) > 0) {
      warning(paste("Adding missing columns in prediction data:", 
                    paste(missing_cols, collapse = ", ")))
      new_cols = matrix(0, nrow = nrow(newdata), ncol = length(missing_cols))
      colnames(new_cols) = missing_cols
      newdata = cbind(newdata, new_cols)
      newdata = newdata[, original_cols, drop = FALSE]
    }
  }
  pred <- predict(object$object, newx = newdata, type = "response", 
                  s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
  return(pred)
}

environment(SL.glmnet.interaction) <- asNamespace("SuperLearner")
environment(predict.SL.glmnet.interaction) <- asNamespace("SuperLearner")

SL.library <- c("SL.glmnet", "SL.glm","SL.glmnet.interaction",
                RANGERlearners$names,
                GAMlearners$names,
                XGBlearners$names)

Nf  <- 20
SLf <-10

set.seed(1234567)
cl <- makeForkCluster(3)
test_sq <- TEVIM_noCV_SL(df$Y^2,df$A,x,SL.library,cov_list,
                        SLfolds=SLf,cluster=cl,SL_type = "both")

saveRDS(test_sq, file = "~/Repo/te_vim/data/test_sq.RDS")

test_sq <-  readRDS("~/Repo/te_vim/data/test_sq.RDS")

to_df_fit <- function(test_sq){
  mua_hat <- ifelse(test_sq$Initial_fit$A == 1,
                    test_sq$Initial_fit$mu1_hat,
                    test_sq$Initial_fit$mu0_hat)
    
  df_fit <- data.frame("Y" = test_sq$Initial_fit$Y,
                       "A" = test_sq$Initial_fit$A,
                       "pi_hat" = test_sq$Initial_fit$pi_hat,
                       "mu1_hat" = test_sq$Initial_fit$mu1_hat,
                       "mu0_hat" = test_sq$Initial_fit$mu0_hat,
                       "mua_hat" = mua_hat,
                       "po" = test_sq$DR_non_discrete$cate_fit$PO,
                       "tau" = test_sq$DR_non_discrete$cate_fit$CATE,
                       "tau_s" = test_sq$DR_non_discrete$maginal_fit$wtkg,
                       "gamma_s" = test_sq$DR_non_discrete$gamma_fit$wtkg)
  return(df_fit)
}

df_fit <- to_df_fit(test_sq)
resEE <- EE_VIM(df_fit)
resTMLE <- TMLE_VIM(df_fit, max_it=1e4, lr=0.0001)



# foldIDs <- getFoldIDs(NROW(df),Nfolds=Nf)
# cl <- makeForkCluster(12)
# test_sq <- TEVIMcrossFit_SL(df$Y^2,df$A,x,SL.library,
#                             foldIDs = foldIDs,cov_list,
#                             SLfolds=SLf,cluster=cl,SL_type = "both",
#                             PS.method="OutofSampleMean")


