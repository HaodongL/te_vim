#Function for estimating CATE
#
#df_fit is a dataframe or tibble containing:
#Y - outcome of interest
#A - treatment indicator
#pi_hat - estimate of the propensity score E(A|X)
#mu1_hat - estimate of E(Y|A=1,X)
#mu0_hat - estimate of E(Y|A=0,X)

# option is "T-Learner" or "DR-Learner"

fit_cate <- function(df, df_fit, option = "T-Learner"){
  
  # get \varphi (the pseudo outcome, or the aipw)
  y = df_fit$Y
  a = df_fit$a
  pi_hat = df_fit$pi_hat
  mu1_hat = df_fit$mu1_hat
  mu0_hat = df_fit$mu0_hat
  mua_hat = df_fit$mua_hat

  varphi = (y - mua_hat)*(a - pi_hat)/(pi_hat*(1-pi_hat)) + mu1_hat - mu0_hat
  
  if (option == "T-Learner"){
    CATE <- mu1_hat - mu0_hat
  }else if (option == "DR-Learner"){
    df_train <- cbind(df, varphi)
    folds <- origami::make_folds(strata_ids = df$A)
    
    # setup sl3
    task_cate <- sl3::make_sl3_Task(
      data = df_train,
      covariates = setdiff(names(df_train), c('Y', 'varphi')),
      outcome = 'varphi',
      folds = folds
    )
    
    # assume we have 'lrnr_stack' and 'ls_metalearner' defined outside
    sl_cate <- Lrnr_sl$new(
      learners = lrnr_stack,
      metalearner = ls_metalearner
    )
    
    cate_fit <- sl_cate$train(task_cate)
    CATE <- cate_fit$predict()
  }
  
  res <- cbind(df_fit, varphi, CATE)
  return(res)
}

