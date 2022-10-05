#Function for estimating CATE
#
#df_fit is a dataframe or tibble containing:
#Y - outcome of interest
#A - treatment indicator
#pi_hat - estimate of the propensity score E(A|X)
#mu1_hat - estimate of E(Y|A=1,X)
#mu0_hat - estimate of E(Y|A=0,X)

# option is "T-Learner" or "DR-Learner"

fit_tau <- function(df, df_fit, tau_bounds = NULL, option = "T-Learner", cv = TRUE){
  
  # get po (the pseudo outcome)
  y = df_fit$Y
  a = df_fit$A
  pi_hat = df_fit$pi_hat
  mu1_hat = df_fit$mu1_hat
  mu0_hat = df_fit$mu0_hat
  mua_hat = df_fit$mua_hat

  po = (y - mua_hat)*(a - pi_hat)/(pi_hat*(1-pi_hat)) + mu1_hat - mu0_hat
  
  if (option == "T-Learner"){
    tau <- mu1_hat - mu0_hat
  }else if (option == "DR-Learner"){
    df_train <- cbind(df, po)
    folds <- origami::folds_vfold(nrow(df))
    
    # setup sl3
    task_cate <- sl3::make_sl3_Task(
      data = df_train,
      covariates = setdiff(names(df_train), c('Y', 'A', 'po')),
      outcome = 'po',
      folds = folds
    )
    
    # assume we have 'lrnr_stack' and 'ls_metalearner' defined outside
    sl_cate <- Lrnr_sl$new(
      learners = lrnr_stack,
      metalearner = ls_metalearner
    )
    
    cate_fit <- sl_cate$train(task_cate)
    
    if (cv){
      tau <- cate_fit$predict_fold(task_cate, "validation")
    }else{
      tau <- cate_fit$predict()
    }
  }
  
  # bound tau in [-1,1]
  if (!is.null(tau_bounds)){
    tau <- bound(tau, tau_bounds)
  }
  
  res <- cbind(df_fit, po, tau)
  return(res)
}

fit_tau_s <- function(df, df_fit, tau_s_bounds = NULL, option = "T-Learner", cv = TRUE, covar = NULL){

  # get po (the pseudo outcome)
  # y = df_fit$Y
  # a = df_fit$A
  # pi_hat = df_fit$pi_hat
  # mu1_hat = df_fit$mu1_hat
  # mu0_hat = df_fit$mu0_hat
  # mua_hat = df_fit$mua_hat

  # po = (y - mua_hat)*(a - pi_hat)/(pi_hat*(1-pi_hat)) + mu1_hat - mu0_hat

  if (option == "DR-Learner"){
    po = df_fit$po
    df_train <- cbind(df, po)
    folds <- origami::folds_vfold(nrow(df))

    # setup sl3
    task_cate <- sl3::make_sl3_Task(
      data = df_train,
      covariates = setdiff(names(df_train), c('Y', 'A', 'po', covar)),
      outcome = 'po',
      folds = folds
    )

    # assume we have 'lrnr_stack' and 'ls_metalearner' defined outside
    sl_cate <- Lrnr_sl$new(
      learners = lrnr_stack,
      metalearner = ls_metalearner
    )

    cate_fit <- sl_cate$train(task_cate)
    
    if (cv){
      tau_s <- cate_fit$predict_fold(task_cate, "validation")
    }else{
      tau_s <- cate_fit$predict()
    }
    
  }else{
    tau = df_fit$tau
    df_train <- cbind(df, tau)
    folds <- origami::folds_vfold(nrow(df))

    # setup sl3
    task_cate <- sl3::make_sl3_Task(
      data = df_train,
      covariates = setdiff(names(df_train), c('Y', 'A', 'tau', covar)),
      outcome = 'tau',
      folds = folds
    )

    # assume we have 'lrnr_stack' and 'ls_metalearner' defined outside
    sl_cate <- Lrnr_sl$new(
      learners = lrnr_stack,
      metalearner = ls_metalearner
    )

    cate_fit <- sl_cate$train(task_cate)
    if (cv){
      tau_s <- cate_fit$predict_fold(task_cate, "validation")
    }else{
      tau_s <- cate_fit$predict()
    }
  }
  
  # bound tau_s in [-1,1]
  if (!is.null(tau_s_bounds)){
    tau_s <- bound(tau_s, tau_s_bounds)
  }

  res <- cbind(df_fit, tau_s)
  return(res)
}



# fit_tau_s <- function(df, df_fit, cv = TRUE, covar = NULL){
# 
#   tau = df_fit$tau
#   df_train <- cbind(df, tau)
#   folds <- origami::folds_vfold(nrow(df))
# 
#   # setup sl3
#   task_cate <- sl3::make_sl3_Task(
#     data = df_train,
#     covariates = setdiff(names(df_train), c('Y', 'A', 'tau', covar)),
#     outcome = 'tau',
#     folds = folds
#   )
# 
#   # assume we have 'lrnr_stack' and 'ls_metalearner' defined outside
#   sl_cate <- Lrnr_sl$new(
#     learners = lrnr_stack,
#     metalearner = ls_metalearner
#   )
# 
#   cate_fit <- sl_cate$train(task_cate)
#   if (cv){
#     tau_s <- cate_fit$predict_fold(task_cate, "validation")
#   }else{
#     tau_s <- cate_fit$predict()
#   }
# 
# 
#   res <- cbind(df_fit, tau_s)
#   return(res)
# }



fit_gamma_s <- function(df, df_fit, gamma_s_bounds = NULL, cv = TRUE, covar = NULL){
  
  tau_2= (df_fit$tau)^2
  df_train <- cbind(df, tau_2)
  folds <- origami::folds_vfold(nrow(df))
  
  # setup sl3
  task_gamma <- sl3::make_sl3_Task(
    data = df_train,
    covariates = setdiff(names(df_train), c('Y', 'A', 'tau_2', covar)),
    outcome = 'tau_2',
    folds = folds
  )
  
  # assume we have 'lrnr_stack' and 'ls_metalearner' defined outside
  sl_gamma <- Lrnr_sl$new(
    learners = lrnr_stack,
    metalearner = ls_metalearner
  )
  
  gamma_fit <- sl_gamma$train(task_gamma)
  
  if (cv){
    gamma_s <- gamma_fit$predict_fold(task_gamma, "validation")
  }else{
    gamma_s <- gamma_fit$predict()
  }
  # bound gamma_s in [0,1]
  if (!is.null(gamma_s_bounds)){
    gamma_s <- bound(gamma_s, gamma_s_bounds)
  }
  
  res <- cbind(df_fit, gamma_s)
  return(res)
}


bound <- function(x, bounds) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  pmin(pmax(x, lower), upper)
}

