#
temp <- fitSL(df)
df_fit <- temp$df_fit
Q_fit <- temp$Q_fit
g_fit <- temp$g_fit
Q_task <- temp$Q_task
g_task <- temp$g_task
tau_bounds = NULL
tau_s_bounds = NULL
gamma_s_bounds = NULL
dr = TRUE
cv = TRUE

fit_tau_gamma <- function(df, 
                          df_fit, 
                          covar, 
                          Q_fit, 
                          g_fit, 
                          Q_task, 
                          g_task, 
                          tau_bounds = NULL, 
                          tau_s_bounds = NULL, 
                          gamma_s_bounds = NULL, 
                          dr = TRUE, 
                          cv = TRUE){
  # Qn <- Q_fit$predict()
  # gn <- g_fit$predict()
  y_full = df$Y
  a_full = df$A
  all_covar = Q_task$nodes$covariates
  all_w = g_task$nodes$covariates
  
  cv_preds_list <- lapply(seq_along(Q_task$folds), function(k){
    # get training dataset for fold k:
    t_index_k <- Q_task$folds[[k]]$training_set
    t_data <- Q_task$data[t_index_k, ]
    t_data0 <- t_data %>%  mutate(A = 0)
    t_data1 <- t_data %>%  mutate(A = 1)
    y <- t_data[["Y"]]
    a <- t_data[["A"]]
    
    t_task <- make_sl3_Task(covariates = all_covar, data = t_data)
    t_task0 <- make_sl3_Task(covariates = all_covar, data = t_data0)
    t_task1 <- make_sl3_Task(covariates = all_covar, data = t_data1)
    t_task_g <- make_sl3_Task(covariates = all_w, data = t_data)
    
    Qn <- Q_fit$fit_object$cv_fit$predict_fold(task = t_task, fold_number = k)
    Q0n <- Q_fit$fit_object$cv_fit$predict_fold(task = t_task, fold_number = k)
    Q1n <- Q_fit$fit_object$cv_fit$predict_fold(task = t_task, fold_number = k)
    gn <- g_fit$fit_object$cv_fit$predict_fold(task = t_task_g, fold_number = k)
    
    
    Qn <- Q_fit$predict_fold(task = t_task, fold_number = k)
    Q0n <- Q_fit$predict_fold(task = t_task0, fold_number = k)
    Q1n <- Q_fit$predict_fold(task = t_task1, fold_number = k)
    
    po = (y - Qn)*(a - gn)/(gn*(1-gn)) + Q1n - Q0n
    
    
    # get validation dataset for fold i:
    v_data <- Q_task$data[Q_task$folds[[k]]$validation_set, ]
    v_task <- make_sl3_Task(covariates = setdiff(all_w, covar), 
                            data = v_data)
    
    if (dr){
      # fit tau, tau_s, gamma_s on training set
      tau_fit <- fit_x(df = t_data, po = po, outcome = 'po', para = 'tau')
      tau_s_fit <- fit_x(df = t_data, po = po, outcome = 'po', para = 'tau_s', covar = covar)
      gamma_s_fit <- fit_x(df = t_data, po = po, outcome = 'po', para = 'gamma_s', covar = covar)
      
      # predict tau, tau_s, gamma_s on validation set
      tau <- tau_fit$fit_object$cv_fit$predict_fold(task = v_task, fold_number = k)
      tau_s <- tau_s_fit$fit_object$cv_fit$predict_fold(task = v_task, fold_number = k)
      gamma_s <- gamma_s_fit$fit_object$cv_fit$predict_fold(task = v_task, fold_number = k)
      
    }else{
      # fit and pred
      Q1n <- Q_fit$fit_object$cv_fit$predict_fold(task = v_task, fold_number = k)
      Q0n <- Q_fit$fit_object$cv_fit$predict_fold(task = v_task, fold_number = k)
      tau <- Q1n - Q0n
      
      # fit and pred
      tau_s_fit <- fit_x(df = t_data, tau = tau, outcome = 'tau', para = 'tau_s', covar = covar)
      tau_s <- tau_s_fit$fit_object$cv_fit$predict_fold(task = v_task, fold_number = k)
      
      # fit and pred
      gamma_s_fit <- fit_x(df = t_data, tau = tau, outcome = 'tau', para = 'gamma_s', covar = covar)
      gamma_s <- gamma_s_fit$fit_object$cv_fit$predict_fold(task = v_task, fold_number = k)
    }
    
    return(list("tau" = tau, 
                "tau_s" = tau_s, 
                "gamma_s" = gamma_s, 
                "v_index" = sl_task$folds[[k]]$validation_set))
  })
  # extract the validation set predictions across all folds
  tau <- do.call(rbind, lapply(cv_preds_list, "[[", "tau"))
  tau_s <- do.call(rbind, lapply(cv_preds_list, "[[", "tau_s"))
  gamma_s <- do.call(rbind, lapply(cv_preds_list, "[[", "gamma_s"))
  
  # extract the indices of validation set observations across all folds
  # then reorder cv_preds_byhand to correspond to the ordering in the data
  row_index_in_data <- unlist(lapply(cv_preds_list, "[[", "v_index"))
  tau <- tau[order(row_index_in_data), ]
  tau_s <- tau_s[order(row_index_in_data), ]
  gamma_s <- gamma_s[order(row_index_in_data), ]
  
  # bound tau in [-1,1]
  if (!is.null(tau_bounds)){
    tau <- bound(tau, tau_bounds)
  }
  # bound tau_s in [-1,1]
  if (!is.null(tau_s_bounds)){
    tau_s <- bound(tau_s, tau_s_bounds)
  }
  # bound gamma_s in [0,1]
  if (!is.null(gamma_s_bounds)){
    gamma_s <- bound(gamma_s, gamma_s_bounds)
  }
  
  df_fit <- cbind(df_fit, tau, tau_s, gamma_s)
  return(df_fit)
}


# 

fit_x <- function(df, po = NULL, tau = NULL, outcome = 'po', para = 'tau', covar = NULL){
  assertthat::assert_that(outcome %in% c('po', 'tau'))
  assertthat::assert_that(para %in% c('tau', 'tau_s', 'gamma_s'))
  if (outcome == 'po'){
    assertthat::assert_that(!is.null(po))
    if (para == 'gamma_s'){
      po <- po^2
    }
    df_train <- cbind(df, po)
  }else{
    if (para == 'gamma_s'){
      tau <- tau^2
    }
    df_train <- cbind(df, tau)
  }
  
  folds <- origami::folds_vfold(nrow(df))
  
  # setup sl3
  task <- sl3::make_sl3_Task(
    data = df_train,
    covariates = setdiff(names(df_train), c('Y', 'A', outcome, covar)),
    outcome = outcome,
    folds = folds
  )
  
  # assume we have 'lrnr_stack' and 'ls_metalearner' defined outside
  sl <- Lrnr_sl$new(
    learners = lrnr_stack,
    metalearner = ls_metalearner
  )
  fit <- sl$train(task_cate)
  return(fit)
}


# -- sl modeling function
fitSL <- function(df, Q_bounds = NULL, g_bounds = c(0.025, 0.975), cv = TRUE){
  
  folds <- origami::make_folds(strata_ids = df$A)
  # folds <- origami::folds_vfold(nrow(df))
  
  # setup sl3
  task_Q <- sl3::make_sl3_Task(
    data = df,
    covariates = setdiff(names(df), c('Y')),
    outcome = 'Y',
    folds = folds
  )
  
  task_g <- sl3::make_sl3_Task(
    data = df,
    covariates = setdiff(names(df), c('Y', 'A')),
    outcome = 'A',
    folds = folds
  )
  
  sl_Q <- Lrnr_sl$new(
    learners = lrnr_stack_Q,
    metalearner = ls_metalearner
  )
  
  sl_g <- Lrnr_sl$new(
    learners = lrnr_stack_g,
    metalearner = lb_metalearner,
    outcome_type = 'binomial'
  )
  
  # fit Q and g
  Q_fit <- sl_Q$train(task_Q)
  g_fit <- sl_g$train(task_g)
  
  pred_Q_cf <- pred_Q_cf(df = df, Q_fit = Q_fit, folds = folds, cv = cv)
  
  # QbarAW <- pred_Q
  Qbar1W <- pred_Q_cf$Qbar1W
  Qbar0W <- pred_Q_cf$Qbar0W
  QbarAW <- ifelse(df$A == 1, Qbar1W, Qbar0W)
  
  # bound Q (is there a better way???)
  if (!is.null(Q_bounds)){
    QbarAW <- bound(QbarAW, Q_bounds)
    Qbar1W <- bound(Qbar1W, Q_bounds)
    Qbar0W <- bound(Qbar0W, Q_bounds)
  }
  
  if (cv){
    pred_g <- g_fit$predict_fold(task_g, "validation")
  }else{
    pred_g <- g_fit$predict()
  }
  
  # bound g
  pred_g <- bound(pred_g, g_bounds)
  
  df_fit <- tibble(
    Y = df$Y, 
    A = df$A,
    pi_hat  = pred_g,
    mu1_hat = Qbar1W,
    mu0_hat = Qbar0W,
    mua_hat = QbarAW)
  
  res <- list('df_fit' = df_fit,
              'Q_task' = task_Q,
              'g_task' = task_g,
              'Q_fit' = Q_fit,
              'g_fit' = g_fit)
  return(res)
}







