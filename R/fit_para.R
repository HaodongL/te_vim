# non-cv fit Q, g, po, tau, tau_s, gamma_s
fit_para <- function(df, 
                     sl_Q, 
                     sl_g,
                     sl_x,
                     ws = NULL, 
                     Q_bounds = NULL, 
                     g_bounds = c(0.025, 0.975),
                     tau_bounds = NULL, 
                     tau_s_bounds = NULL, 
                     gamma_s_bounds = NULL, 
                     dr = TRUE){
  
  all_covar = setdiff(names(df), 'Y')
  all_w = setdiff(all_covar, 'A')
  all_wsc = setdiff(all_w, ws)
  y <- df[["Y"]]
  a <- df[["A"]]
  
  res_Qg <- fitSL(df, sl_Q, sl_g, Q_bounds, g_bounds) # big object rm later
  Q_fit <- res_Qg$Q_fit # big object rm later
  g_fit <- res_Qg$g_fit # big object rm later
  
  Q0_task <- make_sl3_Task(data = df %>% mutate(A=0), covariates = all_covar)
  Q1_task <- make_sl3_Task(data = df %>% mutate(A=1), covariates = all_covar)
  
  Qbar0W <- Q_fit$predict(Q0_task)
  Qbar1W <- Q_fit$predict(Q1_task)
  # bound Q if need
  if (!is.null(Q_bounds)){
    Qbar0W <- bound(Qbar0W, Q_bounds)
    Qbar1W <- bound(Qbar1W, Q_bounds)
  }
  QbarAW <- ifelse(df$A == 1, Qbar1W, Qbar0W)
  gn <- g_fit$predict()
  # bound g
  gn <- bound(gn, g_bounds)
  
  po = (y - QbarAW)*(2*a - 1)/gn + Qbar1W - Qbar0W
  
  if (dr){
    # fit and predict tau, tau_s, gamma_s 
    tau_fit <- fit_x(df = df, sl_x = sl_x, po = po, outcome = 'po', para = 'tau')
    tau <- tau_fit$predict()
    
    tau_s_fit <- fit_x(df = df, sl_x = sl_x, po = po, outcome = 'po', para = 'tau_s', ws = ws)
    tau_s <- tau_s_fit$predict()

    gamma_s_fit <- fit_x(df = df, sl_x = sl_x, po = po, outcome = 'po', para = 'gamma_s', ws = ws)
    gamma_s <- gamma_s_fit$predict()
    
  }else{
    tau <- Qbar1W - Qbar0W
    
    tau_s_fit <- fit_x(df = df, sl_x = sl_x, tau = tau, outcome = 'tau', para = 'tau_s', ws = ws)
    tau_s <- tau_s_fit$predict()
    
    gamma_s_fit <- fit_x(df = df, sl_x = sl_x, tau = tau, outcome = 'tau', para = 'gamma_s', ws = ws)
    gamma_s <- gamma_s_fit$predict()
  }
  
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
  
  df_fit <- data.frame('Y' = df$Y, 
                       'A' = df$A, 
                       'pi_hat' = gn, 
                       'mu0_hat' = Qbar0W,
                       'mu1_hat' = Qbar1W,
                       'mua_hat' = QbarAW,
                       'po' = po,
                       'tau' = tau,
                       'tau_s' = tau_s,
                       'gamma_s' = gamma_s)
  return(df_fit)
}


# cv fit Q, g, po, tau, tau_s, gamma_s
fit_cvpara <- function(df, 
                       sl_Q, 
                       sl_g,
                       sl_x,
                       ws = NULL, 
                       Q_bounds = NULL, 
                       g_bounds = c(0.025, 0.975),
                       tau_bounds = NULL, 
                       tau_s_bounds = NULL, 
                       gamma_s_bounds = NULL, 
                       dr = TRUE){
  N = nrow(df)
  V = 10
  all_covar = setdiff(names(df), 'Y')
  all_w = setdiff(all_covar, 'A')
  all_wsc = setdiff(all_w, ws)
  
  # make folds
  folds = origami::make_folds(n = N, V = V)
  # folds <- origami::make_folds(n = N, V = V, strata_ids = df$A)
  
  # containers for cv preds
  QbarAW = rep(NA, N)
  Qbar1W = rep(NA, N)
  Qbar0W = rep(NA, N)
  gn = rep(NA, N)
  po = rep(NA, N)
  tau = rep(NA, N)
  tau_s = rep(NA, N)
  gamma_s = rep(NA, N)
  
  # fit sl models on training set, predict on validation set
  for (k in 1:V){
    ##==========================##
    ##       Q, g, po           ##
    ##==========================##
    # ---------- T ------------- #
    index_t <- folds[[k]]$training_set
    df_t <- df[index_t, ]
    df_t0 <- df_t %>% mutate(A = 0)
    df_t1 <- df_t %>% mutate(A = 1)
    y_t <- df_t[["Y"]]
    a_t <- df_t[["A"]]
    
    res_Qg <- fitSL(df_t, sl_Q, sl_g, Q_bounds, g_bounds) # big object rm later
    Q_fit <- res_Qg$Q_fit # big object rm later
    g_fit <- res_Qg$g_fit # big object rm later
    
    Q0_task_t <- make_sl3_Task(covariates = all_covar, data = df_t0)
    Q1_task_t <- make_sl3_Task(covariates = all_covar, data = df_t1)
    
    Qbar0W_t <- Q_fit$predict(Q0_task_t)
    Qbar1W_t <- Q_fit$predict(Q1_task_t)
    # bound Q if need
    if (!is.null(Q_bounds)){
      Qbar0W_t <- bound(Qbar0W_t, Q_bounds)
      Qbar1W_t <- bound(Qbar1W_t, Q_bounds)
    }
    QbarAW_t <- ifelse(df_t$A == 1, Qbar1W_t, Qbar0W_t)
    # g_task <- make_sl3_Task(covariates = all_w, data = df_t)
    gn_t <- g_fit$predict()
    # bound g
    gn_t <- bound(gn_t, g_bounds)
    po_t <- (y_t - QbarAW_t)*(2*a_t - 1)/gn_t + Qbar1W_t - Qbar0W_t
    

    # ---------- V ------------- #
    index_v <- folds[[k]]$validation_set
    df_v <- df[index_v, ]
    df_v0 <- df_v %>% mutate(A = 0)
    df_v1 <- df_v %>% mutate(A = 1)
    y_v <- df_v[["Y"]]
    a_v <- df_v[["A"]]
    
    # Q_task <- make_sl3_Task(covariates = all_covar, data = df_v)
    Q0_task_v <- make_sl3_Task(covariates = all_covar, data = df_v0)
    Q1_task_v <- make_sl3_Task(covariates = all_covar, data = df_v1)
    
    # Qn <- Q_fit$predict(Q_task)
    Qbar0W_v <- Q_fit$predict(Q0_task_v)
    Qbar1W_v <- Q_fit$predict(Q1_task_v)
    # bound Q if need
    if (!is.null(Q_bounds)){
      Qbar0W_v <- bound(Qbar0W_v, Q_bounds)
      Qbar1W_v <- bound(Qbar1W_v, Q_bounds)
    }
    QbarAW_v <- ifelse(df_v$A == 1, Qbar1W_v, Qbar0W_v)
    
    g_task_v <- make_sl3_Task(covariates = all_w, data = df_v)
    gn_v <- g_fit$predict(g_task_v)
    # bound g
    gn_v <- bound(gn_v, g_bounds)
    po_v <- (y_v - QbarAW_v)*(2*a_v - 1)/gn_v + Qbar1W_v - Qbar0W_v
    
    # update Q, g, po
    Qbar0W[index_v] <- Qbar0W_v
    Qbar1W[index_v] <- Qbar1W_v
    QbarAW[index_v] <- QbarAW_v
    gn[index_v] <- gn_v
    po[index_v]  <- po_v
    
    ##================================##
    ##       tau, tau_s, gamma_s      ##
    ##================================##
    tau_task_v <- make_sl3_Task(covariates = all_w, data = df_v)
    s_task_v <- make_sl3_Task(covariates = all_wsc, data = df_v)
    tau_t <- Qbar1W_t - Qbar0W_t
    
    if (dr){
      # fit tau, tau_s, gamma_s on training set
      # predict tau, tau_s, gamma_s on validation set
      tau_fit <- fit_x(df = df_t, sl_x = sl_x, po = po_t, outcome = 'po', para = 'tau')
      tau_v <- tau_fit$predict(tau_task_v)
      
      tau_s_fit <- fit_x(df = df_t, sl_x = sl_x, po = po_t, outcome = 'po', para = 'tau_s', ws = ws)
      tau_s_v <- tau_s_fit$predict(s_task_v)

      gamma_s_fit <- fit_x(df = df_t, sl_x = sl_x, po = po_t, outcome = 'po', para = 'gamma_s', ws = ws)
      gamma_s_v <- gamma_s_fit$predict(s_task_v)
      
      # tau_s_fit <- fit_x(df = df_t, sl_x = sl_x, tau = tau_t, outcome = 'tau', para = 'tau_s', ws = ws)
      # tau_s_v <- tau_s_fit$predict(s_task_v)
      # 
      # gamma_s_fit <- fit_x(df = df_t, sl_x = sl_x, tau = tau_t, outcome = 'tau', para = 'gamma_s', ws = ws)
      # gamma_s_v <- gamma_s_fit$predict(s_task_v)
      
    }else{
      # fit and pred
      tau_v <- Qbar1W_v - Qbar0W_v
      
      # fit and pred
      tau_s_fit <- fit_x(df = df_t, sl_x = sl_x, tau = tau_t, outcome = 'tau', para = 'tau_s', ws = ws)
      tau_s_v <- tau_s_fit$predict(s_task_v)
      
      # fit and pred
      gamma_s_fit <- fit_x(df = df_t, sl_x = sl_x, tau = tau_t, outcome = 'tau', para = 'gamma_s', ws = ws)
      gamma_s_v <- gamma_s_fit$predict(s_task_v)
    }
    
    # bound tau in [-1,1]
    if (!is.null(tau_bounds)){
      tau_v <- bound(tau_v, tau_bounds)
    }
    # bound tau_s in [-1,1]
    if (!is.null(tau_s_bounds)){
      tau_s_v <- bound(tau_s_v, tau_s_bounds)
    }
    # bound gamma_s in [0,1]
    if (!is.null(gamma_s_bounds)){
      gamma_s_v <- bound(gamma_s_v, gamma_s_bounds)
    }
    
    # update tau, tau_s, gamma_s
    tau[index_v] <- tau_v
    tau_s[index_v] <- tau_s_v
    gamma_s[index_v] <- gamma_s_v
    
    # remove big object
    rm(res_Qg, Q_fit, g_fit)
  }
  
  df_fit <- data.frame('Y' = df$Y, 
                       'A' = df$A, 
                       'pi_hat' = gn, 
                       'mu0_hat' = Qbar0W,
                       'mu1_hat' = Qbar1W,
                       'mua_hat' = QbarAW,
                       'po' = po,
                       'tau' = tau,
                       'tau_s' = tau_s,
                       'gamma_s' = gamma_s)
  return(df_fit)
}


fitSL <- function(df, sl_Q, sl_g, Q_bounds = NULL, g_bounds = c(0.025, 0.975)){
  
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
  
  # fit Q and g
  Q_fit <- sl_Q$train(task_Q)
  g_fit <- sl_g$train(task_g)
  
  res = list('Q_fit' = Q_fit,
             'g_fit' = g_fit)
  return(res)
}



fit_x <- function(df, sl_x, po = NULL, tau = NULL, outcome = 'po', para = 'tau', ws = NULL){
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
    covariates = setdiff(names(df_train), c('Y', 'A', outcome, ws)),
    outcome = outcome,
    folds = folds
  )
  
  fit <- sl_x$train(task)
  return(fit)
}

fit_gamma <- function(df, sl_x = sl_gamma, po = NULL, tau = NULL, outcome = 'po', para = 'tau', ws = NULL){
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
    covariates = setdiff(names(df_train), c('Y', 'A', outcome, ws)),
    outcome = outcome,
    folds = folds
  )
  
  fit <- sl_x$train(task)
  return(fit)
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

