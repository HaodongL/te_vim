library(haldensify)
# non-cv fit Q, g, po, tau, tau_s, gamma_s with HAL
fit_para_hal <- function(df, 
                     sl_Q, 
                     sl_g,
                     sl_ws,
                     ws = NULL, 
                     Q_bounds = NULL, 
                     g_bounds = c(0.025, 0.975),
                     tau_bounds = NULL, 
                     tau_s_bounds = NULL, 
                     gamma_s_bounds = NULL, 
                     dr = FALSE){
  y <- df[["Y"]]
  a <- df[["A"]]
  df <- df %>% 
    mutate_at(vars(all_of(ws)), ~as.numeric(.)) %>% 
    rename("ws" = all_of(ws))
  ws <- 'ws'
  all_covar = setdiff(names(df), 'Y')
  all_w = setdiff(all_covar, 'A')
  all_wsc = setdiff(all_w, ws)
  
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
  po = (y - QbarAW)*(a/gn - (1-a)/(1-gn)) + Qbar1W - Qbar0W
  
  ws_fit <- fit_ws(df, sl_ws, ws = ws)
  wsn <- ws_fit$predict()
  
  Q00_task <- make_sl3_Task(data = df %>% mutate(A=0, ws=0), covariates = all_covar)
  Q10_task <- make_sl3_Task(data = df %>% mutate(A=1, ws=0), covariates = all_covar)
  Q01_task <- make_sl3_Task(data = df %>% mutate(A=0, ws=1), covariates = all_covar)
  Q11_task <- make_sl3_Task(data = df %>% mutate(A=1, ws=1), covariates = all_covar)
  
  Qbar00 <- Q_fit$predict(Q00_task)
  Qbar10 <- Q_fit$predict(Q10_task)
  Qbar01 <- Q_fit$predict(Q01_task)
  Qbar11 <- Q_fit$predict(Q11_task)
  # bound Q if need
  if (!is.null(Q_bounds)){
    Qbar00 <- bound(Qbar00, Q_bounds)
    Qbar10 <- bound(Qbar10, Q_bounds)
    Qbar01 <- bound(Qbar01, Q_bounds)
    Qbar11 <- bound(Qbar11, Q_bounds)
  }
  
  tau <- Qbar1W - Qbar0W
  tau0 <- Qbar10 - Qbar00
  tau1 <- Qbar11 - Qbar01
  
  # tau_s
  tau_s = tau1 * wsn + tau0 * (1-wsn)
  
  # gamma_s
  gamma_s = tau1^2 * wsn + tau0^2 * (1-wsn)
  
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


fit_ws <- function(df, sl_ws, ws = NULL){
  # folds <- origami::folds_vfold(nrow(df))
  # folds = origami::make_folds(n = nrow(df), V = 10)
  
  # print(setdiff(names(df_train), c('Y', 'A', outcome, ws)))
  # setup sl3
  task <- sl3::make_sl3_Task(
    data = df,
    covariates = setdiff(names(df), c('Y', 'A', ws)),
    outcome = ws
  )
  
  fit <- sl_ws$train(task)
  return(fit)
}







# non-cv fit Q, g, po, tau, tau_s, gamma_s with HAL
fit_para_hal <- function(df, 
                         sl_Q, 
                         sl_g,
                         sl_ws,
                         ws = NULL, 
                         Q_bounds = NULL, 
                         g_bounds = c(0.025, 0.975),
                         tau_bounds = NULL, 
                         tau_s_bounds = NULL, 
                         gamma_s_bounds = NULL, 
                         dr = FALSE,
                         m = 100){
  n <- nrow(df)
  y <- df[["Y"]]
  a <- df[["A"]]
  df <- df %>% 
    mutate_at(vars(all_of(ws)), ~as.numeric(.)) %>% 
    rename("ws" = all_of(ws))
  ws <- 'ws'
  all_covar = setdiff(names(df), 'Y')
  all_w = setdiff(all_covar, 'A')
  all_wsc = setdiff(all_w, ws)
  
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
  po = (y - QbarAW)*(a/gn - (1-a)/(1-gn)) + Qbar1W - Qbar0W
  
  
  # estimate Ws|Wc
  Ws = df %>% select(all_of(ws)) %>% pull()
  Wc = df %>% select(-all_of(c(ws, 'A', 'Y')))
  
  # haldensify_fit <- haldensify(
  #   A = Ws, W = Wc,
  #   n_bins = 10, grid_type = "equal_range",
  #   lambda_seq = exp(seq(-1, -10, length = 100)),
  #   max_degree = 3,
  #   reduce_basis = 1 / sqrt(length(Ws))
  # )
  
  qs = quantile(Ws, seq(0, 1, length = m)) 
  qs = unname(qs)
  
  tau_qj <- matrix(NA,m,n)
  for (j in c(1:m)){
    Q0j_task <- make_sl3_Task(data = df %>% mutate(A=0, ws=qs[j]), covariates = all_covar)
    Q1j_task <- make_sl3_Task(data = df %>% mutate(A=1, ws=qs[j]), covariates = all_covar)
    Qbar0j <- Q_fit$predict(Q0j_task)
    Qbar1j <- Q_fit$predict(Q1j_task)
    # bound Q if need
    if (!is.null(Q_bounds)){
      Qbar0j <- bound(Qbar0j, Q_bounds)
      Qbar1j <- bound(Qbar1j, Q_bounds)
    }
    tau_qj[j,] <- Qbar1j - Qbar0j
  }
  # hn <- predict(haldensify_fit, new_A = qs, new_W = Wc)
  # hn = softmax(hn)
  hn <- rep(1/m, m)
  tau_s <- colSums(tau_qj * hn)
  gamma_s <- colSums(tau_qj^2 * hn)
  tau <- Qbar1W - Qbar0W
  
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

fit_ws <- function(df, sl_ws, ws = NULL){
  # folds <- origami::folds_vfold(nrow(df))
  # folds = origami::make_folds(n = nrow(df), V = 10)
  
  # print(setdiff(names(df_train), c('Y', 'A', outcome, ws)))
  # setup sl3
  task <- sl3::make_sl3_Task(
    data = df,
    covariates = setdiff(names(df), c('Y', 'A', ws)),
    outcome = ws
  )
  
  fit <- sl_ws$train(task)
  return(fit)
}


# library(haldensify)
# N <- 5e2 #size of generated data
# df <- generate_data_simple(N, print_truth = TRUE) 
# ws = c('X2')
# Ws = df %>% select(all_of(ws)) %>% pull()
# Wc = df %>% select(-all_of(c(ws, 'A', 'Y')))
# 
# haldensify_fit <- haldensify(
#   A = Ws, W = Wc,
#   n_bins = 10, grid_type = "equal_range",
#   lambda_seq = exp(seq(-1, -10, length = 100)),
#   max_degree = 3,
#   reduce_basis = 1 / sqrt(length(A))
# )
# 
# m = 11
# qs = quantile(Ws, seq(0, 1, length = m)) 
# qs = unname(qs)
# 
# tau_qj <- rep(NA, m)
# for (j in c(1:m)){
#   Q0j_task <- make_sl3_Task(data = df %>% mutate(A=0, ws=qs[j]), covariates = all_covar)
#   Q1j_task <- make_sl3_Task(data = df %>% mutate(A=1, ws=qs[j]), covariates = all_covar)
#   Qbar0j <- Q_fit$predict(Q0j_task)
#   Qbar1j <- Q_fit$predict(Q1j_task)
#   # bound Q if need
#   if (!is.null(Q_bounds)){
#     Qbar0j <- bound(Qbar0j, Q_bounds)
#     Qbar1j <- bound(Qbar0j, Q_bounds)
#   }
#   tau_qj[j] <- Qbar1j - Qbar0j
# }
# hn <- predict(haldensify_fit, new_A = qs, new_W = Wc)
# hn = softmax(hn)
# tau_s <- sum(tau_qj * hn)
# gamma_s <- sum(tau_qj^2 * hn)




# # test
# x = rnorm(1000000, 0, 1)
# p = dnorm(x)
# p = softmax(p)
# sum(x * p)


#--------------------------------------------------------------------
# softmax function: https://rpubs.com/FJRubio/softmax
#--------------------------------------------------------------------
softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}



