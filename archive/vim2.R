# Parameter: \Theta_s
# Estimator: TMLE and EE
run_VIM2 <- function(df, 
                     sl_Q, 
                     sl_g,
                     sl_x,
                     ws, 
                     cv = FALSE,
                     dr = TRUE,
                     lfm_linear = TRUE,
                     tmle_b = FALSE,
                     max_it = 1e4, 
                     lr = 1e-4,
                     Q_bounds = c(1e-4, 1-1e-4), 
                     g_bounds = c(0.025, 0.975),
                     tau_bounds = c(-1+1e-4, 1-1e-4),
                     tau_s_bounds = c(-1+1e-4, 1-1e-4),
                     gamma_s_bounds = c(1e-8, 1-1e-8)){
  
  if (lfm_linear){
    Q_bounds <- NULL
    tau_bounds <- NULL
    tau_s_bounds <- NULL
    gamma_s_bounds <- NULL
    y_l <- 0
    y_u <- 1
  }else{
    # scale Y between [0,1]
    y_l <- min(df$Y)
    y_u <- max(df$Y)
    df$Y <- scale01(df$Y, y_l, y_u)
  }
  
  # fit Q, g
  # fit tau, tau_s, gamma_s
  if (cv){
    fit_func <- fit_cvpara
  }else{
    fit_func <- fit_para
  }
  df_fit <- fit_func(df = df,
                     sl_Q = sl_Q, 
                     sl_g = sl_g,
                     sl_x = sl_x,
                     ws = ws,
                     dr = dr,
                     Q_bounds = Q_bounds,
                     g_bounds = g_bounds,
                     tau_bounds = tau_bounds,
                     tau_s_bounds = tau_s_bounds,
                     gamma_s_bounds = gamma_s_bounds)
  
  # tmle_vim
  resTMLE <- TMLE_VIM2(df_fit, y_l, y_u, max_it, lr)
  
  # ee
  resEE <- EE_VIM2(df_fit)
  
  # ss
  resSS <- SS_VIM2(df_fit)
  
  res <- list('resTMLE' = resTMLE, 
              'resEE' = resEE,
              'resSS' = resSS)
  return(res)
}


# Estimator: SS
SS_VTE <- function(data){
  N <- NROW(data)
  A <- data$A
  Y <- data$Y
  
  # Q, g, tau, tau_s, gamma_s
  QA_0 <- data$mua_hat
  Q1_0 <- data$mu1_hat
  Q0_0 <- data$mu0_hat
  gn <- data$pi_hat
  
  tau_0 <- data$tau
  tau_s_0 <- data$tau_s
  gamma_s_0 <- data$gamma_s
  ate <- mean(tau_0)
  
  theta_s_0 <- mean((tau_0 - ate)^2)
  ic <- 2*(tau_0 - ate)*(A/gn - (1-A)/(1-gn))*(Y-QA_0) + (tau_0 - ate)^2 - theta_s_0
  ss <- sqrt(var(ic)/N)
  
  coef <- theta_s_0
  std_err <- ss
  names(coef) <- names(std_err) <- c("SS_Theta_s")
  
  out<- list(
    coef = coef,
    std_err = std_err,
    ci_l = coef - 1.96*std_err,
    ci_u = coef + 1.96*std_err,
    ic = ic
  )
  
  return(out)
}

# Estimator: EE
EE_VTE <- function(data){
  N <- NROW(data)
  po <- data$po
  tau <- data$tau
  tau_s <- data$tau_s
  
  ATE <- sum(po)/N
  Sig1 <- sum(po^2)/N - ATE^2
  
  VIM <- mean((po - ATE)^2 - (po - tau)^2)
  ic <- (po - ATE)^2 - (po - tau)^2 - VIM
  ss <- sqrt(var(ic)/N)
  
  coef <- VIM
  std_err <- ss
  names(coef) <- names(std_err) <- c("VIM_Theta_s")
  
  out<- list(
    coef = coef,
    std_err = std_err,
    ci_l = coef - 1.96*std_err,
    ci_u = coef + 1.96*std_err,
    ic = ic
  )
  
  return(out)
}



# Estimator: TMLE
TMLE_VTE <- function(data, y_l, y_u, max_it = 600, lr = 1e-4){
  N <- nrow(data)
  A <- data$A
  Y <- data$Y
  
  # Q, g, tau, tau_s, gamma_s
  QA_0 <- data$mua_hat
  Q1_0 <- data$mu1_hat
  Q0_0 <- data$mu0_hat
  gn <- data$pi_hat
  
  tau_0 <- data$tau
  ate_0 <- rep(mean(tau_0), N)
  
  # sig1
  sig1 <- sd(2*(tau_0 - ate_0)*(A/gn - (1-A)/(1-gn))*(Y - QA_0))
  # eps1
  eps1 <- lr
  
  # i
  i <- 1
  QA_i <- QA_0
  Q1_i <- Q1_0
  Q0_i <- Q0_0
  tau_i <- tau_0
  ate_i <- ate_0
  
  # logging <- matrix(NA, max_it, 3)
  
  while(i <= max_it){
    # step 1. update Q and tau
    if (i == 1){
      HQ_1i <- 2*(tau_i - ate_i)
      H1_1i <- 2*(tau_i - ate_i)*(1/gn)
      H0_1i <- 2*(tau_i - ate_i)*(-1/(1-gn))
      HA_1i <- ifelse(A, H1_1i, H0_1i)
      D_1i <- HA_1i*(Y - QA_i)
      PnD_1i <- mean(D_1i)
    }
    
    Q1_i <- Q1_i + eps1*(HQ_1i)*sign(PnD_1i)
    Q0_i <- Q0_i + eps1*(-HQ_1i)*sign(PnD_1i)
    QA_i <- ifelse(A, Q1_i, Q0_i)
    
    tau_i <- Q1_i - Q0_i
    ate_i <- rep(mean(tau_i), N)
    
    # update D1 D2, check criteria
    HQ_1i <- 2*(tau_i - ate_i)
    H1_1i <- 2*(tau_i - ate_i)*(1/gn)
    H0_1i <- 2*(tau_i - ate_i)*(-1/(1-gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i)
    D_1i <- HA_1i*(Y - QA_i)
    PnD_1i <- mean(D_1i)
    
    c1 <- abs(PnD_1i) <= sig1/(sqrt(N)*log(N))
    
    if (c1){
      break
    }
    i <- i + 1
  }
  
  QA_star <- QA_i
  Q1_star <- Q1_i
  Q0_star <- Q0_i
  tau_star <- tau_i
  ate_star <- rep(mean(tau_i), N)
  gamma_s_star <- rep(mean(tau_i^2), N)
  
  if(i>=max_it) warning("Max iterations reached in TMLE")
  
  # calculate Theta_s, scale back
  theta_s_star <- mean(gamma_s_star - ate_star^2)*(y_u-y_l)^2
  ic <- 2*(tau_star - ate_star)*(A/gn - (1-A)/(1-gn))*(Y-QA_star) + (tau_star - ate_star)^2 - theta_s_star
  ss <- sqrt(var(ic)/N)*(y_u-y_l)^2
  
  coef <- theta_s_star
  std_err <- ss
  names(coef) <- names(std_err) <- c("VTE_Theta_s")
  
  out<- list(
    coef = coef,
    std_err = std_err,
    ci_l = coef - 1.96*std_err,
    ci_u = coef + 1.96*std_err,
    ic = ic
  )
  
  return(out)
}

TMLE_VIM2 <- function(data, y_l, y_u, max_it = 600, lr = 1e-4){
  N <- nrow(data)
  res_vim1 <- TMLE_VIM_a(data, y_l, y_u, max_it, lr)
  res_vim2 <- TMLE_VTE(data, y_l, y_u, max_it, lr)
  
  psi1 <- res_vim1$coef
  psi2 <- res_vim2$coef
  
  ic1 <- res_vim1$ic
  ic2 <- res_vim2$ic
  
  ic <- df_log_rr(x = list(psi2, psi1),
                  dx = list(ic2, ic1))
  
  psi <- log(psi1/psi2)
  se <- sqrt(var(ic)/N)
  
  ci_l <- psi - 1.96*se
  ci_u <- psi + 1.96*se
  
  out<- list(
    coef = exp(psi),
    std_err = se,
    ci_l = exp(ci_l),
    ci_u = exp(ci_u),
    ic = ic
  )
  print(psi1)
  print(psi2)
  return(out)
}

EE_VIM2 <- function(data){
  N <- nrow(data)
  res_vim1 <- EE_VIM(data)
  res_vim2 <- EE_VTE(data)
  
  psi1 <- res_vim1$coef
  psi2 <- res_vim2$coef
  
  ic1 <- res_vim1$ic
  ic2 <- res_vim2$ic
  
  ic <- df_log_rr(x = list(psi2, psi1),
                  dx = list(ic2, ic1))
  
  psi <- log(psi1/psi2)
  se <- sqrt(var(ic)/N)
  
  ci_l <- psi - 1.96*se
  ci_u <- psi + 1.96*se
  
  out<- list(
    coef = exp(psi),
    std_err = se,
    ci_l = exp(ci_l),
    ci_u = exp(ci_u),
    ic = ic
  )
  return(out)
}

SS_VIM2 <- function(data){
  N <- nrow(data)
  res_vim1 <- SS_VIM(data)
  res_vim2 <- SS_VTE(data)
  
  psi1 <- res_vim1$coef
  psi2 <- res_vim2$coef
  
  ic1 <- res_vim1$ic
  ic2 <- res_vim2$ic
  
  ic <- df_log_rr(x = list(psi2, psi1),
                  dx = list(ic2, ic1))
  
  psi <- log(psi1/psi2)
  se <- sqrt(var(ic)/N)
  
  ci_l <- psi - 1.96*se
  ci_u <- psi + 1.96*se
  
  out<- list(
    coef = exp(psi),
    std_err = se,
    ci_l = exp(ci_l),
    ci_u = exp(ci_u),
    ic = ic
  )
  return(out)
}

# Risk Ratio EY1/EY0
f_log_rr <- function(x, dx) {
  log(x[[2]]) - log(x[[1]])
}

df_log_rr <- function(x, dx) {
  dx[[2]] / x[[2]] - dx[[1]] / x[[1]]
}


# TMLE_VIM2 <- function(data, y_l, y_u, max_it = 600){
#   
# }
