# Parameter: \Theta_s
# Estimator: TMLE and EE
run_VIM <- function(df, 
                    sl_Q, 
                    sl_g,
                    sl_x,
                    ws, 
                    cv = FALSE,
                    dr = TRUE,
                    lfm_linear = TRUE, # always true, deprecate the logistic fluctuation now
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
  resTMLE <- TMLE_VIM(df_fit, max_it, lr)
  
  # # TEMP
  # resSS <- TMLE_VIM(df_fit, max_it, lr, logtrans = T)
  
  if (y_l != 0 | y_u != 1){
    df_fit$Y <- rescale(df_fit$Y, y_l, y_u)
    df_fit$tau <- df_fit$tau*(y_u - y_l)
    df_fit$tau_s <- df_fit$tau_s*(y_u - y_l)
    df_fit$po <- df_fit$po*(y_u - y_l)
    df_fit$gamma_s <- df_fit$gamma_s*(y_u - y_l)^2
  }
  
  # ee_vim
  resEE <- EE_VIM(df_fit)
  # ss_vim
  resSS <- SS_VIM(df_fit)
  
  
  res <- list('resTMLE' = resTMLE, 
              'resEE' = resEE,
              'resSS' = resSS)
  return(res)
}


# Estimator: SS
SS_VIM <- function(data){
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
  
  # theta_s_hat <- mean(gamma_s_0 - tau_s_0^2)
  # theta_s_hat <- mean((tau_0 - tau_s_0)^2)
  # theta_s_0 <- mean(gamma_s_0 - tau_s_0^2)
  theta_s_0 <- mean((tau_0 - tau_s_0)^2)
  ic <- 2*(tau_0 - tau_s_0)*(A/gn - (1-A)/(1-gn))*(Y-QA_0) + (tau_0 - tau_s_0)^2 - theta_s_0
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
EE_VIM <- function(data){
  N <- NROW(data)
  po <- data$po
  tau <- data$tau
  tau_s <- data$tau_s
  
  ATE <- sum(po)/N
  Sig1 <- sum(po^2)/N - ATE^2
  
  # VIM <- mean((po - tau_s)^2 - (po - tau)^2)
  # ic <- (po - tau_s)^2 - (po - tau)^2 - VIM
  A <- data$A
  Y <- data$Y
  QA <- data$mua_hat
  gn <- data$pi_hat
  
  
  VIM <- mean(2*(tau - tau_s)*(A/gn - (1-A)/(1-gn))*(Y-QA) + (tau - tau_s)^2)
  ic <- 2*(tau - tau_s)*(A/gn - (1-A)/(1-gn))*(Y-QA) + (tau - tau_s)^2 - VIM
  ss <- sqrt(var(ic)/N)
  
  coef <- VIM
  std_err <- ss
  names(coef) <- names(std_err) <- c("VIM_Theta_s")
  
  out <- list(
    coef = coef,
    std_err = std_err,
    ci_l = coef - 1.96*std_err,
    ci_u = coef + 1.96*std_err,
    ic = ic
  )
  return(out)
}

# Estimator: TMLE
TMLE_VIM <- function(data, max_it = 600, lr = 1e-4, logtrans = FALSE){
  N <- nrow(data)
  A <- data$A
  Y <- data$Y
  po <- data$po
  
  # Q, g, tau, tau_s, gamma_s
  QA_0 <- data$mua_hat
  Q1_0 <- data$mu1_hat
  Q0_0 <- data$mu0_hat
  gn <- data$pi_hat
  
  tau_0 <- data$tau
  tau_s_0 <- data$tau_s
  gamma_s_0 <- data$gamma_s
  
  # sig1, sig2
  sig1 <- sd(2*(tau_0 - tau_s_0)*(A/gn - (1-A)/(1-gn))*(Y - QA_0))
  sig2 <- sd(2*tau_s_0*(tau_0 - tau_s_0))
  
  # eps1,eps2
  eps1 <- lr
  eps2 <- lr
  
  # i
  i <- 1
  QA_i <- QA_0
  Q1_i <- Q1_0
  Q0_i <- Q0_0
  tau_i <- tau_0
  tau_s_i <- tau_s_0
  
  # logging <- matrix(NA, max_it, 3)
  
  while(i <= max_it){
    # step 1. update Q and tau
    if (i == 1){
      HQ_1i <- 2*(tau_i - tau_s_i)
      H1_1i <- 2*(tau_i - tau_s_i)*(1/gn)
      H0_1i <- 2*(tau_i - tau_s_i)*(-1/(1-gn))
      HA_1i <- ifelse(A, H1_1i, H0_1i)
      D_1i <- HA_1i*(Y - QA_i)
      PnD_1i <- mean(D_1i)
    }
    
    Q1_i <- Q1_i + eps1*(HQ_1i)*sign(PnD_1i)
    Q0_i <- Q0_i + eps1*(-HQ_1i)*sign(PnD_1i)
    QA_i <- ifelse(A, Q1_i, Q0_i)
    
    tau_i <- Q1_i - Q0_i
    
    # step 2. update tau_s (think, how to transform/bound)
    H_2i <- 2*tau_s_i
    D_2i <- H_2i*(tau_i - tau_s_i)
    PnD_2i <- mean(D_2i)
    
    tau_s_i <- tau_s_i + eps2*H_2i*sign(PnD_2i)
    
    # update D1 D2, check criteria
    HQ_1i <- 2*(tau_i - tau_s_i)
    H1_1i <- 2*(tau_i - tau_s_i)*(1/gn)
    H0_1i <- 2*(tau_i - tau_s_i)*(-1/(1-gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i)
    D_1i <- HA_1i*(Y - QA_i)
    PnD_1i <- mean(D_1i)
    
    H_2i <- 2*tau_s_i
    D_2i <- H_2i*(tau_i - tau_s_i)
    PnD_2i <- mean(D_2i)
    
    c1 <- abs(PnD_1i) <= sig1/(sqrt(N)*log(N))
    c2 <- abs(PnD_2i) <= sig2/(sqrt(N)*log(N))
    
    if (c1 & c2){
      break
    }
    i <- i + 1
  }
  
  QA_star <- QA_i
  Q1_star <- Q1_i
  Q0_star <- Q0_i
  tau_star <- tau_i
  tau_s_star <- tau_s_i
  
  if(i>=max_it) warning("Max iterations reached in TMLE")
  print(paste0("tmle steps:", i))
  
  # update gamma_s
  suppressWarnings({
    linearUpdate <- glm(tau_star^2 ~ 1,
                        offset = gamma_s_0,
                        family = "gaussian")
  })
  eps3 <- coef(linearUpdate)
  gamma_s_star <- gamma_s_0 + eps3
  
  # calculate Theta_s, scale back
  theta_s_star <- mean(gamma_s_star - tau_s_star^2)
  # ic <- 2*(tau_star - tau_s_star)*(A/gn - (1-A)/(1-gn))*(Y-QA_star) + (tau_star - tau_s_star)^2 - theta_s_star
  # use initial est for ic
  theta_s_0 <- mean(gamma_s_0 - tau_s_0^2)
  ic <- 2*(tau_0 - tau_s_0)*(A/gn - (1-A)/(1-gn))*(Y-QA_0) + (tau_0 - tau_s_0)^2 - theta_s_0
  
  if (logtrans == TRUE){
    ic <- ic/theta_s_star
    se <- sqrt(var(ic)/N)
    ci_l <- exp(log(theta_s_star) - 1.96*se)
    ci_u <- exp(log(theta_s_star) + 1.96*se)
  }else{
    se <- sqrt(var(ic)/N)
    ci_l = theta_s_star - 1.96*se
    ci_u = theta_s_star + 1.96*se
  }
  
  out<- list(
    coef = theta_s_star,
    std_err = se,
    ci_l = ci_l,
    ci_u = ci_u,
    ic = ic
    # df_fit_star = data %>% mutate(mua_hat_star = QA_star,
    #                                 tau_star = tau_star,
    #                                 tau_s_star = tau_s_star,
    #                                 gamma_s_star = gamma_s_star)
  )
  
  # print(out$coef)
  # print(out$ci_l)
  # print(out$ci_u)
  
  return(out)
}




# Parameter: \Psi_s
# Estimator: TMLE, EE and SS
run_VIM2 <- function(df, 
                     sl_Q, 
                     sl_g,
                     sl_x,
                     ws, 
                     cv = FALSE,
                     dr = TRUE,
                     lfm_linear = TRUE, # always true
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
  resTMLE <- TMLE_VIM2(df_fit, max_it, lr)
  
  # ee
  resEE <- EE_VIM2(df_fit)
  
  # ss
  resSS <- SS_VIM2(df_fit)
  
  res <- list('resTMLE' = resTMLE, 
              'resEE' = resEE,
              'resSS' = resSS)
  return(res)
}


# Parameter: VTE
# Estimator: TMLE, EE and SS
run_VTE <- function(df, 
                     sl_Q, 
                     sl_g,
                     sl_x,
                     ws, 
                     cv = FALSE,
                     dr = TRUE,
                     lfm_linear = TRUE, # always true
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
  resTMLE <- TMLE_VTE(df_fit, max_it, lr)
  
  # ee
  resEE <- EE_VTE(df_fit)
  
  # ss
  resSS <- SS_VTE(df_fit)
  
  res <- list('resTMLE' = resTMLE, 
              'resEE' = resEE,
              'resSS' = resSS)
  return(res)
}

# Estimator: SS
SS_VTE <- function(data){
  N <- nrow(data)
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
  
  VTE <- mean((tau_0 - ate)^2)
  ic <- 2*(tau_0 - ate)*(A/gn - (1-A)/(1-gn))*(Y-QA_0) + (tau_0 - ate)^2 - VTE
  se <- sqrt(var(ic)/N)
  
  out<- list(
    coef = VTE,
    std_err = se,
    ci_l = VTE - 1.96*se,
    ci_u = VTE + 1.96*se,
    ic = ic)
  return(out)
}

# Estimator: EE
EE_VTE <- function(data){
  N <- nrow(data)
  po <- data$po
  tau <- data$tau
  tau_s <- data$tau_s
  
  ate <- sum(po)/N
  Sig1 <- sum(po^2)/N - ate^2
  
  # VTE <- mean((po - ate)^2 - (po - tau)^2)
  # ic <- (po - ate)^2 - (po - tau)^2 - VTE
  A <- data$A
  Y <- data$Y
  QA <- data$mua_hat
  gn <- data$pi_hat
  
  VTE <- mean(2*(tau - ate)*(A/gn - (1-A)/(1-gn))*(Y-QA) + (tau - ate)^2)
  ic <- 2*(tau - ate)*(A/gn - (1-A)/(1-gn))*(Y-QA) + (tau - ate)^2 - VTE
  se <- sqrt(var(ic)/N)
  
  out<- list(
    coef = VTE,
    std_err = se,
    ci_l = VTE - 1.96*se,
    ci_u = VTE + 1.96*se,
    ic = ic
  )
  return(out)
}



# Estimator: TMLE
TMLE_VTE <- function(data, max_it = 600, lr = 1e-4){
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
  # theta_s_star <- mean(gamma_s_star - ate_star^2)
  theta_s_star <- mean((tau_star - ate_star)^2)
  # ic <- 2*(tau_star - ate_star)*(A/gn - (1-A)/(1-gn))*(Y-QA_star) + (tau_star - ate_star)^2 - theta_s_star
  # use initial est for ic
  theta_s_0 <- mean((tau_0 - ate_0)^2)
  ic <- 2*(tau_0 - ate_0)*(A/gn - (1-A)/(1-gn))*(Y-QA_0) + (tau_0 - ate_0)^2 - theta_s_0
  se <- sqrt(var(ic)/N)
  
  out<- list(
    coef = theta_s_star,
    std_err = se,
    ci_l = theta_s_star - 1.96*se,
    ci_u = theta_s_star + 1.96*se,
    ic = ic
  )
  return(out)
}

TMLE_VIM2 <- function(data, max_it = 600, lr = 1e-4){
  N <- nrow(data)
  res_vim1 <- TMLE_VIM(data, max_it, lr)
  res_vim2 <- TMLE_VTE(data, max_it, lr)
  
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
  # print(psi1)
  # print(psi2)
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



# helper function
scale01 <- function(x, x_l, x_u){
  return((x - x_l)/(x_u - x_l))
}

rescale <- function(x_scaled, x_l, x_u){
  return(x_scaled*(x_u - x_l) + x_l)
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

# Risk Ratio EY1/EY0
f_log_rr <- function(x, dx) {
  log(x[[2]]) - log(x[[1]])
}

df_log_rr <- function(x, dx) {
  dx[[2]] / x[[2]] - dx[[1]] / x[[1]]
}


# VIM <- function(data, method="AIPW",...){
#   
#   cnames <- c("Y","A","mu1_hat","mu0_hat","pi_hat")
#   if(any(!cnames %in% colnames(data))) stop("Data must contain estimated outcome and propensity score values")
#   
#   ## Use appropriate fitting method
#   if (method == "TMLE_a"){
#     res <- TMLE_VIM_a(data, ...)
#   }else if (method == "TMLE_b"){
#     res <- TMLE_VIM_b(data, ...)
#   }else if (method == "AIPW"){
#     res <- EE_VIM(data, ...)
#   }else{
#     stop("Method not recognized")
#   }
#   
#   res$Method <- method
#   res$call <- match.call(expand.dots = FALSE)
#   class(res) <- "VIM"
#   return(res)
# }
# 
# 
# print.VIM <- function(object){
#   a <- object
#   
#   cat("\nCall:\n", paste(deparse(a$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
#   cat("Estimator:",a$Method)
#   
#   cat("\nEstimate Values:\n")
#   df <- data.frame(Estimate = a$coef,
#                    Std_Error = a$std_err,
#                    CI_l = a$ci_l,
#                    CI_u = a$ci_u)
#   
#   printCoefmat(df, digits = 6, signif.stars = T, na.print = "NA")
# }

