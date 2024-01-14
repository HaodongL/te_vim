# Target VTE and VIM simultaneously

TMLE_VIMnVTE <- function(data, max_it = 600, lr = 1e-4){
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
  ate_0 <- rep(mean(tau_0), N)
  
  # sig1, sig2
  sig1 <- sd(2*(tau_0 - tau_s_0)*(A/gn - (1-A)/(1-gn))*(Y - QA_0))
  sig2 <- sd(2*tau_s_0*(tau_0 - tau_s_0))
  
  # sig1 VTE
  sig1_vte <- sd(2*(tau_0 - ate_0)*(A/gn - (1-A)/(1-gn))*(Y - QA_0))
  
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
  ate_i <- ate_0
  
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
      
      HQ_1i_vte <- 2*(tau_i - ate_i)
      H1_1i_vte <- 2*(tau_i - ate_i)*(1/gn)
      H0_1i_vte <- 2*(tau_i - ate_i)*(-1/(1-gn))
      HA_1i_vte <- ifelse(A, H1_1i_vte, H0_1i_vte)
      D_1i_vte <- HA_1i_vte*(Y - QA_i)
      PnD_1i_vte <- mean(D_1i_vte)
    }
    
    norm_pd <- sqrt(PnD_1i^2 + PnD_1i_vte^2)
    Q1_i <- Q1_i + eps1*(HQ_1i*PnD_1i+ HQ_1i_vte*PnD_1i_vte)/norm_pd 
    Q0_i <- Q0_i + eps1*(-HQ_1i*PnD_1i- HQ_1i_vte*PnD_1i_vte)/norm_pd 
    QA_i <- ifelse(A, Q1_i, Q0_i)
  
    tau_i <- Q1_i - Q0_i
    ate_i <- rep(mean(tau_i), N)
    
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
    
    # update D1 D2, check criteria VTE
    HQ_1i_vte <- 2*(tau_i - ate_i)
    H1_1i_vte <- 2*(tau_i - ate_i)*(1/gn)
    H0_1i_vte <- 2*(tau_i - ate_i)*(-1/(1-gn))
    HA_1i_vte <- ifelse(A, H1_1i_vte, H0_1i_vte)
    D_1i_vte <- HA_1i_vte*(Y - QA_i)
    PnD_1i_vte <- mean(D_1i_vte)
    
    c1_vte <- abs(PnD_1i_vte) <= sig1_vte/(sqrt(N)*log(N))
    
    if (c1 & c2 & c1_vte){
      break
    }
    i <- i + 1
  }
  
  QA_star <- QA_i
  Q1_star <- Q1_i
  Q0_star <- Q0_i
  tau_star <- tau_i
  tau_s_star <- tau_s_i
  ate_star <- rep(mean(tau_i), N)
  
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
  
  se <- sqrt(var(ic)/N)
  ci_l = theta_s_star - 1.96*se
  ci_u = theta_s_star + 1.96*se
  
  # vte
  vte_star <- mean((tau_star - ate_star)^2)
  vte_0 <- mean((tau_0 - ate_0)^2)
  ic_vte <- 2*(tau_0 - ate_0)*(A/gn - (1-A)/(1-gn))*(Y-QA_0) + (tau_0 - ate_0)^2 - vte_0
  se_vte <- sqrt(var(ic_vte)/N)
  ci_l_vte = vte_star - 1.96*se_vte
  ci_u_vte = vte_star + 1.96*se_vte
  
  out<- list(
    "VIM" = list(coef = theta_s_star,
                 std_err = se,
                 ci_l = ci_l,
                 ci_u = ci_u),
    "VTE" = list(coef = vte_star,
                 std_err = se_vte,
                 ci_l = ci_l_vte,
                 ci_u = ci_u_vte)
  )
  
  return(out)
}


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
  sig2 <- sd((A/gn - (1-A)/(1-gn))*(Y - QA_0))
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
      # HQ_1i <- 2*(tau_i - ate_i)
      H1_1i <- 2*(tau_i - ate_i)*(1/gn)
      H0_1i <- 2*(tau_i - ate_i)*(-1/(1-gn))
      HA_1i <- ifelse(A, H1_1i, H0_1i)
      D_1i <- HA_1i*(Y - QA_i)
      PnD_1i <- mean(D_1i)
      
      # HQ_1i_ate <- 2*(tau_i - ate_i)
      H1_1i_ate <- 1/gn
      H0_1i_ate <- -1/(1-gn)
      HA_1i_ate <- ifelse(A, H1_1i_ate, H0_1i_ate)
      D_1i_ate <- HA_1i_ate*(Y - QA_i)
      PnD_1i_ate <- mean(D_1i_ate)
    }
    
    # Q1_i <- Q1_i + eps1*(HQ_1i)*sign(PnD_1i)
    # Q0_i <- Q0_i + eps1*(-HQ_1i)*sign(PnD_1i)
    # QA_i <- ifelse(A, Q1_i, Q0_i)
    # Q1_i <- plogis(qlogis(Q1_i) + eps1*H1_1i*sign(PnD_1i))
    # Q0_i <- plogis(qlogis(Q0_i) + eps1*H0_1i*sign(PnD_1i))
    norm_pd <- sqrt(PnD_1i^2 + PnD_1i_ate^2)
    Q1_i <- Q1_i + eps1*(H1_1i*PnD_1i+ H1_1i_ate*PnD_1i_ate)/norm_pd 
    Q0_i <- Q0_i + eps1*(H0_1i*PnD_1i + H0_1i_ate*PnD_1i_ate)/norm_pd 
    QA_i <- ifelse(A, Q1_i, Q0_i)
    
    tau_i <- Q1_i - Q0_i
    ate_i <- rep(mean(tau_i), N)
    
    # update D1 D2, check criteria
    # HQ_1i <- 2*(tau_i - ate_i)
    H1_1i <- 2*(tau_i - ate_i)*(1/gn)
    H0_1i <- 2*(tau_i - ate_i)*(-1/(1-gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i)
    D_1i <- HA_1i*(Y - QA_i)
    PnD_1i <- mean(D_1i)
    
    H1_1i_ate <- 1/gn
    H0_1i_ate <- -1/(1-gn)
    HA_1i_ate <- ifelse(A, H1_1i_ate, H0_1i_ate)
    D_1i_ate <- HA_1i_ate*(Y - QA_i)
    PnD_1i_ate <- mean(D_1i_ate)
    
    c1 <- abs(PnD_1i) <= sig1/(sqrt(N)*log(N))
    c2 <- abs(PnD_1i_ate) <= sig2/(sqrt(N)*log(N))
    
    # c1 <- abs(PnD_1i) <= sqrt(sig1)/N  # sig1/(sqrt(N)*log(N))
    # c2 <- abs(PnD_1i_ate) <= sqrt(sig2)/N # sig2/(sqrt(N)*log(N))
    
    if (c1 & c2){
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
  sig2 <- sd((A/gn - (1-A)/(1-gn))*(Y - QA_0))
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
      
      HQ_1i_ate <- 1
      H1_1i_ate <- 1/gn
      H0_1i_ate <- -1/(1-gn)
      HA_1i_ate <- ifelse(A, H1_1i_ate, H0_1i_ate)
      D_1i_ate <- HA_1i_ate*(Y - QA_i)
      PnD_1i_ate <- mean(D_1i_ate)
    }
    
    # Q1_i <- Q1_i + eps1*(HQ_1i)*sign(PnD_1i)
    # Q0_i <- Q0_i + eps1*(-HQ_1i)*sign(PnD_1i)
    # QA_i <- ifelse(A, Q1_i, Q0_i)
    # Q1_i <- plogis(qlogis(Q1_i) + eps1*H1_1i*sign(PnD_1i))
    # Q0_i <- plogis(qlogis(Q0_i) + eps1*H0_1i*sign(PnD_1i))
    norm_pd <- sqrt(PnD_1i^2 + PnD_1i_ate^2)
    Q1_i <- Q1_i + eps1*(HQ_1i*PnD_1i+ HQ_1i_ate*PnD_1i_ate)/norm_pd 
    Q0_i <- Q0_i - eps1*(HQ_1i*PnD_1i + HQ_1i_ate*PnD_1i_ate)/norm_pd 
    QA_i <- ifelse(A, Q1_i, Q0_i)
    
    tau_i <- Q1_i - Q0_i
    ate_i <- rep(mean(tau_i), N)
    
    # update D1 D2, check criteria
    # HQ_1i <- 2*(tau_i - ate_i)
    H1_1i <- 2*(tau_i - ate_i)*(1/gn)
    H0_1i <- 2*(tau_i - ate_i)*(-1/(1-gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i)
    D_1i <- HA_1i*(Y - QA_i)
    PnD_1i <- mean(D_1i)
    
    H1_1i_ate <- 1/gn
    H0_1i_ate <- -1/(1-gn)
    HA_1i_ate <- ifelse(A, H1_1i_ate, H0_1i_ate)
    D_1i_ate <- HA_1i_ate*(Y - QA_i)
    PnD_1i_ate <- mean(D_1i_ate)
    
    c1 <- abs(PnD_1i) <= sig1/(sqrt(N)*log(N))
    c2 <- abs(PnD_1i_ate) <= sig2/(sqrt(N)*log(N))
    
    # c1 <- abs(PnD_1i) <= sqrt(sig1)/N  # sig1/(sqrt(N)*log(N))
    # c2 <- abs(PnD_1i_ate) <= sqrt(sig2)/N # sig2/(sqrt(N)*log(N))
    
    if (c1 & c2){
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
  sig2 <- sd((A/gn - (1-A)/(1-gn))*(Y - QA_0))
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
      
      HQ_1i_ate <- 1
      H1_1i_ate <- 1/gn
      H0_1i_ate <- -1/(1-gn)
      HA_1i_ate <- ifelse(A, H1_1i_ate, H0_1i_ate)
      D_1i_ate <- HA_1i_ate*(Y - QA_i)
      PnD_1i_ate <- mean(D_1i_ate)
    }
    
    # Q1_i <- Q1_i + eps1*(HQ_1i)*sign(PnD_1i)
    # Q0_i <- Q0_i + eps1*(-HQ_1i)*sign(PnD_1i)
    # QA_i <- ifelse(A, Q1_i, Q0_i)
    # Q1_i <- plogis(qlogis(Q1_i) + eps1*H1_1i*sign(PnD_1i))
    # Q0_i <- plogis(qlogis(Q0_i) + eps1*H0_1i*sign(PnD_1i))
    norm_pd <- sqrt(PnD_1i^2 + PnD_1i_ate^2)
    # Q1_i <- Q1_i + eps1*(HQ_1i*PnD_1i+ HQ_1i_ate*PnD_1i_ate)/norm_pd
    # Q0_i <- Q0_i - eps1*(HQ_1i*PnD_1i + HQ_1i_ate*PnD_1i_ate)/norm_pd

    Q1_i <- plogis(qlogis(Q1_i) + eps1*(H1_1i*PnD_1i+ H1_1i_ate*PnD_1i_ate)/norm_pd)
    Q0_i <- plogis(qlogis(Q0_i) + eps1*(H0_1i*PnD_1i + H0_1i_ate*PnD_1i_ate)/norm_pd)
    
    QA_i <- ifelse(A, Q1_i, Q0_i)
    
    tau_i <- Q1_i - Q0_i
    ate_i <- rep(mean(tau_i), N)
    
    # update D1 D2, check criteria
    # HQ_1i <- 2*(tau_i - ate_i)
    H1_1i <- 2*(tau_i - ate_i)*(1/gn)
    H0_1i <- 2*(tau_i - ate_i)*(-1/(1-gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i)
    D_1i <- HA_1i*(Y - QA_i)
    PnD_1i <- mean(D_1i)
    
    H1_1i_ate <- 1/gn
    H0_1i_ate <- -1/(1-gn)
    HA_1i_ate <- ifelse(A, H1_1i_ate, H0_1i_ate)
    D_1i_ate <- HA_1i_ate*(Y - QA_i)
    PnD_1i_ate <- mean(D_1i_ate)
    
    c1 <- abs(PnD_1i) <= sig1/(sqrt(N)*log(N))
    c2 <- abs(PnD_1i_ate) <= sig2/(sqrt(N)*log(N))
    
    # c1 <- abs(PnD_1i) <= sqrt(sig1)/N  # sig1/(sqrt(N)*log(N))
    # c2 <- abs(PnD_1i_ate) <= sqrt(sig2)/N # sig2/(sqrt(N)*log(N))
    
    if (c1 & c2){
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
