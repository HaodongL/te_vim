library(glmnet)

data = df_fit
y_l = 0
y_u = 1
max.it= 1e4
lfm_linear = T

TMLE_VIM <- function(data, y_l, y_u, max.it=600, lfm_linear = FALSE){
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
  
  # sig1, sig2
  sig1 <- sd(2*(tau_0 - tau_s_0)*((2*A-1)/gn)*(Y - QA_0))
  sig2 <- sd(2*tau_s_0*(tau_0 - tau_s_0))
  
  # eps1,eps2
  eps1 <- 1e-4
  eps2 <- 1e-4
  
  # i 
  i <- 1
  QA_i <- QA_0
  Q1_i <- Q1_0
  Q0_i <- Q0_0
  tau_i <- tau_0
  tau_s_i <- tau_s_0
  
  # logging <- matrix(NA, max.it, 3)
  
  while(i <= max.it){
    # step 1. update Q and tau
    if (i == 1){
      H1_1i <- 2*(tau_i - tau_s_i)*((2*1-1)/gn)
      H0_1i <- 2*(tau_i - tau_s_i)*((2*0-1)/gn)
      HA_1i <- ifelse(A, H1_1i, H0_1i)
      D_1i <- HA_1i*(Y - QA_i)
      PnD_1i <- mean(D_1i)
    }
    
    if (lfm_linear){
      # X <- matrix(1, N, 2)
      # update_fit1 <- cv.glmnet(x = X, 
      #                          y = Y, 
      #                          offset = QA_i, 
      #                          weights = HA_1i, 
      #                          intercept = FALSE, 
      #                          family = "gaussian", 
      #                          lower.limits = -eps1, 
      #                          upper.limits = eps1)
      # eps1_n <- coef(update_fit1, s = "lambda.min")[1]
      eps1_n <- eps1
      loss1_a <- mean(HA_1i*(Y - QA_i)^2)
      loss1_b <- mean(HA_1i*(Y - (QA_i + eps1_n))^2)
      eps1_n <- ifelse(loss1_a >= loss1_b, eps1_n, -eps1_n)
      Q1_i <- Q1_i + eps1_n
      Q0_i <- Q0_i + eps1_n
    }else{
      Q1_i <- plogis(qlogis(Q1_i) + eps1*H1_1i*sign(PnD_1i))
      Q0_i <- plogis(qlogis(Q0_i) + eps1*H0_1i*sign(PnD_1i))
    }
    QA_i <- ifelse(A, Q1_i, Q0_i)
    tau_i <- Q1_i - Q0_i
    
    # step 2. update tau_s (think, how to transform/bound)
    H_2i <- 2*tau_s_i
    D_2i <- H_2i*(tau_i - tau_s_i)
    PnD_2i <- mean(D_2i)
    
    # logging[i, 1] <- i
    # logging[i, 2] <- PnD_1i
    # logging[i, 3] <- PnD_2i
    
    # tau_s_i <- plogis(qlogis(tau_s_i) - eps2*H_2i*sign(PnD_2i))
    if (lfm_linear){
      # X <- rep(1, N)
      # update_fit2 <- cv.glmnet(x = X, y = tau_i, offset = tau_s_i, weights = H_2i, intercept = FALSE, family = "gaussian", lower.limits = -eps2, upper.limits = eps2)
      # eps2_n <- coef(update_fit2, s = "lambda.min")[1]
      eps2_n <- eps2
      loss2_a <- sum(H_2i*(tau_i - tau_s_i)^2)
      loss2_b <- sum(H_2i*(tau_i - (tau_s_i + eps2_n))^2)
      eps2_n <- ifelse(loss2_a >= loss2_b, eps2_n, -eps2_n)
      tau_s_i <- tau_s_i + eps2_n
    }else{
      maxabs <- max(abs(tau_s_i))
      eps2_new <- min(0.5*(1/maxabs - 1), eps2)
      tau_s_i <- tau_s_i + eps2_new*H_2i*sign(PnD_2i)
    }
    
    
    # update D1 D2, check creteria
    H1_1i <- 2*(tau_i - tau_s_i)*((2*1-1)/gn)
    H0_1i <- 2*(tau_i - tau_s_i)*((2*0-1)/gn)
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
    print(eps1_n)
    print(abs(PnD_1i))
    print(abs(PnD_2i))
  }
  
  QA_star <- QA_i
  Q1_star <- Q1_i
  Q0_star <- Q0_i
  tau_star <- tau_i
  tau_s_star <- tau_s_i
  
  # if(i>=max.it) stop("Max iterations reached in TMLE")
  
  if(i>=max.it) warning("Max iterations reached in TMLE")
  
  # update gamma_s
  if (lfm_linear){
    suppressWarnings({
      linearUpdate <- glm(tau_star^2 ~ 1, 
                          offset = qlogis(gamma_s_0),
                          family = "gaussian")
    })
    eps3 <- coef(linearUpdate)
    gamma_s_star <- gamma_s_0 + eps3
  }else{
    suppressWarnings({
      logitUpdate <- glm(tau_star^2 ~ 1, 
                         offset = qlogis(gamma_s_0),
                         family = "quasibinomial")
    })
    eps3 <- coef(logitUpdate)
    gamma_s_star <- plogis(qlogis(gamma_s_0) + eps3)
  }
  
  
  # calculate Theta_s, scale back
  if (lfm_linear){
    y_u = 1
    y_l = 0
  }
  theta_s_star <- mean(gamma_s_star - tau_s_star^2)*(y_u-y_l)^2
  # theta_s_star <- mean((tau_star - tau_s_star)^2)*(y_u-y_l)^2
  ic <- 2*(tau_star - tau_s_star)*((2*A-1)/gn)*(Y-QA_star) + (tau_star - tau_s_star)^2 - theta_s_star
  ss <- sqrt(var(ic)/N)*(y_u-y_l)^2
  
  coef <- theta_s_star
  std_err <- ss
  names(coef) <- names(std_err) <- c("VIM_Theta_s")
  
  out<- list(
    coef = coef,
    std_err = std_err,
    ci_l = coef - 1.96*std_err,
    ci_u = coef + 1.96*std_err
  )
  
  return(out)
}

# Estimator: TMLE and EE
run_VIM_Theta <- function(df, 
                          sl_Q, 
                          sl_g,
                          sl_x,
                          ws, 
                          cv = TRUE,
                          dr = TRUE,
                          lfm_linear = FALSE,
                          max.it=600, 
                          Q_bounds = c(0.001, 0.999), 
                          g_bounds = c(0.025, 0.975),
                          tau_bounds = c(-1+1e-3, 1-1e-3),
                          tau_s_bounds = c(-1+1e-3, 1-1e-3),
                          gamma_s_bounds = c(1e-6, 1-1e-6)){
  # transform Y into [0,1]
  y_l <- min(df$Y)
  y_u <- max(df$Y)
  if (!lfm_linear){
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
  resTMLE <- TMLE_VIM(df_fit, y_l, y_u, max.it, lfm_linear)
  
  # ee_vim
  if (!lfm_linear){
    df_fit_sl$Y <- rescale(df_fit_sl$Y, y_l, y_u)
    df_fit_sl$tau <- df_fit_sl$tau*(y_u - y_l)
    df_fit_sl$tau_s <- df_fit_sl$tau_s*(y_u - y_l)
    df_fit_sl$po <- df_fit_sl$po*(y_u - y_l)
  }
  resEE <- EE_VIM(df_fit)
  
  res <- list('resTMLE' = resTMLE, 
              'resEE' = resEE)
  return(res)
}
