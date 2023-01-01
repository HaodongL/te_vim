TMLE_VIM_a <- function(data, y_l, y_u, max.it = 600, lr = 1e-4){
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

  # logging <- matrix(NA, max.it, 3)

  while(i <= max.it){
    # step 1. update Q and tau
    if (i == 1){
      H1_1i <- 2*(tau_i - tau_s_i)*(1/gn)
      H0_1i <- 2*(tau_i - tau_s_i)*(-1/(1-gn))
      HA_1i <- ifelse(A, H1_1i, H0_1i)
      D_1i <- HA_1i*(Y - QA_i)
      PnD_1i <- mean(D_1i)
    }

    # Q1_i <- plogis(qlogis(Q1_i) + eps1*H1_1i*sign(PnD_1i))
    # Q0_i <- plogis(qlogis(Q0_i) + eps1*H0_1i*sign(PnD_1i))
    # QA_i <- ifelse(A, Q1_i, Q0_i)

    Q1_i <- Q1_i + eps1*H1_1i*sign(PnD_1i)
    Q0_i <- Q0_i + eps1*H0_1i*sign(PnD_1i)
    QA_i <- ifelse(A, Q1_i, Q0_i)

    tau_i <- Q1_i - Q0_i

    # step 2. update tau_s (think, how to transform/bound)
    H_2i <- 2*tau_s_i
    D_2i <- H_2i*(tau_i - tau_s_i)
    PnD_2i <- mean(D_2i)

    # maxabs <- max(abs(tau_s_i))
    # eps2_new <- min(0.5*(1/maxabs - 1), eps2)
    # tau_s_i <- tau_s_i + eps2_new*H_2i*sign(PnD_2i)
    tau_s_i <- tau_s_i + eps2*H_2i*sign(PnD_2i)

    # update D1 D2, check creteria
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

  if(i>=max.it) warning("Max iterations reached in TMLE")

  # update gamma_s
  # suppressWarnings({
  #   logitUpdate <- glm(tau_star^2 ~ 1,
  #                      offset = qlogis(gamma_s_0),
  #                      family = "quasibinomial")
  # })
  # eps3 <- coef(logitUpdate)
  # gamma_s_star <- plogis(qlogis(gamma_s_0) + eps3)
  suppressWarnings({
    linearUpdate <- glm(tau_star^2 ~ 1,
                       offset = gamma_s_0,
                       family = "gaussian")
  })
  eps3 <- coef(linearUpdate)
  gamma_s_star <- gamma_s_0 + eps3

  # calculate Theta_s, scale back
  theta_s_star <- mean(gamma_s_star - tau_s_star^2)*(y_u-y_l)^2
  # theta_s_star <- mean((tau_star - tau_s_star)^2)*(y_u-y_l)^2
  ic <- 2*(tau_star - tau_s_star)*(A/gn - (1-A)/(1-gn))*(Y-QA_star) + (tau_star - tau_s_star)^2 - theta_s_star
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


TMLE_VIM_b <- function(data, y_l, y_u, max.it = 600, lr = 1e-4){
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

  # logging <- matrix(NA, max.it, 3)

  po <- (Y - QA_i)*(A/gn - (1-A)/(1-gn)) + Q1_i - Q0_i
  H_1i <- 2*(tau_i - tau_s_i)
  D_1i <- H_1i*(po - tau_i)
  PnD_1i <- mean(D_1i)
  a <- -1
  b <- 1

  while(i <= max.it){

    # step 1. update tau and Q
    # tau_i_scale <- (tau_i - a)/(b - a)
    # H_1i_scale <- H_1i/(b - a)
    # tau_i_scale <- plogis(qlogis(tau_i_scale) + eps1*H_1i_scale*sign(PnD_1i))
    # tau_i <- tau_i_scale * (b - a) + a

    tau_i <- tau_i + eps1*H_1i*sign(PnD_1i)
    Q1_i <- tau_i + Q0_0
    QA_i[which(A == 1)] <- Q1_i[which(A == 1)]

    # step 2. update tau_s (think, how to transform/bound)
    H_2i <- 2*tau_s_i
    D_2i <- H_2i*(tau_i - tau_s_i)
    PnD_2i <- mean(D_2i)
    # maxabs <- max(abs(tau_s_i))
    # eps2_new <- min(0.5*(1/maxabs - 1), eps2)
    # tau_s_i <- tau_s_i + eps2_new*H_2i*sign(PnD_2i)
    tau_s_i <- tau_s_i + eps2*H_2i*sign(PnD_2i)

    # update D1 D2, check creteria
    po <- (Y - QA_i)*(A/gn - (1-A)/(1-gn)) + Q1_i - Q0_i
    H_1i <- 2*(tau_i - tau_s_i)
    D_1i <- H_1i*(po - tau_i)
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

  if(i>=max.it) warning("Max iterations reached in TMLE")

  # update gamma_s
  # suppressWarnings({
  #   logitUpdate <- glm(tau_star^2 ~ 1,
  #                      offset = qlogis(gamma_s_0),
  #                      family = "quasibinomial")
  # })
  # eps3 <- coef(logitUpdate)
  # gamma_s_star <- plogis(qlogis(gamma_s_0) + eps3)
  suppressWarnings({
    linearUpdate <- glm(tau_star^2 ~ 1,
                       offset = gamma_s_0,
                       family = "gaussian")
  })
  eps3 <- coef(linearUpdate)
  gamma_s_star <- gamma_s_0 + eps3

  # calculate Theta_s, scale back
  theta_s_star <- mean(gamma_s_star - tau_s_star^2)*(y_u-y_l)^2
  # theta_s_star <- mean((tau_star - tau_s_star)^2)*(y_u-y_l)^2
  ic <- 2*(tau_star - tau_s_star)*(A/gn - (1-A)/(1-gn))*(Y-QA_star) + (tau_star - tau_s_star)^2 - theta_s_star
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

