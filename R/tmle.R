#Function for computing ATE and VTE using TMLE
#
#data is a dataframe or tibble containing:
#Y - outcome of interest
#A - treatment indicator
#pi_hat - estimate of the propensity score E(A|X)
#mu1_hat - estimate of E(Y|A=1,X)
#mu0_hat - estimate of E(Y|A=0,X)
#
#ab - (optional) Min and Max of range of Y. Defaults to min/max of (Y,mu0_hat,mu1_hat)

TMLE_VTE <- function(data,ab=NULL){
  ## Some params
  max.it <- 600 #maximum number of iterations in targetting step
  eps <- 0.0001 #TMLE target step size
  
  a <- with(data,min(Y,mu1_hat,mu0_hat))
  b <- with(data,max(Y,mu1_hat,mu0_hat)) - a
  
  if(!is.null(ab)){
    if(ab[1]<=a & ab[2]>=(a+b)){
      a <- ab[1]
      b <- ab[2]-ab[1]
    }else warning("ab values incorrect - using default method")
  }

  N = NROW(data)
  A <- data$A
  beta1 <- with(data,  1/pi_hat )
  beta0 <- with(data,  1/(pi_hat-1) )
  beta <- ifelse(A,beta1,beta0)
  y <- (data$Y -a)/b

  expitq1 <- qlogis((data$mu1_hat-a)/b)
  expitq0 <- qlogis((data$mu0_hat-a)/b)
  
  retarget <- TRUE
  it <- 1
  while(retarget & it <= max.it){
    q_1 <- plogis(expitq1)
    q_0 <- plogis(expitq0)
    cate <- q_1 - q_0
    #q_m   <- q_0 + A*cate
    #loss  <- -sum(y*log(q_m) + (1-y)*log(1-q_m))/N  
    po <- beta*(y-q_0 - A*cate) + cate
    
    a0 <- sum(cate)/N
    a1 <- sum(po)/N
    b0 <- sum(cate^2)/N
    b1 <- sum(po^2)/N
    b2 <- sum(po*cate)/N
    
    ATE <- a0
    VTE <- b0 - a0^2
    pD1 <- a1 - a0
    pD2 <- 2*(b2-b0+a0^2-a0*a1)
    pD_norm <- sqrt(pD1^2 + pD2^2)
    
    Sig1 <- b1 - 2*a1*a0 + a0^2
    Sig2 <- sum(((po-a0)^2 -(po-cate)^2 - VTE)^2  )/N 
    
    if( (N^2*pD1^2 >= Sig1) | (N^2*pD2^2 >= Sig2)){
      expitq1 <- expitq1  + eps*beta1*(pD1+2*(cate-ATE)*pD2)/pD_norm
      expitq0 <- expitq0  + eps*beta0*(pD1+2*(cate-ATE)*pD2)/pD_norm 
      
    }else{
      retarget <- FALSE
    }
    it <- it+1

  }
  if(it>=max.it) warning("Max iterations reached in TMLE")
  
  # print(paste0("it = ", it))
  print(paste0("eic1 = ", N^2*pD1^2))
  print(paste0("cutoff1 = ", Sig1))
  print(paste0("eic2 = ", N^2*pD2^2))
  print(paste0("cutoff2 = ", Sig2))

  
  rootV <- ifelse(VTE>=0,sqrt(VTE),NA)
  ss <- sqrt(Sig2/N)

  coef <- c(ATE*b,
            rootV*b,
            VTE*b*b)
  std.err <- c(sqrt(Sig1/N)*b,
               ss*b/rootV,
               ss*b*b)
  names(coef) <- names(std.err) <- c("ATE","rootVTE","VTE")
  names(cate) <- NULL
  
  out<- list(
    coef = coef,
    std.err = std.err,
    CATE = cate*b
  )
  
  return(out)
}



# 
# TMLE_VIMD <- function(df, df_fit, option, covar){
#   ## Some params
#   max.it <- 600 #maximum number of iterations in targeting step
#   eps <- 0.0001 #TMLE target step size
#   
#   y_l <- with(df_fit,min(Y,mu1_hat,mu0_hat))
#   y_u <- with(df_fit,max(Y,mu1_hat,mu0_hat))
#   
#   # Q,Y,g,A
#   n = NROW(df_fit)
#   A <- df_fit$A
#   g_hat <- df_fit$pi_hat
#   Y_t <- (df_fit$Y - y_l)/(y_u - y_l)
#   Q_t <- (df_fit$mua_hat - y_l)/(y_u - y_l)
#   Q1_t <- (df_fit$mu1_hat - y_l)/(y_u - y_l)
#   Q0_t <- (df_fit$mu0_hat - y_l)/(y_u - y_l)
#   
#   # cate
#   tau <- df_fit$CATE
#   tau_s <- df_fit$tau_s
#   
#   # psi and ic
#   psi <- mean((tau - tau_s)^2)
#   HAW = 2*(tau - tau_s)*(2*A - 1)/(g_hat) 
#   H1W = 2*(tau - tau_s)*(2*1 - 1)/(g_hat)
#   H0W = 2*(tau - tau_s)*(2* - 1)/(g_hat)
#   D = HAW*(Y_t - Q_t) + (tau - tau_s)^2 - psi
#   
#   sig_hat <- sd(D)
#   criterion <- abs(mean(D)) < sig_hat/(sqrt(n)*log(n))
#   
#   it <- 1
#   while(it <= max.it){
#     if(criterion){
#       break
#     }else{
#       
#       # update Q
#       Q1_t <- qlogis(plogis(Q1_t) - eps*H1W*sign(D))
#       Q0_t <- qlogis(plogis(Q0_t) - eps*H0W*sign(D))
#       Q_t <- ifelse(A, Q1_t, Q0_t)
#       
#       # update tau, tau_s
#       df_fit <- df_fit %>% select('Y', 'A', "pi_hat",  "mu1_hat", "mu0_hat", "mua_hat")
#       df_fit$mu1_hat <- Q1_t
#       df_fit$mu0_hat <- Q0_t
#       df_fit$mua_hat <- Q_t
#       
#       df_fit <- fit_tau(df = df, df_fit = df_fit, option = option)
#       df_fit <- fit_tau_s(df = df, df_fit = df_fit, covar = covar)
#       
#       tau <- df_fit$CATE
#       tau_s <- df_fit$tau_s
#       
#       # update psi, ic
#       psi <- mean((tau - tau_s)^2)
#       HAW = 2*(tau - tau_s)*(2*A - 1)/(g_hat) 
#       H1W = 2*(tau - tau_s)*(2*1 - 1)/(g_hat)
#       H0W = 2*(tau - tau_s)*(2* - 1)/(g_hat)
#       D = HAW*(Y_t - Q_t) + (tau - tau_s)^2 - psi
#       
#       # re-evaluate criterion
#       criterion <- abs(mean(D)) < sig_hat/(sqrt(n)*log(n))
#       it <- it + 1
#     }
#   }
#   if(it>=max.it) warning("Max iterations reached in TMLE")
#   
#   # scale back
#   
#   psi <- c(psi*(y_u - y_l)^2)
#   se <- sd(D)*(y_u - y_l)^2/sqrt(n)
#   ci_l <- psi - 1.96*se
#   ci_u <- psi + 1.96*se
#   
#   out<- list(
#     psi = psi,
#     ci_l = ci_l,
#     ci_u = ci_u
#   )
#   
#   return(out)
# }

