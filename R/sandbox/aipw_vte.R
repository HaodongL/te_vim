#Function for computing ATE and VTE using one-step estimators
#
#data is a dataframe or tibble containing:
#Y - outcome of interest
#A - treatment indicator
#pi_hat - estimate of the propensity score E(A|X)
#mu1_hat - estimate of E(Y|A=1,X)
#mu0_hat - estimate of E(Y|A=0,X)
#CATE - (optional) Separate estimate of CATE. Defaults to CATE = mu1_hat-mu0_hat

AIPW_VTE <- function(data){
  N <- NROW(data)
  
  if("CATE" %in% names(data)) {
    CATE <- data$CATE
  } else{
    CATE <- with(data,mu1_hat - mu0_hat)  
  }
  PO <- with(data,(Y-mu0_hat-A*(mu1_hat - mu0_hat))*(A-pi_hat)/(pi_hat*(1-pi_hat))+ mu1_hat - mu0_hat)
  
  ATE <- sum(PO)/N
  Sig1 <- sum(PO^2)/N - ATE^2
  x <- PO-CATE
  VTE <- Sig1 - sum(x^2)/N
  Sig2 <- sum((CATE-ATE)^2*(CATE-ATE+2*x)^2)/N -VTE^2
  
  rootV <- ifelse(VTE>=0,sqrt(VTE),NA)
  ss <- sqrt(Sig2/N)

  coef <- c(ATE,rootV,VTE)
  std.err <- c(sqrt(Sig1/N),ss/rootV,ss)
  names(coef) <- names(std.err) <- c("ATE","rootVTE","VTE")
  
  out<- list(
    coef = coef,
    std.err = std.err,
    CATE = CATE
  )
  
  return(out)
}


