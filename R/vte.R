VTE <- function(data,method="TMLE",...){

  cnames <- c("Y","A","mu1_hat","mu0_hat","pi_hat")
  if(any(!cnames %in% colnames(data))) stop("Data must contain estimated outcome and propensity score values")

  ## Use appropriate fitting method
  if (method=="TMLE"){
    a <- TMLE_VTE(data,...)
  }else if (method == "AIPW"){
    a <- AIPW_VTE(data,...)
  }else{
    stop("Method not recognized")
  }
  
  a$Method <- method
  a$call <- match.call(expand.dots = FALSE)
  class(a) <- "VTE"
  a
}


print.VTE <- function(object){
  a <- object
  
  cat("\nCall:\n", paste(deparse(a$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Estimator:",a$Method)
  
  cat("\nEstimand Values:\n")
  
  wald <- (a$coef/a$std.err)^2
  df <- data.frame(Estimate = a$coef, Std.Error = a$std.err ,
                   Wald.value = wald,
                   Wald.pval = pchisq(wald,df=1,lower.tail = F) )
  
  printCoefmat(df, digits = 6, signif.stars = T,na.print = "NA",
               tst.ind = 3,P.values=T,has.Pvalue=T)
}

