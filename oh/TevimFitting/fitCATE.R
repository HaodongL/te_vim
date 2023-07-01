

DRLearn_SL <- function(df,X,SL.library,SLfolds=10,SL_type = "discrete"){
  ind <- df$inTrain
  
  if(SL_type =="discrete"){
    PO <- with(df,(A-pi_hat)*(Y-A*mu1_hat - (1-A)*mu0_hat)/(pi_hat*(1-pi_hat)) +  mu1_hat-mu0_hat  )
  }else if(SL_type =="non_discrete"){
    PO <- with(df,(A-pi_hat_SL)*(Y-A*mu1_hat_SL - (1-A)*mu0_hat_SL)/(pi_hat*(1-pi_hat_SL)) +  mu1_hat_SL-mu0_hat_SL  )
  }else stop("SL_type not recognised")

  mod <- SuperLearner::SuperLearner(PO[ind],X[ind,],newX=X,
                                    family = gaussian(),
                                    SL.library = SL.library,
                                    cvControl = list(V=SLfolds))
  
  if(SL_type =="discrete"){
    CATE <- mod$library.predict[, which.min(mod$cvRisk)] #discrete superlearner
  }else if(SL_type =="non_discrete"){
    CATE <- mod$SL.predict[,1] 
  }
  
  out<- tibble(
    PO=PO,
    CATE=CATE,
    ATE = sum(ind*CATE)/sum(ind),
    inTrain = ind
  ) 
  return(out)
}

TLearn_SL <- function(df,SL_type = "discrete"){
  ind <- df$inTrain
  
  if(SL_type =="discrete"){
    PO   <- with(df,(A-pi_hat)*(Y-A*mu1_hat - (1-A)*mu0_hat)/(pi_hat*(1-pi_hat)) +  mu1_hat-mu0_hat  )
    CATE <- with(df,mu1_hat-mu0_hat)
  }else if(SL_type =="non_discrete"){
    PO   <- with(df,(A-pi_hat_SL)*(Y-A*mu1_hat_SL - (1-A)*mu0_hat_SL)/(pi_hat*(1-pi_hat_SL)) +  mu1_hat_SL-mu0_hat_SL)
    CATE <- with(df,mu1_hat_SL-mu0_hat_SL)
  }else stop("SL_type not recognised")
  
  out<- tibble(
    PO=PO,
    CATE=CATE,
    ATE = sum(ind*CATE)/sum(ind),
    inTrain = ind
  ) 
  return(out)
}


marginal_CATEs <- function(df,X,cov,SL.library,SLfolds=10,SL_type = "discrete"){
  ind <- df$inTrain
  X<- select(X,!matches(cov))
  
  if(all(ind)){
    newX <- X
  }else{
    newX <- X[!ind,]
  }
  
  mod <- SuperLearner::SuperLearner(df$CATE[ind],X[ind,],newX=newX,
                                    family = gaussian(),
                                    SL.library = SL.library,
                                    cvControl = list(V=SLfolds))
  
  if(SL_type =="discrete"){
    return(mod$library.predict[, which.min(mod$cvRisk)]) #discrete superlearner
  }else if(SL_type =="non_discrete"){
    return(mod$SL.predict) 
  }else stop("SL_type not recognised")
  
}

gamma_fit <- function(df,X,cov,SL.library,SLfolds=10,SL_type = "discrete"){
  ind <- df$inTrain
  X<- select(X,!matches(cov))
  
  if(all(ind)){
    newX <- X
  }else{
    newX <- X[!ind,]
  }
  
  mod <- SuperLearner::SuperLearner((df$CATE[ind])^2,X[ind,],newX=newX,
                                    family = gaussian(),
                                    SL.library = SL.library,
                                    cvControl = list(V=SLfolds))
  
  if(SL_type =="discrete"){
    return(mod$library.predict[, which.min(mod$cvRisk)]) #discrete superlearner
  }else if(SL_type =="non_discrete"){
    return(mod$SL.predict) 
  }else stop("SL_type not recognised")
  
}