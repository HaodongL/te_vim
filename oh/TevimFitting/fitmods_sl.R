#require(SuperLearner)
require(tidyverse)

fitMods_SL <- function(Y,A,X,SL.library,PS.method="SuperLearn",
                       SLfolds=10,Y_family=gaussian(),
                       fold=NULL){
  
  if(is.numeric(fold)){
    Y_train <- Y[-fold]
    A_train <- A[-fold]
    X_train <- X[-fold,]
  }else{
    Y_train <- Y
    A_train <- A
    X_train <- X
  }
  
  if(PS.method=="SuperLearn"){
    pi_fit <- SuperLearner::SuperLearner(A_train,X_train,newX=X,
                                         family = binomial(),
                                         SL.library = SL.library,
                                         cvControl = list(V=SLfolds))
    pi_hat_disc <- pi_fit$library.predict[, which.min(pi_fit$cvRisk)] #discrete superlearner
    pi_hat_SL   <- pi_fit$SL.predict #non-discrete SL
    
  }else if(PS.method=="SampleMean"){
    pi_hat_disc  <- pi_hat_SL  <- mean(A)
  }else if(PS.method=="OutofSampleMean"){
    pi_hat_disc  <- pi_hat_SL  <- mean(A_train)
  }
  
  mu_fit <- SuperLearner::SuperLearner(Y_train,cbind(A=A_train,X_train),
                                       newX = rbind(cbind(A=1,X),cbind(A=0,X)),
                                       family = Y_family,
                                       SL.library = SL.library,
                                       cvControl = list(V=SLfolds)) 
  
  mu_hat_disc <- mu_fit$library.predict[, which.min(mu_fit$cvRisk)] #discrete superlearner
  mu_hat_SL   <- mu_fit$SL.predict #non-discrete SL
  
  N <- NROW(Y)
  
  df_discrete <- tibble(
    Y = Y, A=A,
    pi_hat  = pi_hat_disc,
    mu1_hat = mu_hat_disc[1:N],
    mu0_hat = mu_hat_disc[(N+1):(2*N)],
    pi_hat_SL  = pi_hat_SL,
    mu1_hat_SL = mu_hat_SL[1:N],
    mu0_hat_SL = mu_hat_SL[(N+1):(2*N)],
    inTrain    = TRUE
  )
  df_discrete[fold,"inTrain"]<- FALSE
  
  return(df_discrete)
}







