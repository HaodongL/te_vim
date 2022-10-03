require(tidyverse)
require(mgcv) #we use gam fitting for this example

generate_data_simple <- function(N, print_truth = FALSE){
  # A simple data generating process
  X1 <- runif(N,-1,1)
  X2 <- runif(N,-1,1)
  pi0 <- plogis(0.1*X1*X2-0.4*X1)
  A <- rbinom(N,1,prob=pi0)
  muY0 <- X1*X2 + 2*X2^2 -X1
  CATE <- X1^2*(X1+7/5) + (5*X2/3)^2
  muY = muY0+A*CATE
  Y <- rnorm(N,sd=1,mean= muY)
  
  if (print_truth){
    VTE <- var(CATE)
    tau_s <- X1^2*(X1+7/5)
    VTE_s <- var(tau_s)
    VIM_Theta_s <- VTE - VTE_s
    print(paste0("VIM_Theta_s: ", VIM_Theta_s))
    print(paste0("VTE: ", VTE))
  }
  return(tibble(X1=X1,X2=X2,A=A,Y=Y))
}


# fitMods <- function(df,df_train=df,outcome_family = gaussian()){
#   #Here we use a gam model as an example
#   #Generally one could use any fitting procedure of interest
#   mod.ps <- gam(A~s(X1)+s(X2)+ti(X1,X2),
#                 family = binomial(),data=df_train) 
#   mod.m <- gam(Y~s(X1) + s(X2) + ti(X1,X2)+s(X1,by=A) + s(X2,by=A) + ti(X1,X2,by=A),
#                family = outcome_family,data=df_train)
#   tibble(
#     Y = df$Y, A=df$A,
#     pi_hat  = predict(mod.ps,df,type="response"),
#     mu1_hat = predict(mod.m,mutate(df,A=1),type="response"),
#     mu0_hat = predict(mod.m,mutate(df,A=0),type="response")
#   )%>% return()
# }
# 
# 
# crossFit <- function(df,foldIDs,Nfolds=max(foldIDs)){
#   fold_list = lapply(1:Nfolds, function(i,foldIDs){which(foldIDs==i)},foldIDs = foldIDs)
#   
#   a <- sapply(fold_list,function(fold,df){
#     fitMods(df[fold,],df[-1*fold,]) 
#   },df=df,simplify = FALSE) %>% bind_rows()
# 
#   a <- a[order(unlist(fold_list)),] %>%
#     mutate(FoldID = foldIDs)
# 
#   return(a)
# }
# 
# getFoldIDs <- function(N,Nfolds,shuffle=TRUE){
#   if (shuffle){
#     return(sample(rep_len(seq_len(Nfolds),length.out=N)))
#   } else{
#     return(rep_len(seq_len(Nfolds),length.out=N))
#   }
# }