TEVIM_noCV_SL <- function(Y,A,X,SL.library,cov_list,
                          SLfolds=10,cluster=NULL,SL_type = "both",...){
  
  iFit <- fitMods_SL(Y,A,X,SL.library,fold=NULL,SLfolds=SLfolds,...)
  
  if(SL_type=="both"){
    SL_types <- c("discrete","non_discrete")
  }else{
    SL_types <- c(SL_type)
  }
  
  a <- expand_grid(x2=c("T","DR"),x3=SL_types) %>%
    as.list() %>%
    purrr::pmap(function(x2,x3)list(Learner=x2,SLt=x3))
  
  cate_fits <- parallel::parLapply(cl=cluster,a,function(z){
    if(z$Learner=="T"){
      z$cFit <- TLearn_SL(iFit,SL_type = z$SLt)
    }else if(z$Learner=="DR"){
      z$cFit <- DRLearn_SL(iFit,X,SL.library,SL_type = z$SLt)
    }
    (z)
  })
  
  b <- expand_grid(x1=cate_fits,x2=names(cov_list)) %>%
    as.list() %>%
    purrr::pmap(function(x1,x2){c(x1,list(cov=x2))})
  
  margcate_fits <- parallel::parLapply(cl=cluster,b,function(z){
    (list( mfit = marginal_CATEs(z$cFit,X,z$cov,SL.library,SLfolds,SL_type = z$SLt),
           Learner=z$Learner,SLt=z$SLt,cov=z$cov))
  })
  
  gamma_fits <- parallel::parLapply(cl=cluster,b,function(z){
    (list( mfit = gamma_fit(z$cFit,X,z$cov,SL.library,SLfolds,SL_type = z$SLt),
           Learner=z$Learner,SLt=z$SLt,cov=z$cov))
  })
  
  a <- expand_grid(x1=c("T","DR"),x2=SL_types) %>%
    as.list() %>%
    purrr::pmap(function(x1,x2)list(Learner=x1,SL_type=x2))
  
  out <- lapply(a,function(z){
    c_out   <- sapply(cate_fits,function(y){
      if(y$Learner==z$Learner  & y$SLt==z$SL_type) y$cFit
    }) %>% bind_rows()
    
    
    b <- sapply(margcate_fits,function(y){
      if(y$Learner==z$Learner  & y$SLt==z$SL_type) y$mfit
    }) 
    b <- do.call(cbind,b)
    
    colnames(b) <- names(cov_list)
    m_out <- as_tibble(b)
    
    
    d <- sapply(gamma_fits,function(y){
      if(y$Learner==z$Learner  & y$SLt==z$SL_type) y$mfit
    })
    d <- do.call(cbind,d)
    colnames(d) <- names(cov_list)
    gamma_out <- as_tibble(d)
    
    
    (list(cate_fit= c_out,
          maginal_fit = m_out,
          gamma_fit = gamma_out)) # 
  })
  
  names(out) <- lapply(a,function(z)paste0(z$Learner,"_",z$SL_type))
  # print(paste0("length out: ", length(out)))
  return(c(list(Initial_fit=iFit),out))
}



