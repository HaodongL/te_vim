



TEVIMcrossFit_SL <- function(Y,A,X,SL.library,foldIDs,cov_list,
                        SLfolds=10,cluster=NULL,SL_type = "both",...){
  Nfolds <- max(foldIDs)
  Ncov <- length(cov_list)
  
  initial_fits <- parallel::parSapply(cl=cluster,1:Nfolds,function(i){
    fold <- which(foldIDs==i)
    (fitMods_SL(Y,A,X,SL.library,fold=fold,SLfolds=SLfolds,...))
    },simplify = FALSE)
  
  if(SL_type=="both"){
    SL_types <- c("discrete","non_discrete")
  }else{
    SL_types <- c(SL_type)
  }

  #Get arguments for initial cate fits
  a <- expand_grid(x1=1:Nfolds,x2=c("T","DR"),x3=SL_types) %>%
    as.list() %>%
    purrr::pmap(function(x1,x2,x3)list(i=x1,Learner=x2,SLt=x3))

  cate_fits <- parallel::parLapply(cl=cluster,a,function(z){
      fit <- initial_fits[[z$i]]
      if(z$Learner=="T"){
        z$cFit <- TLearn_SL(fit,SL_type = z$SLt)
      }else if(z$Learner=="DR"){
        z$cFit <- DRLearn_SL(fit,X,SL.library,SL_type = z$SLt)
      }
      (z)
    })
  
  #Get arguments for submodel CATE fits
  b <- expand_grid(x1=cate_fits,x2=names(cov_list)) %>%
    as.list() %>%
    purrr::pmap(function(x1,x2){c(x1,list(cov=x2))})
  
  margcate_fits <- parallel::parLapply(cl=cluster,b,function(z){
    (list( mfit = marginal_CATEs(z$cFit,X,z$cov,SL.library,SLfolds,SL_type = z$SLt),
          i=z$i, Learner=z$Learner,SLt=z$SLt,cov=z$cov))
  })

  
  ##Now to aggregate results accross folds
  ord <- lapply(1:Nfolds, function(i){which(foldIDs==i)}) %>%
    unlist() %>% order()
  
  iFit <- lapply(initial_fits,function(fit){
    (fit %>% filter(inTrain==FALSE) %>% select(-inTrain) )
    }) %>% bind_rows()
  iFit <- iFit[ord,] #initial pi,mu1,mu0 fits with correct ordering

  cate_fits <- lapply(cate_fits,function(fit){
    fit$cFit <- filter(fit$cFit,inTrain==FALSE) %>% select(-inTrain)
    (fit)
  })
  
  a <- expand_grid(x1=c("T","DR"),x2=SL_types) %>%
    as.list() %>%
    purrr::pmap(function(x1,x2)list(Learner=x1,SL_type=x2))
  
  out <- lapply(a,function(z){
    
    c_out   <- sapply(cate_fits,function(y){
      if(y$Learner==z$Learner  & y$SLt==z$SL_type) y$cFit
      }) %>% bind_rows()
    c_out <- c_out[ord,]

    m_out <- sapply(1:Nfolds, function(i){
      b <- sapply(margcate_fits,function(y){
        if(y$Learner==z$Learner  & y$SLt==z$SL_type & y$i==i) y$mfit
      }) 
      b <- do.call(cbind,b)
      colnames(b) <- names(cov_list)
      as_tibble(b)
    },simplify=FALSE) %>% bind_rows()
    m_out <- m_out[ord,]
    (list(cate_fit=mutate(c_out,FoldID=foldIDs),
          maginal_fit = m_out))
  })
  
  names(out) <- lapply(a,function(z)paste0(z$Learner,"_",z$SL_type))
  
  return(c(list(Initial_fit=iFit),out))
}
