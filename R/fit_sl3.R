library(sl3)

# -- sl setup
lrnr_lm <- Lrnr_glm$new()
lrnr_lasso <- Lrnr_glmnet$new()
lrnr_xgb <- Lrnr_xgboost$new()
lrnr_gam <- Lrnr_gam$new()
lrnr_earth <- Lrnr_earth$new()
lrnr_ranger <- Lrnr_ranger$new()

lrnr_stack <- make_learner("Stack",
                           lrnr_lm,
                           lrnr_earth,
                           lrnr_xgb)


ls_metalearner <- make_learner(Lrnr_nnls)
lb_metalearner <- make_learner(Lrnr_solnp,
                               learner_function = metalearner_logistic_binomial,
                               loss_function = loss_loglik_binomial)

# -- sl modeling function
fitSL <- function(df, Q_bounds = NULL, g_bounds = c(0.025, 0.975)){
  
  folds <- origami::make_folds(strata_ids = df$A)
  
  # setup sl3
  task_Q <- sl3::make_sl3_Task(
    data = df,
    covariates = setdiff(names(df), c('Y')),
    outcome = 'Y',
    folds = folds
  )
  
  task_g <- sl3::make_sl3_Task(
    data = df,
    covariates = setdiff(names(df), c('Y', 'A')),
    outcome = 'A',
    folds = folds
  )
  
  sl_Q <- Lrnr_sl$new(
    learners = lrnr_stack,
    metalearner = ls_metalearner,
    outcome_type = 'continuous'
  )
  
  sl_g <- Lrnr_sl$new(
    learners = lrnr_stack,
    metalearner = lb_metalearner,
    outcome_type = 'binomial'
  )
  
  # fit Q and g
  Q_fit <- sl_Q$train(task_Q)
  g_fit <- sl_g$train(task_g)
  
  # preds Q and g
  task_Q_pred <- sl3::make_sl3_Task(
    data = df,
    covariates = setdiff(names(df), c('Y')),
    outcome = 'Y',
    folds = folds
  )
  # pred_Q <- Q_fit$predict(task_Q_pred)
  pred_Q <- Q_fit$predict_fold(task_Q_pred, "validation")
  pred_Q_cf <- pred_Q_cf(df = df, Q_fit = Q_fit, folds = folds)
  
  QbarAW <- pred_Q
  Qbar1W <- pred_Q_cf$Qbar1W
  Qbar0W <- pred_Q_cf$Qbar0W
  
  # bound Q (is there a better way???)
  if (!is.null(Q_bounds)){
    QbarAW <- bound(QbarAW, Q_bounds)
    Qbar1W <- bound(Qbar1W, Q_bounds)
    Qbar0W <- bound(Qbar0W, Q_bounds)
  }
  
  task_g_pred <- sl3::make_sl3_Task(
    data = df,
    covariates = setdiff(names(df), c('Y', 'A')),
    outcome = 'A',
    folds = folds
  )
  pred_g <- g_fit$predict_fold(task_g_pred, "validation")
  
  # bound g
  pred_g <- bound(pred_g, g_bounds)
  
  # mod.m <- gam(Y~s(X1) + s(X2) + ti(X1,X2)+s(X1,by=A) + s(X2,by=A) + ti(X1,X2,by=A),
  #              family = gaussian(),data=df_train)
  # mod.ps <- gam(A~s(X1)+s(X2)+ti(X1,X2),
  #               family = binomial(),data=df_train) 
  
  # mod.m <- ranger::ranger(Y ~ .,data = df_train, num.trees =  200)
  
  tibble(
    Y = df$Y, 
    A = df$A,
    pi_hat  = pred_g,
    mu1_hat = Qbar1W,
    mu0_hat = Qbar0W,
    mua_hat = QbarAW
    # mu1_hat = predict(mod.m,mutate(df,A=1),type="response"),
    # mu0_hat = predict(mod.m,mutate(df,A=0),type="response")
  ) %>% return()
}


# helper function to predict Qbar1W and Qbar0W
pred_Q_cf <- function(df, Q_fit, folds){
  
  # cf data
  df1 <- df
  df1$A <- 1
  df0 <- df
  df0$A <- 0
  
  # cf tasks
  Q_task1 <- sl3::make_sl3_Task(
    data = df1,
    covariates = setdiff(names(df1), c('Y')),
    outcome = 'Y',
    folds = folds
  )
  
  Q_task0 <- sl3::make_sl3_Task(
    data = df0,
    covariates = setdiff(names(df0), c('Y')),
    outcome = 'Y',
    folds = folds
  )
  
  # cf Q
  Qbar1W <- Q_fit$predict_fold(Q_task1,"validation")
  Qbar0W <- Q_fit$predict_fold(Q_task0,"validation")
  
  return(list('Qbar1W' = Qbar1W,
              'Qbar0W' = Qbar0W))
}

# crossFitSL <- function(df,foldIDs,Nfolds=max(foldIDs)){
#   fold_list = lapply(1:Nfolds, function(i,foldIDs){which(foldIDs==i)},foldIDs = foldIDs)
#   
#   a <- sapply(fold_list,function(fold,df){
#     fitSL(df[fold,],df[-1*fold,]) 
#   },df=df,simplify = FALSE) %>% bind_rows()
#   
#   a <- a[order(unlist(fold_list)),] %>%
#     mutate(FoldID = foldIDs)
#   
#   return(a)
# }

# getFoldIDs <- function(N,Nfolds,shuffle=TRUE){
#   if (shuffle){
#     return(sample(rep_len(seq_len(Nfolds),length.out=N)))
#   } else{
#     return(rep_len(seq_len(Nfolds),length.out=N))
#   }
# }

bound <- function(x, bounds) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  pmin(pmax(x, lower), upper)
}
