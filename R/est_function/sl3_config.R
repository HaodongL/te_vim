library(sl3)

# -- sl setup
lrnr_lm <- Lrnr_glm$new()
lrnr_mean <- Lrnr_mean$new()
lrnr_earth <- Lrnr_earth$new()
lrnr_lasso <- Lrnr_glmnet$new()
lrnr_xgb <- Lrnr_xgboost$new()
lrnr_ranger <- Lrnr_ranger$new()
lrnr_gam_Q <- Lrnr_gam$new('Y ~ s(X1) + s(X2) + ti(X1,X2) + s(X1,by=A) + s(X2,by=A) + ti(X1,X2,by=A)')
lrnr_gam_g <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# lrnr_gam_tau <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# lrnr_gam_tau_s <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# lrnr_gam_gamma_s <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# lrnr_gam <- Lrnr_gam$new()
# lrnr_polspline <- Lrnr_polspline$new()
lrnr_hal <- Lrnr_hal9001$new(num_knots = c(50, 25, 15))
# lrnr_grf <- Lrnr_grf$new()
# lrnr_lm_inter<- Lrnr_glm$new(formula = "~.^2")

# lrnr_stack <- make_learner("Stack",
#                            lrnr_lm,
#                            lrnr_earth)

lrnr_stack_Q <- make_learner("Stack",
                           lrnr_lm,
                           lrnr_lasso,
                           lrnr_xgb,
                           lrnr_earth)

lrnr_stack_g <- make_learner("Stack",
                             lrnr_lm,
                             lrnr_lasso,
                             lrnr_xgb,
                             lrnr_earth)

lrnr_stack_x <- make_learner("Stack",
                             lrnr_lm,
                             lrnr_lasso,
                             lrnr_xgb,
                             lrnr_earth)


# ls_metalearner <- make_learner(Lrnr_nnls)
# lb_metalearner <- make_learner(Lrnr_solnp,
#                                learner_function = metalearner_logistic_binomial,
#                                loss_function = loss_loglik_binomial)
ls_metalearner <- Lrnr_cv_selector$new(loss_squared_error)
lb_metalearner <- Lrnr_cv_selector$new(loss_loglik_binomial)

# sl_Q <- Lrnr_sl$new(
#   learners = lrnr_stack_Q,
#   metalearner = Lrnr_cv_selector$new(loss_loglik_binomial)
# )
# 
# sl_g <- Lrnr_sl$new(
#   learners = lrnr_stack_g,
#   metalearner = lb_metalearner,
#   outcome_type = 'binomial'
# )
# 
# sl_x <- Lrnr_sl$new(
#   learners = lrnr_stack_x,
#   metalearner = ls_metalearner
# )

sl_Q <- lrnr_earth
sl_g <- lrnr_earth
sl_x <- lrnr_earth
# sl_x_t <- lrnr_lm

# # -- sl modeling function
# fitSL <- function(df, Q_bounds = NULL, g_bounds = c(0.025, 0.975), cv = TRUE){
#   
#   folds <- origami::make_folds(strata_ids = df$A)
#   # folds <- origami::folds_vfold(nrow(df))
#   
#   # setup sl3
#   task_Q <- sl3::make_sl3_Task(
#     data = df,
#     covariates = setdiff(names(df), c('Y')),
#     outcome = 'Y',
#     folds = folds
#   )
#   
#   task_g <- sl3::make_sl3_Task(
#     data = df,
#     covariates = setdiff(names(df), c('Y', 'A')),
#     outcome = 'A',
#     folds = folds
#   )
#   
#   sl_Q <- Lrnr_sl$new(
#     learners = lrnr_stack_Q,
#     metalearner = ls_metalearner
#   )
#   
#   sl_g <- Lrnr_sl$new(
#     learners = lrnr_stack_g,
#     metalearner = lb_metalearner,
#     outcome_type = 'binomial'
#   )
#   
#   # fit Q and g
#   Q_fit <- sl_Q$train(task_Q)
#   g_fit <- sl_g$train(task_g)
#   
#   # preds Q and g
#   # task_Q_pred <- sl3::make_sl3_Task(
#   #   data = df,
#   #   covariates = setdiff(names(df), c('Y')),
#   #   outcome = 'Y',
#   #   folds = folds
#   # )
#   # pred_Q <- Q_fit$predict_fold(task_Q, "validation")
#   
#   # pred_Q <- Q_fit$predict()
#   pred_Q_cf <- pred_Q_cf(df = df, Q_fit = Q_fit, folds = folds, cv = cv)
#   
#   # QbarAW <- pred_Q
#   Qbar1W <- pred_Q_cf$Qbar1W
#   Qbar0W <- pred_Q_cf$Qbar0W
#   QbarAW <- ifelse(df$A == 1, Qbar1W, Qbar0W)
#   
#   # bound Q (is there a better way???)
#   if (!is.null(Q_bounds)){
#     QbarAW <- bound(QbarAW, Q_bounds)
#     Qbar1W <- bound(Qbar1W, Q_bounds)
#     Qbar0W <- bound(Qbar0W, Q_bounds)
#   }
#   
#   # task_g_pred <- sl3::make_sl3_Task(
#   #   data = df,
#   #   covariates = setdiff(names(df), c('Y', 'A')),
#   #   outcome = 'A',
#   #   folds = folds
#   # )
#   
#   if (cv){
#     pred_g <- g_fit$predict_fold(task_g, "validation")
#   }else{
#     pred_g <- g_fit$predict()
#   }
#   
#   # bound g
#   pred_g <- bound(pred_g, g_bounds)
#   
#   # mod.m <- gam(Y~s(X1) + s(X2) + ti(X1,X2)+s(X1,by=A) + s(X2,by=A) + ti(X1,X2,by=A),
#   #              family = gaussian(),data=df_train)
#   # mod.ps <- gam(A~s(X1)+s(X2)+ti(X1,X2),
#   #               family = binomial(),data=df_train) 
#   
#   # mod.m <- ranger::ranger(Y ~ .,data = df_train, num.trees =  200)
#   
#   tibble(
#     Y = df$Y, 
#     A = df$A,
#     pi_hat  = pred_g,
#     mu1_hat = Qbar1W,
#     mu0_hat = Qbar0W,
#     mua_hat = QbarAW
#     # mu1_hat = predict(mod.m,mutate(df,A=1),type="response"),
#     # mu0_hat = predict(mod.m,mutate(df,A=0),type="response")
#   ) %>% return()
# }


# # helper function to predict Qbar1W and Qbar0W
# pred_Q_cf <- function(df, Q_fit, folds, cv){
#   
#   # cf data
#   df1 <- df
#   df1$A <- 1
#   df0 <- df
#   df0$A <- 0
#   
#   # cf tasks
#   Q_task1 <- sl3::make_sl3_Task(
#     data = df1,
#     covariates = setdiff(names(df1), c('Y')),
#     outcome = 'Y',
#     folds = folds
#   )
#   
#   Q_task0 <- sl3::make_sl3_Task(
#     data = df0,
#     covariates = setdiff(names(df0), c('Y')),
#     outcome = 'Y',
#     folds = folds
#   )
#   
#   # cf Q
#   if (cv){
#     Qbar1W <- Q_fit$predict_fold(Q_task1, "validation")
#     Qbar0W <- Q_fit$predict_fold(Q_task0, "validation")
#   }else{
#     Qbar1W <- Q_fit$predict(Q_task1)
#     Qbar0W <- Q_fit$predict(Q_task0)
#   }
#   
#   return(list('Qbar1W' = Qbar1W,
#               'Qbar0W' = Qbar0W))
# }
# 
# # crossFitSL <- function(df,foldIDs,Nfolds=max(foldIDs)){
# #   fold_list = lapply(1:Nfolds, function(i,foldIDs){which(foldIDs==i)},foldIDs = foldIDs)
# #   
# #   a <- sapply(fold_list,function(fold,df){
# #     fitSL(df[fold,],df[-1*fold,]) 
# #   },df=df,simplify = FALSE) %>% bind_rows()
# #   
# #   a <- a[order(unlist(fold_list)),] %>%
# #     mutate(FoldID = foldIDs)
# #   
# #   return(a)
# # }
# 
# # getFoldIDs <- function(N,Nfolds,shuffle=TRUE){
# #   if (shuffle){
# #     return(sample(rep_len(seq_len(Nfolds),length.out=N)))
# #   } else{
# #     return(rep_len(seq_len(Nfolds),length.out=N))
# #   }
# # }
# 
# bound <- function(x, bounds) {
#   lower <- bounds[[1]]
#   if (length(bounds) > 1) {
#     upper <- bounds[[2]]
#   } else {
#     upper <- 1 - lower
#   }
#   pmin(pmax(x, lower), upper)
# }
