

# -- sl setup
lrnr_lm <- Lrnr_glm$new()
lrnr_mean <- Lrnr_mean$new()
lrnr_earth <- Lrnr_earth$new()
lrnr_lasso <- Lrnr_glmnet$new()
lrnr_xgb <- Lrnr_xgboost$new()
lrnr_xgb100 <- Lrnr_xgboost$new(nrounds=100)
lrnr_xgb2 <- Lrnr_xgboost$new(nrounds = 2000, max_depth = 2)
lrnr_xgb3 <- Lrnr_xgboost$new(nrounds = 2000, max_depth = 3)
lrnr_ranger <- Lrnr_ranger$new()
lrnr_ranger2 <- Lrnr_ranger$new(num.trees = 2000, mtry = 2)
lrnr_ranger3 <- Lrnr_ranger$new(num.trees = 2000, mtry = 3)
lrnr_gam <- Lrnr_gam$new()
lrnr_gam_Q <- Lrnr_gam$new('Y ~ s(X1) + s(X2) + ti(X1,X2) + s(X1,by=A) + s(X2,by=A) + ti(X1,X2,by=A)')
lrnr_gam_g <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# lrnr_gam_tau <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# lrnr_gam_tau_s <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# lrnr_gam_gamma_s <- Lrnr_gam$new('A ~ s(X1) + s(X2) + ti(X1,X2)')
# lrnr_gam <- Lrnr_gam$new()
# lrnr_polspline <- Lrnr_polspline$new()
lrnr_hal <- Lrnr_hal9001$new()
lrnr_hal_fast <- Lrnr_hal9001$new(num_knots = c(40, 15, 10))
lrnr_hal_xfast <- Lrnr_hal9001$new(num_knots = c(20, 10, 5))
# lrnr_uhal_fast <- Lrnr_uhal9001$new(num_knots = c(50, 25, 15))
# lrnr_grf <- Lrnr_grf$new()
# lrnr_lm_inter<- Lrnr_glm$new(formula = "~.^2")
grid_params <- list(
  max_depth = c(3, 5, 8),
  eta = c(0.05, 0.1, 0.3),
  nrounds = 100
)

grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)

xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_xgboost$new, as.list(tuning_params))
})

# lrnr_stack <- make_learner("Stack",
#                            lrnr_lm,
#                            lrnr_earth)
lrnr_stack_g <- make_learner("Stack",
                             lrnr_lm,
                             lrnr_xgb,
                             lrnr_earth)

# lrnr_stack_Q <- make_learner("Stack",
#                              lrnr_xgb2,
#                              lrnr_xgb3,
#                              lrnr_ranger2,
#                              lrnr_ranger3,
#                              lrnr_lm,
#                              lrnr_hal_fast)
# 
# lrnr_stack_x <- make_learner("Stack",
#                              lrnr_xgb2,
#                              lrnr_xgb3,
#                              lrnr_ranger2,
#                              lrnr_ranger3,
#                              lrnr_lm)

lrnr_stack_Q <- make_learner("Stack",
                             lrnr_xgb,
                             lrnr_ranger,
                             lrnr_lm,
                             lrnr_earth,
                             lrnr_hal_fast)

lrnr_stack_x <- make_learner("Stack",
                             lrnr_xgb,
                             lrnr_ranger,
                             lrnr_lm,
                             lrnr_earth,
                             lrnr_hal_fast)

ls_metalearner <- make_learner(Lrnr_nnls)
lb_metalearner <- make_learner(Lrnr_solnp,
                               learner_function = metalearner_logistic_binomial,
                               loss_function = loss_loglik_binomial)
mn_metalearner <- make_learner(Lrnr_solnp, 
                               metalearner_linear_multinomial,
                               loss_loglik_multinomial)

ls_metalearner <- Lrnr_cv_selector$new(loss_squared_error)
lb_metalearner <- Lrnr_cv_selector$new(loss_loglik_binomial)

sl_Q <- Lrnr_sl$new(
  learners = lrnr_stack_Q,
  metalearner = lb_metalearner,
  outcome_type = 'binomial'
)

sl_Q <- Lrnr_sl$new(
  learners = lrnr_stack_Q,
  metalearner = ls_metalearner
)

sl_g <- Lrnr_sl$new(
  learners = lrnr_stack_g,
  metalearner = lb_metalearner,
  outcome_type = 'binomial'
)

sl_g_mn <- Lrnr_sl$new(
  learners = list(xgb_learners[[1]], 
                  xgb_learners[[3]],
                  xgb_learners[[5]],
                  lrnr_earth),
  metalearner = mn_metalearner
)

sl_x <- Lrnr_sl$new(
  learners = lrnr_stack_x,
  metalearner = ls_metalearner
)

# sl_Q <- lrnr_earth
# sl_g <- lrnr_earth
# sl_x <- lrnr_earth

sl_Q_hal <- lrnr_hal_fast
sl_g_hal <- lrnr_hal_fast
sl_ws_hal <- lrnr_hal_fast
