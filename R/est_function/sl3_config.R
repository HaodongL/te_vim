library(sl3)

# -- sl setup
lrnr_lm <- Lrnr_glm$new()
lrnr_mean <- Lrnr_mean$new()
lrnr_earth <- Lrnr_earth$new()
lrnr_lasso <- Lrnr_glmnet$new()
lrnr_xgb <- Lrnr_xgboost$new()
lrnr_ranger <- Lrnr_ranger$new()
lrnr_gam <- Lrnr_gam$new()
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

sl_Q <- Lrnr_sl$new(
  learners = lrnr_stack_Q,
  metalearner = lb_metalearner
)

# sl_Q <- Lrnr_sl$new(
#   learners = lrnr_stack_Q,
#   metalearner = ls_metalearner
# )

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

# sl_Q <- lrnr_hal
sl_g <- lrnr_earth
sl_x <- lrnr_earth


