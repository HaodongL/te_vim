
library(here)
library(R6)
library(hal9001)
library(glmnet)

rm(list = ls())
source(paste0(here(), "/R/est_function/relax_hal.R"))


# fit relax HAL
lrnr_rhal <- Lrnr_rhal9001$new()


# load example data set
data(cpp)
cpp <- cpp %>%
  dplyr::filter(!is.na(haz)) %>%
  mutate_all(~ replace(., is.na(.), 0))

# use covariates of intest and the outcome to build a task object
covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs",
            "sexn")
task <- sl3_Task$new(
  data = cpp,
  covariates = covars,
  outcome = "haz"
)

# stack learners into a model (including screeners and pipelines)
rhal_fit <- lrnr_rhal$train(task)
preds <- rhal_fit$predict()




X <- as.matrix(task$X)
outcome_type <- lrnr_rhal$get_outcome_type(task)
Y <- outcome_type$format(task$Y)

hal_init <-
relax_init(X = X,
           Y = Y,
           formula = NULL,
           X_unpenalized = NULL,
           max_degree = ifelse(ncol(X) >= 20, 2, 3),
           smoothness_orders = 0,
           num_knots = num_knots_generator(
             max_degree = 3,
             smoothness_orders = 0,
             base_num_knots_0 = 200,
             base_num_knots_1 = 50
           ),
           reduce_basis = 1 / sqrt(length(Y)),
           family = c("gaussian", "binomial", "poisson", "cox"),
           lambda = NULL,
           id = NULL,
           offset = NULL,
           fit_control = list(
             cv_select = TRUE,
             n_folds = 10,
             foldid = NULL,
             use_min = TRUE,
             lambda.min.ratio = 1e-4,
             prediction_bounds = "default"
           ),
           basis_list = NULL,
           return_lasso = TRUE,
           return_x_basis = TRUE,
           yolo = FALSE)

x_basis <- hal_init$fit_init$x_basis
nonzero_col <- hal_init$nonzero_col
pf <- rep(Inf, ncol(x_basis))
pf[nonzero_col] <- 0

rhal_fit1 <- glmnet(x_basis, 
                   Y, 
                   lambda = 1e-5,
                   family = "gaussian", 
                   standardize = FALSE, 
                   penalty.factor  = pf)
coefs1 <- stats::coef(rhal_fit1)


rhal_fit2 <- glmnet(x_basis, 
                    Y, 
                    lambda = 1e5,
                    family = "gaussian", 
                    standardize = FALSE, 
                    penalty.factor  = pf)
coefs2 <- stats::coef(rhal_fit2)


print(paste0("zero: ", coefs1[-1][-nonzero_col]))
print(paste0("non-zero: ", coefs1[-1][nonzero_col]))


print(paste0("zero: ", coefs2[-1][-nonzero_col]))
print(paste0("non-zero: ", coefs2[-1][nonzero_col]))

