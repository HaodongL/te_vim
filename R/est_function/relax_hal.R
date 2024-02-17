# library(here)
# library(sl3)
# source(paste0(here(), "/R/est_function/Lrnr_base.R"))
#' Scalable Highly Adaptive Lasso (HAL)
#'
#' The Highly Adaptive Lasso (HAL) is a nonparametric regression function that
#' has been demonstrated to optimally estimate functions with bounded (finite)
#' variation norm. The algorithm proceeds by first building an adaptive basis
#' (i.e., the HAL basis) based on indicator basis functions (or higher-order
#' spline basis functions) representing covariates and interactions of the
#' covariates up to a pre-specified degree. The fitting procedures included in
#' this learner use \code{\link[hal9001]{fit_hal}} from the \pkg{hal9001}
#' package. For details on HAL regression, consider consulting the following
#' \insertCite{benkeser2016hal;textual}{sl3}),
#' \insertCite{coyle2020hal9001-rpkg;textual}{sl3}),
#' \insertCite{hejazi2020hal9001-joss;textual}{sl3}).
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom origami folds2foldvec
#' @importFrom stats predict quasibinomial
#'
#' @export
#'
#' @keywords data
#'
#' @return A learner object inheriting from \code{\link{Lrnr_base}} with
#'  methods for training and prediction. For a full list of learner
#'  functionality, see the complete documentation of \code{\link{Lrnr_base}}.
#'
#' @format An \code{\link[R6]{R6Class}} object inheriting from
#'  \code{\link{Lrnr_base}}.
#'
#' @family Learners
#'
#' @section Parameters:
#'   - \code{...}: Arguments passed to \code{\link[hal9001]{fit_hal}}. See
#'    it's documentation for details.
#'
#' @examples
#' data(cpp_imputed)
#' covs <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs")
#' task <- sl3_Task$new(cpp_imputed, covariates = covs, outcome = "haz")
#'
#' # instantiate with max 2-way interactions, 0-order splines, and binning
#' # (i.e., num_knots) that decreases with increasing interaction degree
#' hal_lrnr <- Lrnr_hal9001$new(
#'   max_degree = 2, num_knots = c(20, 10), smoothness_orders = 0
#' )
#' hal_fit <- hal_lrnr$train(task)
#' hal_preds <- hal_fit$predict()
Lrnr_rhal9001 <- R6Class(
  classname = "Lrnr_rhal9001",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "weights", "ids"),
    .train = function(task) {
      args <- self$params
      
      args$X <- as.matrix(task$X)
      
      outcome_type <- self$get_outcome_type(task)
      args$Y <- outcome_type$format(task$Y)
      
      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family()
      }
      
      if (!any(grepl("fit_control", names(args)))) {
        args$fit_control <- list()
      }
      args$fit_control$foldid <- origami::folds2foldvec(task$folds)
      
      if (task$has_node("id")) {
        args$id <- task$id
      }
      
      if (task$has_node("weights")) {
        args$fit_control$weights <- task$weights
      }
      
      if (task$has_node("offset")) {
        args$offset <- task$offset
      }
      
      # fit HAL, allowing glmnet-fitting arguments
      other_valid <- c(
        names(formals(glmnet::cv.glmnet)), names(formals(glmnet::glmnet))
      )
      
      fit_object <- call_with_args(
        fit_rhal, args,
        other_valid = other_valid
      )
      
      return(fit_object)
    },
    .predict = function(task = NULL) {
      predictions <- stats::predict(
        self$fit_object,
        new_data = data.matrix(task$X)
      )
      if (!is.na(safe_dim(predictions)[2])) {
        p <- ncol(predictions)
        colnames(predictions) <- sprintf("lambda_%0.3e", self$params$lambda)
      }
      return(predictions)
    },
    .required_packages = c("hal9001", "glmnet")
  )
)


#' relax HAL
#'
#' Note: Current relax HAL use a global criterion.
#' Future work need to be done for user-specified criterion driven by the target parameter.
#'
#'

###############################################################################
#'  fit relax HAL
#'
#' @details a wrapper function that execuate relax HAL and return the fit. See
#' more parameter definitions in \code{fit_hal}.
#' @param Nlam A \code{integer} scalar of number of lambda candidates

fit_rhal <- function(X,
                     Y,
                     formula = NULL,
                     X_unpenalized = NULL,
                     max_degree = ifelse(ncol(X) >= 20, 2, 3),
                     smoothness_orders = 0,
                     num_knots = num_knots_generator(
                       max_degree = max_degree,
                       smoothness_orders = smoothness_orders,
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
                     yolo = FALSE){
  
  # initialize the undersmoothing procedure
  args_f <- as.list(environment())
  hal_init <- do.call(relax_init,
                      args = args_f)
  # do relaxed HAL
  # return relaxed HAL fit
  x_basis <- hal_init$fit_init$x_basis
  nonzero_col <- hal_init$nonzero_col
  
  # nonzero_basis <- hal_init$basis_mat
  # if (length(family) > 1) {family <- family[1]}
  # 
  # if (family == "gaussian"){
  #   family = gaussian()
  # }else if (family == "binomial"){
  #   family = binomial()
  # }else if (family == "poisson"){
  #   family = poisson()
  # }
  # rhal_fit <- glm.fit(x = nonzero_basis, y = Y, family = family)
  
  pf <- rep(Inf, ncol(x_basis))
  pf[nonzero_col] <- 0
  
  rhal_fit <- glmnet(x_basis, 
                     Y, 
                     lambda = 1e-5,
                     family = family, 
                     standardize = FALSE, 
                     penalty.factor  = pf)
  
  coefs <- stats::coef(rhal_fit)
  print(paste0("L1 norm: ", sum(abs(coefs[-1]))))
  # print(paste0("zero: ", coefs[-nonzero_col]))
  # print(paste0("non-zero: ", coefs[nonzero_col]))
  unpenalized_covariates <- hal_init$fit_init$unpenalized_covariates
  X_colnames <- hal_init$fit_init$X_colnames
  copy_map <- hal_init$fit_init$copy_map
  times <- hal_init$fit_init$times
  basis_list <- hal_init$fit_init$basis_list
  
  # construct output object via lazy S3 list
  fit <- list(
    x_basis =
      if (return_x_basis) {
        x_basis
      } else {
        NULL
      },
    basis_list = basis_list,
    X_colnames = X_colnames,
    copy_map = copy_map,
    coefs = as.matrix(coefs),
    times = times,
    reduce_basis = reduce_basis,
    family = family,
    lasso_fit =
      if (return_lasso) {
        rhal_fit
      } else {
        NULL
      },
    unpenalized_covariates = unpenalized_covariates,
    prediction_bounds = fit_control$prediction_bounds
  )
  class(fit) <- "hal9001"
  return(fit)
}

###############################################################################
#' initialize relax HAL
#'
#' @details fit a regular 0-order HAL(\code{smoothness_orders = 0}) with binning
#'  (\code{num_knots}), select the set of basis functions with non-zero coefs
#'  which will be used in following procedure. Note this is not the only valid initializing
#'  procedure, the user can define their own initialization.
#' @keywords internal

relax_init <- function(X, Y, ...) {
  
  # get initial fit
  args_f <- c(as.list(environment()), list(...))
  fit_init <- do.call(fit_hal, args = args_f)
  
  # select the non-zero directions/basis
  init_coef <-fit_init$coefs[-1]
  nonzero_col <- which(init_coef != 0)
  basis_mat <- as.matrix(fit_init$x_basis)
  basis_mat <- as.matrix(basis_mat[, nonzero_col])
  
  res <- list("fit_init" = fit_init,
              "basis_mat" = basis_mat,
              "nonzero_col" = nonzero_col)
  return(res)
}

###############################################################################

#' A default generator for the \code{num_knots} argument for each degree of
#' interactions and the smoothness orders.
#'
#' @param d interaction degree.
#' @param smoothness_orders see \code{\link{fit_hal}}.
#' @param base_num_knots_0 The base number of knots for zeroth-order smoothness
#'  basis functions. The number of knots by degree interaction decays as
#'  `base_num_knots_0/2^(d-1)` where `d` is the interaction degree of the basis
#'  function.
#' @param base_num_knots_1 The base number of knots for 1 or greater order
#'  smoothness basis functions. The number of knots by degree interaction
#'  decays as `base_num_knots_1/2^(d-1)` where `d` is the interaction degree of
#'  the basis function.
#'
#' @keywords internal

num_knots_generator <- function(max_degree, smoothness_orders, base_num_knots_0 = 500,
                                base_num_knots_1 = 200) {
  if (all(smoothness_orders > 0)) {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_1 / 2^(d - 1))
    }))
  } else {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_0 / 2^(d - 1))
    }))
  }
}


################################################################################

#' Call with filtered argument list
#'
#' Call a function with a list of arguments, eliminating any that aren't
#' matched in the function prototype
#'
#' @param fun A \code{function} whose signature will be used to reduce the
#' @param args A \code{list} of function arguments to use.
#' @param other_valid A \code{list} of function arguments names that are valid,
#'   but not formals of \code{fun}.
#' @param keep_all A \code{logical} don't drop arguments, even if they aren't
#'   matched in either the function prototype or other_valid.
#' @param silent A \code{logical} indicating whether to pass \code{message}s
#'  when arguments not found in \code{formals} are passed to \code{fun}.
#' @param ignore A \code{character} vector indicating which arguments should be dropped
#'
#' @keywords internal
call_with_args <- function(fun, args, other_valid = list(), keep_all = FALSE,
                           silent = FALSE, ignore = c()) {
  
  # drop ignore args
  args <- args[!(names(args) %in% ignore)]
  if (!keep_all) {
    # catch arguments to be kept
    formal_args <- names(formals(fun))
    all_valid <- c(formal_args, other_valid)
    
    # find invalid arguments based on combination of formals and other_valid
    invalid <- names(args)[which(!(names(args) %in% all_valid))]
    
    # subset arguments to pass
    args <- args[which(names(args) %in% all_valid)]
    
    # return warnings when dropping arguments
    if (!silent & length(invalid) > 0) {
      message(sprintf(
        "Learner called function %s with unknown args: %s. These will be dropped.\nCheck the params supported by this learner.",
        as.character(substitute(fun)), paste(invalid, collapse = ", ")
      ))
    }
  }
  do.call(fun, args)
}
