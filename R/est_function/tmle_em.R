
# cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
#               "vkantag", "caantag", "thiazide", "loopdiur")
# cm = cm_names[1]

tmle_em <- function(df, cm_names){
  df_res <- data.frame()
  for (cm in cm_names){
    df_train <- df %>% 
      rename("cm" = all_of(cm)) %>% 
      mutate(A = case_when(A==1 & cm == TRUE ~ "11",
                           A==0 & cm == TRUE ~ "01",
                           A==1 & cm == FALSE ~ "10",
                           A==0 & cm == FALSE ~ "00"))
    
    df_train$A <- as.factor(df_train$A)
    df_train <- df_train %>% select(-cm)
    
    node_list <- list(
      W = setdiff(names(df_train), c("Y", "A")),
      A = "A",
      Y = "Y"
    )
    tsm_spec <- tmle_TSM_all()
    tmle_task <- tsm_spec$make_tmle_task(df_train, node_list)
    
    learner_list <- list(A = sl_g_mn, Y = sl_Q)
    
    initial_likelihood <- tsm_spec$make_initial_likelihood(tmle_task, learner_list)
    targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
    
    all_tsm_params <- tsm_spec$make_params(tmle_task, targeted_likelihood)
    ate_param_true <- define_param(
      Param_delta, targeted_likelihood,
      delta_param_ATE,
      list(all_tsm_params[[2]], all_tsm_params[[4]])
    )
    ate_param_false <- define_param(
      Param_delta, targeted_likelihood,
      delta_param_ATE,
      list(all_tsm_params[[1]], all_tsm_params[[3]])
    )
    ate_param_diff <- define_param(
      Param_delta, targeted_likelihood,
      delta_param_ATE,
      list(ate_param_false, ate_param_true)
    )
    
    all_params <- c(all_tsm_params, 
                    ate_param_true, 
                    ate_param_false,
                    ate_param_diff)
    
    tmle_fit_multiparam <- fit_tmle3(
      tmle_task, targeted_likelihood, all_params,
      targeted_likelihood$updater
    )
    
    psi_hat <- tail(tmle_fit_multiparam$summary$tmle_est, 3)
    se <- tail(tmle_fit_multiparam$summary$se, 3)
    
    df_res_i <- data.frame("tau" = psi_hat,
                           "se" = se,
                           "group" = c(paste0("T_", str_sub(cm, 1,6)),
                                       paste0("F_", str_sub(cm, 1,6)),
                                       paste0("Diff_", str_sub(cm, 1,6))))
    df_res <- rbind(df_res, df_res_i)
  }
  return(df_res)
}


run_tmle3 <- function(df){
  node_list <- list(
    W = setdiff(names(df), c("Y", "A")),
    A = "A",
    Y = "Y"
  )
  
  ate_spec <- tmle_ATE(
    treatment_level = 1,
    control_level = 0
  )
  
  learner_list <- list(A = sl_g, Y = sl_Q)
  
  tmle_task <- ate_spec$make_tmle_task(df, node_list)
  initial_likelihood <- ate_spec$make_initial_likelihood(tmle_task, learner_list)
  # without CV-TMLE
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood,
                                                 updater = list(cvtmle = FALSE))
  ate_params <- ate_spec$make_params(tmle_task, targeted_likelihood) 
  tmle3_fit <- fit_tmle3(tmle_task, targeted_likelihood, 
                         ate_params, targeted_likelihood$updater)
  return(tmle3_fit)
}






















