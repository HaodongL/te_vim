not_constant <- function(x){
  res = !all(duplicated(x)[-1L])
  return(res)
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
                                                 updater = list(cvtmle = TRUE))
  ate_params <- ate_spec$make_params(tmle_task, targeted_likelihood) 
  tmle3_fit <- fit_tmle3(tmle_task, targeted_likelihood, 
                         ate_params, targeted_likelihood$updater)
  return(tmle3_fit)
}


tmle_stratified <- function(df, cm_names){
  n_cm <- length(cm_names)
  df_res <- foreach(i = 1:n_cm, .combine = 'rbind') %dopar% {
    # print(paste0("CM name: ", cm_names[i]))
    cm <- cm_names[i]
    df_true <- df %>% 
      filter(if_any(all_of(cm), ~ .x == "TRUE")) %>% 
      select(-all_of(cm)) %>% 
      select(where(not_constant))
    
    df_false <- df %>% 
      filter(if_any(all_of(cm), ~ .x == "FALSE")) %>% 
      select(-all_of(cm)) %>% 
      select(where(not_constant))
    
    tmle_fit1 <- run_tmle3(df_true)
    psi_hat1 <- tmle_fit1$summary$tmle_est
    se1 <- tmle_fit1$summary$se
    
    tmle_fit0 <- run_tmle3(df_false)
    psi_hat0 <- tmle_fit0$summary$tmle_est
    se0 <- tmle_fit0$summary$se
    
    df_res_i <- data.frame("tau" = c(psi_hat1, psi_hat0),
                           "se" = c(se1, se0),
                           "group" = c(paste0("T_", stringr::str_sub(cm, 1,6)),
                                       paste0("F_", stringr::str_sub(cm, 1,6))))
  }
  return(df_res)
}
  