run_all_simu <- function(B, N, truth, covar = c('X2'), max.it = 1e3){
  results_cols <- c('i', 'truth', 'cvtmle', 'cvtmle_se',
                    'cvtmle_lower', 'cvtmle_upper', 
                    'cvaiptw', 'cvaiptw_se', 'cvaiptw_lower', 
                    'cvaiptw_upper')
  
  results_df <- data.frame(matrix(NA, nrow = B, ncol = length(results_cols)))
  colnames(results_df) <- results_cols
  
  run_bootstrap <- foreach(b = 1:B, .combine = 'rbind') %dopar% {
    
    print("may the power be with you")
    
    results_df_row <- data.frame(matrix(NA, nrow = 1, ncol = length(results_cols)))
    colnames(results_df_row) <- results_cols
    
    results_df_row$i <- b
    results_df_row$truth <- truth
    
    df <- generate_data_simple(N)
    res_ee <- run_EE_VIM(df, covar)
    res_tmle <- run_TMLE_VIM(df, covar, max.it)
    
    # CVTMLE 
    results_df_row$cvtmle <- res_tmle$coef
    results_df_row$cvtmle_se <- res_tmle$std_err
    results_df_row$cvtmle_lower <- res_tmle$ci_l
    results_df_row$cvtmle_upper <- res_tmle$ci_u
    
    # CVAIPW 
    results_df_row$cvaiptw <- res_ee$coef
    results_df_row$cvaiptw_se <- res_ee$std_err
    results_df_row$cvaiptw_upper <- res_ee$ci_l
    results_df_row$cvaiptw_lower <- res_ee$ci_u
    
    # print(se_cvaiptw)
    # print(se_aiptw)
    return_list <- c('i' = results_df_row$i, 
                     'truth' = results_df_row$truth,
                     'cvtmle' = results_df_row$cvtmle, 
                     'cvtmle_se' = results_df_row$cvtmle_se,
                     'cvtmle_lower' = results_df_row$cvtmle_lower, 
                     'cvtmle_upper' = results_df_row$cvtmle_upper, 
                     'cvaiptw' = results_df_row$cvaiptw , 
                     'cvaiptw_se' = results_df_row$cvaiptw_se, 
                     'cvaiptw_lower' = results_df_row$cvaiptw_lower, 
                     'cvaiptw_upper' = results_df_row$cvaiptw_upper)
    
    return(return_list)
  }
  
  results_df <- as.data.frame(run_bootstrap)
  results_df[, 1] <- sub('.*\\.', '', results_df[, 1])
  
  traceback()
  gc()
  return(results_df)
}
