run_all_simu <- function(B, N, truth, 
                         cv = FALSE, 
                         dr = TRUE, 
                         lfm_linear = FALSE, 
                         tmle_dr_update = FALSE,
                         ws = c('X2'), 
                         max.it = 1e4,
                         lr = 1e-4){
  results_cols <- c('i', 'truth', 'cvtmle', 'cvtmle_se',
                    'cvtmle_lower', 'cvtmle_upper', 
                    'cvaiptw', 'cvaiptw_se', 'cvaiptw_lower', 
                    'cvaiptw_upper', 'ss')
  
  results_df <- data.frame(matrix(NA, nrow = B, ncol = length(results_cols)))
  colnames(results_df) <- results_cols
  
  run_bootstrap <- foreach(b = 1:B, .combine = 'rbind') %dopar% {
    
    print(paste0("may the power be with you! ", b))
    set.seed(1234 + b)
    
    results_df_row <- data.frame(matrix(NA, nrow = 1, ncol = length(results_cols)))
    colnames(results_df_row) <- results_cols
    
    results_df_row$i <- b
    results_df_row$truth <- truth
    
    df <- generate_data_simple(N)
    # res_ee <- run_EE_VIM(df, ws)
    # res_tmle <- run_TMLE_VIM(df, ws, max.it)
    res <- run_VIM_Theta(df = df, 
                         sl_Q = sl_Q, 
                         sl_g = sl_g,
                         sl_x = sl_x,
                         ws = ws, 
                         cv = cv,
                         dr = dr,
                         lfm_linear = lfm_linear,
                         max.it = max.it, 
                         lr = lr,
                         tmle_dr_update = tmle_dr_update,
                         Q_bounds = c(1e-4, 1-1e-4), 
                         g_bounds = c(0.025, 0.975),
                         tau_bounds = c(-1+1e-4, 1-1e-4),
                         tau_s_bounds = c(-1+1e-4, 1-1e-4),
                         gamma_s_bounds = c(1e-8, 1-1e-8)
                         )
    res_ee <- res$resEE
    res_tmle <- res$resTMLE
    res_ss <- res$resSS
    
    # SS
    results_df_row$ss <- res_ss$coef
    
    # CVTMLE 
    results_df_row$cvtmle <- res_tmle$coef
    results_df_row$cvtmle_se <- res_tmle$std_err
    results_df_row$cvtmle_lower <- res_tmle$ci_l
    results_df_row$cvtmle_upper <- res_tmle$ci_u
    
    # CVAIPW 
    results_df_row$cvaiptw <- res_ee$coef
    results_df_row$cvaiptw_se <- res_ee$std_err
    results_df_row$cvaiptw_lower <- res_ee$ci_l
    results_df_row$cvaiptw_upper <- res_ee$ci_u
    
    
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
                     'cvaiptw_upper' = results_df_row$cvaiptw_upper,
                     'ss' = results_df_row$ss)
    
    return(return_list)
  }
  
  results_df <- as.data.frame(run_bootstrap)
  results_df[, 1] <- sub('.*\\.', '', results_df[, 1])
  
  traceback()
  gc()
  return(results_df)
}



# 
sum_metric <- function(res){
  table_results_data <- res %>%
    dplyr::mutate(cvtmle_proportion = truth <= cvtmle_upper & truth >= cvtmle_lower,
                  cvaiptw_proportion = truth <= cvaiptw_upper & truth >= cvaiptw_lower,
                  
                  cvtmle_widthCI = cvtmle_upper-cvtmle_lower,
                  cvaiptw_widthCI = cvaiptw_upper-cvaiptw_lower,
    ) %>%
    dplyr::group_by(n, truth) %>%
    summarize(cvtmle_coverage = mean(cvtmle_proportion),
              cvaiptw_coverage = mean(cvaiptw_proportion),
              
              cvtmle_bias = mean(cvtmle) - mean(truth),
              cvaiptw_bias = mean(cvaiptw) - mean(truth),
              
              cvtmle_var = var(cvtmle),
              cvaiptw_var = var(cvaiptw),
              
              cvtmle_mse = cvtmle_bias^2 + var(cvtmle),
              cvaiptw_mse = cvaiptw_bias^2 + var(cvaiptw),
              
              
              # 2020-02-01 coverage of oracle CI
              cvtmle_oracle = mean(truth <= cvtmle + 1.96*sd(cvtmle) & truth >= cvtmle - 1.96*sd(cvtmle)),
              cvaiptw_oracle = mean(truth <= cvaiptw + 1.96*sd(cvaiptw) & truth >= cvaiptw - 1.96*sd(cvaiptw)),
              
              cvtmle_meanwidthCI = mean(cvtmle_widthCI),
              cvaiptw_meanwidthCI = mean(cvaiptw_widthCI)
              
    ) %>%  ungroup()
  
  return(table_results_data)
}