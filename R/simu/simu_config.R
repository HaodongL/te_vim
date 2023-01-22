run_all_simu <- function(B, N, truth, 
                         cv = FALSE, 
                         dr = FALSE, 
                         ws = c('X2'), 
                         max_it = 1e4,
                         lr = 1e-4,
                         do_sshal = FALSE,
                         target_para = "Theta"){
  
  results_cols <- c('i', 'truth', 'cvtmle', 'cvtmle_se',
                    'cvtmle_lower', 'cvtmle_upper', 
                    'cvaiptw', 'cvaiptw_se', 'cvaiptw_lower', 
                    'cvaiptw_upper', 'ss')
  
  results_df <- data.frame(matrix(NA, nrow = B, ncol = length(results_cols)))
  colnames(results_df) <- results_cols
  
  run_bootstrap <- foreach(b = 1:B, .combine = 'rbind') %dopar% {
  
    print(paste0("may the power be with you! ", b))
    set.seed(123 + b)
    
    results_df_row <- data.frame(matrix(NA, nrow = 1, ncol = length(results_cols)))
    colnames(results_df_row) <- results_cols
    
    results_df_row$i <- b
    results_df_row$truth <- truth
    
    df <- generate_data_simple(N)
    # df <- generate_data_v2(N)
    
    if (target_para == "Theta"){
      run_fun <- run_VIM
    }else if (target_para == "Psi"){
      run_fun <- run_VIM2
    }else if (target_para == "VTE"){
      run_fun <- run_VTE
    }

    res <- run_fun(df = df, 
                   sl_Q = sl_Q, 
                   sl_g = sl_g,
                   sl_x = sl_x,
                   ws = ws, 
                   cv = cv,
                   dr = dr,
                   max_it = max_it, 
                   lr = lr,
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
    results_df_row$ss_se <- res_ss$std_err
    results_df_row$ss_lower <- res_ss$ci_l
    results_df_row$ss_upper <- res_ss$ci_u
    
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
    
    # plug-in HAL
    if (do_sshal) {
      df_fit_hal <- fit_para_hal(df = df, 
                                 sl_Q = sl_Q_hal, 
                                 sl_g = sl_g_hal,
                                 sl_ws = sl_ws_hal,
                                 ws = ws)
      res_hal <- HAL_VIM(df_fit_hal)
      
      results_df_row$sshal <- res_hal$coef
      results_df_row$sshal_se <- res_hal$std_err
      results_df_row$sshal_lower <- res_hal$ci_l
      results_df_row$sshal_upper <- res_hal$ci_u
    }else{
      results_df_row$sshal <- 0
      results_df_row$sshal_se <- 0
      results_df_row$sshal_lower <- 0
      results_df_row$sshal_upper <- 0
    }
    
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
                     'ss' = results_df_row$ss,
                     'ss_se' = results_df_row$ss_se, 
                     'ss_lower' = results_df_row$ss_lower, 
                     'ss_upper' = results_df_row$ss_upper,
                     'sshal' = results_df_row$sshal , 
                     'sshal_se' = results_df_row$sshal_se, 
                     'sshal_lower' = results_df_row$sshal_lower, 
                     'sshal_upper' = results_df_row$sshal_upper)
    
    return(return_list)
  }
  
  results_df <- as.data.frame(run_bootstrap)
  results_df[, 1] <- sub('.*\\.', '', results_df[, 1])
  
  traceback()
  gc()
  return(results_df)
}


