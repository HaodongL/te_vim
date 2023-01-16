run_all_simu <- function(B, N, truth, 
                         cv = FALSE, 
                         dr = FALSE, 
                         tmle_b = FALSE,
                         ws = c('X2'), 
                         max_it = 1e4,
                         lr = 1e-4,
                         do_sshal = FALSE){
  results_cols <- c('i', 'truth', 'cvtmle', 'cvtmle_se',
                    'cvtmle_lower', 'cvtmle_upper', 
                    'cvaiptw', 'cvaiptw_se', 'cvaiptw_lower', 
                    'cvaiptw_upper', 'ss')
  
  results_df <- data.frame(matrix(NA, nrow = B, ncol = length(results_cols)))
  colnames(results_df) <- results_cols
  
  run_bootstrap <- foreach(b = 1:B, .combine = 'rbind') %dopar% {
    
    # for (b in c(208:500)) {
    
    print(paste0("may the power be with you! ", b))
    set.seed(123 + b)
    
    results_df_row <- data.frame(matrix(NA, nrow = 1, ncol = length(results_cols)))
    colnames(results_df_row) <- results_cols
    
    results_df_row$i <- b
    results_df_row$truth <- truth
    
    df <- generate_data_simple(N)
    # df <- generate_data_v2(N)
    
    res <- run_VIM2(df = df,
                    sl_Q = sl_Q, 
                    sl_g = sl_g,
                    sl_x = sl_x,
                    ws = ws, 
                    cv = cv,
                    dr = dr,
                    tmle_b = tmle_b, 
                    max_it = max_it,
                    lr = lr)
    
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
    
    # }
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
                     'ss' = results_df_row$ss,
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



require(tidyverse)
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))
source(paste0(repo_path, "R/est_function/vim2.R"))

library('foreach')
library('doParallel')
library('tictoc')

###### set up the number of cores for parallel loop
# parallel
uname <- system("uname -a", intern = TRUE)
os <- sub(" .*", "", uname)
if(os=="Darwin"){
  cpus_logical <- as.numeric(system("sysctl -n hw.logicalcpu", intern = TRUE))
} else if(os=="Linux"){
  cpus_logical <- system("lscpu | grep '^CPU(s)'", intern = TRUE)
  cpus_logical <- as.numeric(gsub("^.*:[[:blank:]]*","", cpus_logical))
} else {
  stop("unsupported OS")
}
ncore <- floor(cpus_logical/2)


######  calculate the truth 
# set.seed(1234)
# N <- 1e6 #size of generated data
# df <- generate_data_simple(N, print_truth = TRUE) # VIM_Theta_s: 0.686

###### run simu
for (N in c(2e2, 5e2, 1e3, 3e3, 5e3, 1e4, 2e4)){
  print(N)
  B <- 500 #rounds of simu
  registerDoParallel(10)
  tic()
  bootstrap_results <- run_all_simu(B = B, 
                                    N = N, 
                                    cv = F, 
                                    dr = F, 
                                    tmle_b = F, 
                                    truth = 0.68,
                                    do_sshal = F)
  toc()
  
  output_filename <- paste0('~/Repo/te_vim/simu_res/psi_s/',"hal_t_", N, "_", Sys.Date(),'.csv')
  write.csv(bootstrap_results, output_filename)
}