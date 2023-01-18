library(data.table)
read_hp <- function(filename, 
                    savedate, 
                    Ns = c(2e2, 5e2, 1e3, 3e3, 5e3, 1e4, 2e4),
                    ss = FALSE,
                    target_para = "theta_s"){
  res <- data.frame()
  for (N in Ns){
    output_filename <- paste0('~/Repo/te_vim/simu_res/',target_para,'/',filename, N, "_", savedate,'.csv')
    res1 <- read_csv(output_filename) %>% mutate(n = N)
    res <- rbind(res, res1)
  }
  return(res)
}



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


proc_df_tbl <- function(res, N_len = 8){
  
  table_results_data <- sum_metric(res)
  
  #wide to long
  data_long_coverage <- gather(table_results_data, method, coverage, cvtmle_coverage:cvaiptw_coverage, factor_key=TRUE) %>%select(c(n, truth, method, coverage))
  
  data_long_oracle <- gather(table_results_data, method, oracle_coverage, cvtmle_oracle:cvaiptw_oracle, factor_key=TRUE) %>% select(c(oracle_coverage))
  
  data_long_bias <- gather(table_results_data, method, bias, cvtmle_bias:cvaiptw_bias, factor_key=TRUE) %>% select(c(bias))
  
  data_long_var<- gather(table_results_data, method, var, cvtmle_var:cvaiptw_var, factor_key=TRUE) %>%select(c(var))
  
  data_long_mse <- gather(table_results_data, method, mse, cvtmle_mse:cvaiptw_mse, factor_key=TRUE) %>%select(c(mse))
  
  data_long_CIwidth <- gather(table_results_data, method, CIwidth, cvtmle_meanwidthCI:cvaiptw_meanwidthCI, factor_key=TRUE) %>%select(c(CIwidth))
  
  data_long = cbind(data_long_coverage,
                    data_long_oracle, 
                    data_long_bias,
                    data_long_var,
                    data_long_mse,
                    data_long_CIwidth) %>% 
    select(c(n, method,truth,var,bias,mse,coverage,oracle_coverage, CIwidth))
  
  setnames(data_long, old = c('n' ,'method','truth','var','bias','mse','coverage','oracle_coverage','CIwidth'), 
           new = c ('n',"Method",'True_Theta','Variance','Bias','MSE','Coverage','Coverage_or','CI_width'))
  
  data_long$Method = c(rep('TMLE',N_len), rep('EE',N_len))
  data_long <- data_long %>% arrange(n)
  return(data_long)
}






# read ss hal res
sum_metric <- function(res){
  table_results_data <- res %>%
    dplyr::mutate(cvtmle_proportion = truth <= cvtmle_upper & truth >= cvtmle_lower,
                  cvaiptw_proportion = truth <= cvaiptw_upper & truth >= cvaiptw_lower,
                  ss_proportion = truth <= ss_upper & truth >= ss_lower,
                  sshal_proportion = truth <= sshal_upper & truth >= sshal_lower,


                  cvtmle_widthCI = cvtmle_upper-cvtmle_lower,
                  cvaiptw_widthCI = cvaiptw_upper-cvaiptw_lower,
                  ss_widthCI = ss_upper-ss_lower,
                  sshal_widthCI = sshal_upper-sshal_lower,
    ) %>%
    dplyr::group_by(n, truth) %>%
    summarize(cvtmle_coverage = mean(cvtmle_proportion),
              cvaiptw_coverage = mean(cvaiptw_proportion),
              ss_coverage = mean(ss_proportion),
              sshal_coverage = mean(sshal_proportion),

              cvtmle_bias = mean(cvtmle) - mean(truth),
              cvaiptw_bias = mean(cvaiptw) - mean(truth),
              ss_bias = mean(ss) - mean(truth),
              sshal_bias = mean(sshal) - mean(truth),

              cvtmle_var = var(cvtmle),
              cvaiptw_var = var(cvaiptw),
              ss_var = var(ss),
              sshal_var = var(sshal),

              cvtmle_mse = cvtmle_bias^2 + var(cvtmle),
              cvaiptw_mse = cvaiptw_bias^2 + var(cvaiptw),
              ss_mse = ss_bias^2 + var(ss),
              sshal_mse = sshal_bias^2 + var(sshal),

              # 2020-02-01 coverage of oracle CI
              cvtmle_oracle = mean(truth <= cvtmle + 1.96*sd(cvtmle) & truth >= cvtmle - 1.96*sd(cvtmle)),
              cvaiptw_oracle = mean(truth <= cvaiptw + 1.96*sd(cvaiptw) & truth >= cvaiptw - 1.96*sd(cvaiptw)),
              ss_oracle = mean(truth <= ss + 1.96*sd(ss) & truth >= ss - 1.96*sd(ss)),
              sshal_oracle = mean(truth <= sshal + 1.96*sd(sshal) & truth >= sshal - 1.96*sd(sshal)),

              cvtmle_meanwidthCI = mean(cvtmle_widthCI),
              cvaiptw_meanwidthCI = mean(cvaiptw_widthCI),
              ss_meanwidthCI = mean(ss_widthCI),
              sshal_meanwidthCI = mean(sshal_widthCI)

    ) %>%  ungroup()

  return(table_results_data)
}


proc_df_tbl <- function(res, N_len = 7){

  table_results_data <- sum_metric(res)

  #wide to long
  data_long_coverage <- gather(table_results_data, method, coverage, cvtmle_coverage:sshal_coverage, factor_key=TRUE) %>%select(c(n, truth, method, coverage))

  data_long_oracle <- gather(table_results_data, method, oracle_coverage, cvtmle_oracle:sshal_oracle, factor_key=TRUE) %>% select(c(oracle_coverage))

  data_long_bias <- gather(table_results_data, method, bias, cvtmle_bias:sshal_bias, factor_key=TRUE) %>% select(c(bias))

  data_long_var<- gather(table_results_data, method, var, cvtmle_var:sshal_var, factor_key=TRUE) %>%select(c(var))

  data_long_mse <- gather(table_results_data, method, mse, cvtmle_mse:sshal_mse, factor_key=TRUE) %>%select(c(mse))

  data_long_CIwidth <- gather(table_results_data, method, CIwidth, cvtmle_meanwidthCI:sshal_meanwidthCI, factor_key=TRUE) %>%select(c(CIwidth))

  data_long = cbind(data_long_coverage,
                    data_long_oracle,
                    data_long_bias,
                    data_long_var,
                    data_long_mse,
                    data_long_CIwidth) %>%
    select(c(n, method,truth,var,bias,mse,coverage,oracle_coverage, CIwidth))

  setnames(data_long, old = c('n' ,'method','truth','var','bias','mse','coverage','oracle_coverage','CIwidth'),
           new = c ('n',"Method",'True_Theta','Variance','Bias','MSE','Coverage','Coverage_or','CI_width'))

  data_long$Method = c(rep('TMLE',N_len), rep('EE',N_len), rep('SS',N_len), rep('SS-HAL',N_len))
  data_long <- data_long %>% arrange(n)
  return(data_long)
}
