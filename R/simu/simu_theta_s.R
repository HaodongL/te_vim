
rm(list = ls())
require(tidyverse)
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/vim.R"))
source(paste0(repo_path, "R/est_function/fit_hal.R"))

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
                                    truth = 0.685,
                                    do_sshal = F,
                                    target_para = "Theta")
  toc()

  output_filename <- paste0('~/Repo/te_vim/simu_res/theta_s/',"hal_t_", N, "_", Sys.Date(),'.csv')
  write.csv(bootstrap_results, output_filename)
}


# Check results
# bootstrap_results = read_csv(paste0('~/Repo/te_vim/simu_res/theta_s/',"hal_dr_", N, "_", '2022-12-27','.csv'))
# res <- bootstrap_results %>% mutate(n = N)
# table_results_data <- sum_metric(res)

# res <- read_hp(filename = 'earth_t_', savedate = '2023-01-12')
# data_long <- proc_df_tbl(res, N_len = 7)



