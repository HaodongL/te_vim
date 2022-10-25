# .libPaths("~/Repo/Rlib_backup/library")

require(tidyverse)
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/vim.R"))

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

# for (N in c(500, 1e3, 2e3,  3e3, 4e3, 5e3, 1e4, 2e4)){
#   print(N)
#   B <- 500 #rounds of simu
# 
#   registerDoParallel(10)
#   tic()
#   bootstrap_results <- run_all_simu(B = B, N = N, cv = F, dr = F, max.it = 1e4, truth = 0.686)
#   toc()
# 
#   output_filename <- paste0('~/Repo/te_vim/simu_res/theta_s/',"local_earth_nocv_t_", N, "_", Sys.Date(),'.csv')
#   write.csv(bootstrap_results, output_filename)
# }



set.seed(1234)
B <- 500 #rounds of simu
N <- 5e2

registerDoParallel(10)
tic()
bootstrap_results <- run_all_simu(B = B, 
                                  N = N, 
                                  truth = 0.686,
                                  cv = F, 
                                  dr = T, 
                                  max.it = 1e4, 
                                  lr = 1e-3, 
                                  tmle_dr_update = T)
toc()
# 
# 
# 
# output_filename <- paste0('~/Repo/te_vim/simu_res/theta_s/',"local_earth_nocv_", N, "_", Sys.Date(),'.csv')
# write.csv(bootstrap_results, output_filename)



# local
res <- bootstrap_results %>% mutate(n = N)
table_results_data <- sum_metric(res)
# 
# res <- res[-which(is.na(res$cvtmle)),]
# 
# temp <- bootstrap_results
# 
# 
# covar = c('X2')
# cv = TRUE
# max.it = 2e3
# Q_bounds = c(0.001, 0.999)
# g_bounds = c(0.025, 0.975)
# tau_bounds = c(-1+1e-3, 1-1e-3)
# tau_s_bounds = c(-1+1e-3, 1-1e-3)
# gamma_s_bounds = c(1e-6, 1-1e-6)
# cate_option = "DR-Learner"
