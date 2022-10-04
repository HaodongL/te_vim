.libPaths(c(.libPaths(), 
            "~/Repo/Rlib_backup/library"))

require(tidyverse)
repo_path = "/Users/haodongli/Repo/te_vim/"
source(paste0(repo_path, "R/fit_sl3.R"))
source(paste0(repo_path, "R/example_helpers.R")) #Used for the current examples
source(paste0(repo_path, "R/fit_cate.R"))
source(paste0(repo_path, "R/vim.R"))
source(paste0(repo_path, "R/simu_config.R"))

library('foreach')
library('doParallel')


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
set.seed(1234)
B <- 500 #rounds of simu
N <- 1e4 #size of data

registerDoParallel(5)
bootstrap_results <- run_all_simu(B = B, N = N, truth = 0.686)


output_filename <- paste0('~/Repo/te_vim/simu_res/theta_s/',"local_", N, "_", Sys.Date(),'.csv')
write.csv(bootstrap_results, output_filename)




# local


