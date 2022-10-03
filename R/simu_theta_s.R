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
#ncore = 1

#nCoresPerNode <- 1
#nodeNames <-strsplit(Sys.getenv("SLURM_NODELIST"), ",")[[1]]
#machines=rep(nodeNames, each = nCoresPerNode)
#cl = makeCluster(machines, type = "SOCK")


######  calculate the truth 
# set.seed(1234)
# N <- 1e6 #size of generated data
# df <- generate_data_simple(N, print_truth = TRUE) # VIM_Theta_s: 0.686

###### run simu
set.seed(1234)
B <- 2 #rounds of simu
N <- 5e2 #size of data

#registerDoParallel(ncore)
bootstrap_results <- run_all_simu(B = B, N = N, truth = 0.686)


output_filename <- paste0('~/Repo/te_vim/simu_res/theta_s',"theta_s_", N, "_", Sys.Date(),'.csv')

write.csv(bootstrap_results, output_filename)






