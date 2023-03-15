library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(foreach)
library(doParallel)
# library(ggplot2)
# library(table1)
# library(ggpubr)

rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/tmle_em.R"))

### ------------  Part 1. import data  ------------ ###
outcome = 'diab'; t = 24
df <- get_data(outcome, t, rm_baseIns=T)

# outcome = 'diab2'; t = 24
# df <- get_data(outcome, t)



# df <- get_data(outcome = 'diab', t = 24, rm_baseIns = FALSE, drop_censor = FALSE)
# df2 <- get_data(outcome = 'diab', t = 24, rm_baseIns = TRUE, drop_censor = FALSE)
# df3 <- get_data(outcome = 'diab', t = 24, rm_baseIns = TRUE, drop_censor = TRUE)




nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data


### ------------  Part 2. Estimation ------------ ###
set.seed(1)

# TMLE effect modification
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

# cm_names <- c("statin_use")
registerDoParallel(10)
tic()
df_em <- tmle_em(df, cm_names)
toc()

saveRDS(df_em, file = "~/Repo/te_vim/data/df_em.RDS")



fit_em <- tmle_em_one(df, c("statin_use"))

fit_em$summary



# df_em <- readRDS("~/Repo/te_vim/data/df_em.RDS")
# p_cate_strat <- plot_tmle_em(cm_names, df_em)




