library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(ggplot2)
library(table1)
library(ggpubr)

rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/tmle_stratified.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))

### ------------  Part 1. import data  ------------ ###
outcome = 'diab'; t = 24
df <- get_data(outcome, t, rm_baseIns=T)

# outcome = 'diab2'; t = 24
# df <- get_data(outcome, t)

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data


### ------------  Part 2. Estimation ------------ ###
set.seed(1)

# stratified TMLE
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

registerDoParallel(9)
df_strat <- tmle_stratified(df, cm_names)
saveRDS(df_strat, file = "~/Repo/te_vim/data/df_strat.RDS")


# p_cate_strat <- plot_tmle_strat(cm_names, df_strat)
# ggsave("tnp/plot/p_cate_strat_ci.png", p_cate_strat, width = 5, height = 5)