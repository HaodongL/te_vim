rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))

# Data
outcome = 'diab2'; t = 24
df <- get_data(outcome, t, rm_baseIns=F)

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data


table("Y" = df$Y, "A" = df$A)

# Initial Fit
fit_path <- paste0("~/Repo/te_vim/data/df_fit_", outcome,".RDS")
df_fit <- readRDS(fit_path)

# ATE
res_ate <- TMLE_ATE(df_fit[[1]])
res_ate$coef


# VTE
res_vte_ee <- EE_VTE(data = df_fit[[1]])
res_vte_tmle <- TMLE_VTE(data = df_fit[[1]], max_it = 1e4, lr = 1e-4)

res_vte_ee$coef
res_vte_tmle$coef

# VIM
set.seed(1)
df_vim <- foreach(i = 1:length(list_ws), .combine = 'rbind') %do% {
  
  # vim
  res_ee <- EE_VIM(data = df_fit[[i]])
  res_tmle <- TMLE_VIM(data = df_fit[[i]], max_it = 1e4, lr = 1e-4)
  
  df_vim_i <- data.frame("varname" = c(names(list_ws)[i], names(list_ws)[i]), 
                         "importance" = c(res_ee$coef, res_tmle$coef), 
                         "ci_l" = c(res_ee$ci_l, res_tmle$ci_l), 
                         "ci_u" = c(res_ee$ci_u, res_tmle$ci_u), 
                         "method" = c('EE', 'TMLE'))
  return(df_vim_i)
}

library(ggplot2)
library(ggpubr)
p_ee <- plot_vim(df_vim, "EE", res_vte_ee$coef, res_vte_ee$ci_l, res_vte_ee$ci_u)
p_tmle <- plot_vim(df_vim, "TMLE", res_vte_tmle$coef, res_vte_tmle$ci_l, res_vte_tmle$ci_u)

p_vim <-
  ggarrange(p_ee + ggtitle("EE VIM estimates"),
            p_tmle + ggtitle("TMLE VIM estimates"))

ggsave(paste0(repo_path, "tnp/plot/p_vim.png"), p_vim, width = 12, height = 7)

# Stratified ATE
source(paste0(repo_path, "R/est_function/tmle_stratified.R"))
set.seed(1)
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

registerDoParallel(3)
df_strat <- tmle_stratified(df, cm_names)

p_strat <- plot_tmle_strat(cm_names, df_strat)

# saveRDS(df_strat, file = "~/Repo/te_vim/data/df_strat_diab2.RDS")

# EM
source(paste0(repo_path, "R/est_function/tmle_em.R"))
registerDoParallel(3)
set.seed(1)
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")
df_em <- tmle_em(df, cm_names)


df_em <- readRDS("~/Repo/te_vim/data/df_em.RDS")

p_em <- plot_tmle_em(cm_names, 
                     df_em, 
                     ate_hat = res_ate$coef, 
                     ate_l = res_ate$ci_l,
                     ate_u = res_ate$ci_u)
p_em  

ggsave(paste0(repo_path, "tnp/plot/p_em.png"), p_em)


drug_levels = c(
  `Statin use` = "statin", 
  `Anti-hypertensives` = "antihy",
  `Beta blockers` = "betab",
  `Mineralocorticoid receptor antagonists` = "minera", 
  `ADP receptor inhibitors`= "adp",
  `Vitamin K antagonists` = "vkanta",
  `Ca antagonists` = "caanta", 
  Thiazide = "thiazi",
  `Loop diuretic` = "loopdi")

# saveRDS(df_em, file = "~/Repo/te_vim/data/df_em.RDS")

  
