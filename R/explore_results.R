

rm(list = ls())
source(paste0(here::here(),"/R/0_config.R"))




res_vte_t_a1c <- readRDS(paste0(here::here(),"/analy_res/res_vte_t_a1c.RDS"))
res_vte_t_cv <- readRDS(paste0(here::here(),"/analy_res/res_vte_t_cv.RDS"))
res_vte_t_diab <- readRDS(paste0(here::here(),"/analy_res/res_vte_t_diab.RDS"))
res_vte_t_diab2 <- readRDS(paste0(here::here(),"/analy_res/res_vte_t_diab2.RDS"))

df_strat_diab <- readRDS(paste0(here::here(),"/analy_res/df_strat_diab.RDS"))
df_strat_diab2 <- readRDS(paste0(here::here(),"/analy_res/df_strat_diab2.RDS"))
df_strat_cv <- readRDS(paste0(here::here(),"/analy_res/df_strat_cv.RDS"))
df_strat_a1c <- readRDS(paste0(here::here(),"/analy_res/df_strat_a1c.RDS"))



df_vim_t_a1c <- readRDS(paste0(here::here(),"/analy_res/df_vim_t_a1c.RDS"))
df_vim_t_cv <- readRDS(paste0(here::here(),"/analy_res/df_vim_t_cv.RDS"))
df_vim_t_diab <- readRDS(paste0(here::here(),"/analy_res/df_vim_t_diab.RDS"))
df_vim_t_diab2 <- readRDS(paste0(here::here(),"/analy_res/df_vim_t_diab2.RDS"))

fit_em_su_a1c <- readRDS(paste0(here::here(),"/analy_res/fit_em_su_a1c.RDS"))
fit_em_su_cv <- readRDS(paste0(here::here(),"/analy_res/fit_em_su_cv.RDS"))
fit_em_su_diab <- readRDS(paste0(here::here(),"/analy_res/fit_em_su_diab.RDS"))
fit_em_su_diab2 <- readRDS(paste0(here::here(),"/analy_res/fit_em_su_diab2.RDS"))
