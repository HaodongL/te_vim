
rm(list = ls())
source(paste0(here::here(),"/R/0_config.R"))

source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
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

# ATE
# Diff in mean
mean(df$Y[which(df$A == 1)]) - mean(df$Y[which(df$A == 0)])


y_l <- 0
y_u <-1

# CATE
ws = c('statin_use')
cv = F
dr = F
max.it = 1e4
Q_bounds = NULL
g_bounds = c(0.025, 0.975)
tau_bounds = NULL
tau_s_bounds = NULL
gamma_s_bounds = NULL
add_tau_sc = F

tic()
df_fit <- fit_para(df = df,
                   sl_Q = sl_Q, 
                   sl_g = sl_g,
                   sl_x = sl_x,
                   ws = ws,
                   dr = dr,
                   Q_bounds = Q_bounds,
                   g_bounds = g_bounds,
                   tau_bounds = tau_bounds,
                   tau_s_bounds = tau_s_bounds,
                   gamma_s_bounds = gamma_s_bounds)
toc()


# saveRDS(df_fit, file = "~/Repo/te_vim/data/df_fit.RDS")


cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

p_cate_t <- plot_cate(cm_names, cbind(df, "tau" = df_fit$mu1_hat - df_fit$mu0_hat))

p_cate_dr <- plot_cate(cm_names, cbind(df, "tau" = df_fit$tau))

# ggsave("tnp/plot/p_cate_t.png", p_cate_t, width = 5, height = 5)
# ggsave("tnp/plot/p_cate_dr.png", p_cate_dr, width = 5, height = 5)

# stratified TMLE
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")
df_strat <- tmle_stratified(df, cm_names)
# saveRDS(df_strat, file = "~/Repo/te_vim/data/df_strat.RDS")

p_cate_strat <- plot_tmle_strat(cm_names, df_strat)

# ggsave("tnp/plot/p_cate_strat_ci.png", p_cate_strat, width = 5, height = 5)

# GRF cate
library(grf)

df_W <- df %>% select(setdiff(names(df), c("Y", "A")))
W <- model.matrix(~. , data=df_W)
Y <- df$Y
A <- df$A

c_forest <- causal_forest(X = W , Y = Y, W = A)
c_forest_hat <- predict(c_forest, estimate.variance = TRUE)
p_cate_grf <- plot_cate(cm_names, 
                        cbind(df, 
                              "tau" = c_forest_hat$predictions,
                              "var_tau" = c_forest_hat$variance.estimates))
# ggsave("tnp/plot/p_cate_grf.png", p_cate_grf, width = 5, height = 5)

p_cate_all <- ggarrange(p_cate_t + ggtitle("SL CATE estimates (T-learner)"), 
                        p_cate_dr + ggtitle("SL CATE estimates (DR-learner)"), 
                        p_cate_grf + ggtitle("GRF CATE estimates"), 
                        p_cate_strat + ggtitle("TMLE Stratified ATE estimates"))
# ggsave("tnp/plot/p_cate_all.png", p_cate_all, width = 10, height = 10)


