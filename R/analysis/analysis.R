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
source(paste0(repo_path, "R/est_function/vim.R"))
source(paste0(repo_path, "R/est_function/vte.R"))
source(paste0(repo_path, "R/est_function/tmle_stratified.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/fit_paraloop.R"))

### ------------  Part 0. eda ------------ ###
# df_w <- read_csv(file = "data/supp/df_w.csv")
# 
# # Table of Baseline Characteristics 
# df_w_summary <- df_w %>% 
#   filter(INSNVFL == FALSE) %>% 
#   select(-c("USUBJID", "INSNVFL")) %>% 
#   mutate(A = ifelse(A == 1, "Liraglutide", "Placebo"))
# 
# df_w_summary <- labelled::remove_labels(df_w_summary)
# tbl <- table1(~. | A, data = df_w_summary)
# tbl

### ------------  Part 1. import data  ------------ ###
# outcome = 'diab'; t = 24
# df <- get_data(outcome, t, rm_baseIns=T)

outcome = 'diab2'; t = 24
df <- get_data(outcome, t)

# df_old <- get_data_old(outcome, t)
# # df_old$EGFREPI <- df$EGFMDRBC
# df <- df_old

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data



# df_old = df
# df_fit_old = df_fit

### ------------  Part 2. Estimation ------------ ###
set.seed(1)
if (outcome == "a1c"){
  sl_Q <- Lrnr_sl$new(
    learners = lrnr_stack_Q,
    metalearner = ls_metalearner
  )
}


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

hist(df_fit$pi_hat)

res_tmle = TMLE_VIM_a(df_fit, y_l, y_u, max_it = 1e4, lr = 1e-4)
res_ee = EE_VIM(df_fit)

# saveRDS(df_fit, file = "~/Repo/te_vim/data/df_fit.RDS")

# visualize CATE
# df_fit <- readRDS("~/Repo/te_vim/data/df_fit.RDS")

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


# VTE
# df_fit$tau <- df_fit$mu1_hat - df_fit$mu0_hat
aipw_vte <- VTE(df_fit, method = "AIPW")
tmle_vte <- VTE(df_fit, method = "TMLE")


# VIM
# tmle vim
theta_TMLE_a <- VIM(df_fit, method = "TMLE_a", y_l = 0, y_u = 1, max.it = 1e4, lr = 1e-3)
theta_TMLE_b <- VIM(df_fit, method = "TMLE_b", y_l = 0, y_u = 1, max.it = 1e4, lr = 1e-3)

tic()
res <- run_VIM_Theta(df = df,
                     sl_Q = sl_Q, 
                     sl_g = sl_g,
                     sl_x = sl_x,
                     ws = c('statin_use'), 
                     cv = F,
                     dr = F,
                     tmle_b = F, 
                     max_it = 1e4)
toc()

# ee vim
# df_fit$Y <- rescale(df_fit$Y, y_l, y_u)
# df_fit$tau <- df_fit$tau*(y_u - y_l)
# df_fit$tau_s <- df_fit$tau_s*(y_u - y_l)
# df_fit$po <- df_fit$po*(y_u - y_l)
theta_EE <- VIM(df_fit, method = "AIPW")

# vim loop over all covariates
set.seed(1)
ws = c("AGE", "SEX", "RACE", "COUNTRY", "SMOKER", "NYHACLAS",
       "DIABDUR", "ANTDBFL", "AHYPERFL", "INCPASSN", 
       "BMIBL", "PULSEBL", "SYSBPBL", "DIABPBL", "HBA1CBL", "HDL1BL", "LDL1BL",
       "CHOL1BL", "RC", "RCoverHDL","TRIG1BL", "CREATBL", "EGFMDRBC", 
       "RETINSEV", "GERDBLFL", "PPIFL", "H2BLFL", 
       "MIFL", "STROKEFL", "REVASFL", "STENFL", "CHDFL", "IHDFL", "CHFFL",
       "KIDFL", "MICFL", "HYPFL", "LVSDFL", "PADFL", "CVRISK", "HBA1CGRN", "DDURGRN", 
       "statin_use", "antihypertensives", "betab", "minera", "adp",
       "vkantag", "caantag", "thiazide", "loopdiur")

# ws = c("statin_use", "antihypertensives")
cv = F
dr = T
max.it = 1e4
Q_bounds = c(1e-4, 1-1e-4)
g_bounds = c(0.025, 0.975)
tau_bounds = c(-1+1e-4, 1-1e-4)
tau_s_bounds = c(-1+1e-4, 1-1e-4)
gamma_s_bounds = c(1e-8, 1-1e-8)
add_tau_sc = F

tic()
df_fit_list <- fit_paraloop(df = df,
                   sl_Q = sl_Q, 
                   sl_g = sl_g,
                   sl_x = sl_x,
                   ws = ws,
                   dr = dr,
                   Q_bounds = Q_bounds,
                   g_bounds = g_bounds,
                   tau_bounds = tau_bounds,
                   tau_s_bounds = tau_s_bounds,
                   gamma_s_bounds = gamma_s_bounds,
                   add_tau_sc = add_tau_sc)
toc()

# saveRDS(df_fit_list, file = "~/Repo/te_vim/data/df_fit_list.RDS")

# temp <- df_fit_list[[1]]
vim_list <- list()

for (i in c(1:length(df_fit_list))){
  theta_TMLE_a <- VIM(df_fit_list[[i]], method = "TMLE_a", y_l = 0, y_u = 1, max.it = 1e4, lr = 1e-3)
  theta_TMLE_b <- VIM(df_fit_list[[i]], method = "TMLE_b", y_l = 0, y_u = 1, max.it = 1e4, lr = 1e-3)
  theta_EE <- VIM(df_fit_list[[i]], method = "AIPW")
  res = list("theta_EE" = theta_EE,
             "theta_TMLE_a" = theta_TMLE_a,
             "theta_TMLE_b" = theta_TMLE_b)
  vim_list[[i]] = res
}

df_theta <- data.frame("varname" = NULL, "importance" = NULL, "ci_l" = NULL, "ci_u" = NULL, "method" = NULL)
for (i in c(1:length(df_fit_list))){
  theta_EE <- vim_list[[i]]$theta_EE
  theta_TMLE_a <- vim_list[[i]]$theta_TMLE_a
  theta_TMLE_b <- vim_list[[i]]$theta_TMLE_b
  
  importance <- c(theta_EE$coef, theta_TMLE_a$coef, theta_TMLE_b$coef)
  ci_l <- c(theta_EE$ci_l, theta_TMLE_a$ci_l, theta_TMLE_b$ci_l)
  ci_u <- c(theta_EE$ci_u, theta_TMLE_a$ci_u, theta_TMLE_b$ci_u)
  varname <- rep(ws[i], 3)
  method <- c(theta_EE$Method, theta_TMLE_a$Method, theta_TMLE_b$Method)
  
  df_theta_i <- data.frame("varname" = varname, "importance" = importance, 
                           "ci_l" = ci_l, "ci_u" = ci_u, "method" = method)
  df_theta <- rbind(df_theta, df_theta_i)
}
# df_theta <- df_theta %>% filter(varname != "INSNVFL")

p_theta_ee <- plot_theta(df_theta, estimator = "AIPW")
p_theta_tmle_a <- plot_theta(df_theta, estimator = "TMLE_a")
p_theta_tmle_b <- plot_theta(df_theta, estimator = "TMLE_b")

p_theta_all <-
  ggarrange(p_theta_ee + ggtitle("EE VIM estimates (DR-learner)"),
            p_theta_tmle_a + ggtitle("TMLE-a VIM estimates (DR-learner)"),
            p_theta_tmle_b + ggtitle("TMLE-b VIM estimates (DR-learner)"))

# ggsave("tnp/plot/p_theta_ee.png", p_theta_ee, width = 7, height = 7)
# ggsave("tnp/plot/p_theta_tmle_a.png", p_theta_tmle_a, width = 7, height = 7)
# ggsave("tnp/plot/p_theta_tmle_b.png", p_theta_tmle_b, width = 7, height = 7)

# ggsave("tnp/plot/p_theta_all.png", p_theta_all, width = 12, height = 12)

plot_theta <- function(df_theta, estimator = "AIPW"){
  p <- ggplot(data=df_theta %>% 
           filter(method == estimator) %>% 
           arrange(importance) %>% 
           mutate(varname = factor(varname, levels = varname)),
         aes(x=varname, y=importance, ymin=ci_l, ymax=ci_u)) +
    geom_pointrange() + 
    # geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
    # geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Variable") + ylab("Importance") +
    theme_bw()  # use a white background
  return(p)
}

df_theta <-  readRDS("~/Repo/te_vim/data/df_vim_dr.RDS")

p_theta_ee <- plot_theta(df_theta, estimator = "EE")
p_theta_tmle_a <- plot_theta(df_theta, estimator = "TMLE-a")
p_theta_ss_hal <- plot_theta(df_theta, estimator = "SS-HAL")

ggsave("tnp/p_theta_ee.png", p_theta_ee, width = 7, height = 7)
ggsave("tnp/p_theta_tmle_a.png", p_theta_tmle_a, width = 7, height = 7)

p_theta_all <-
  ggarrange(p_theta_ee + ggtitle("EE VIM estimates (T-learner)"),
            p_theta_tmle_a + ggtitle("TMLE-a VIM estimates (T-learner)"))

ggsave("tnp/p_theta_all.png", p_theta_all, width = 14, height = 10)

# # aggregating up ---------------------------------------------------------------
# mod = c_forest
# model_hat <- predict(mod, estimate.variance = TRUE)
# model_sigma <- sqrt(model_hat$variance.estimates)
# 
# # saving estimated treatment effect and upper/lower bounds to data frame
# # save these back to the original data frame and calculate confidence intervals
# dat <- df_W
# dat$pred_est <- c(model_hat$predictions)
# dat$pred_var <- c(model_sigma)
# dat$pred_est_lb <- dat$pred_est - 1.96 * dat$pred_var
# dat$pred_est_ub <- dat$pred_est + 1.96 * dat$pred_var
# dat <- dat %>% mutate_if(is.logical, as.factor)
# 
# # for each factor, get the average at each level
# # get results for every level of every variable by aggregating up
# cates <- lapply(names(dat[, c(32:40)]), function(x) {
#   tmp <- dat %>% 
#     group_by(.dots = x) %>% 
#     transmute(
#       variable = x,
#       ate = round(mean(pred_est) * 1, 4),
#       ate_lb = round(mean(pred_est_lb) * 1, 4),
#       ate_ub = round(mean(pred_est_ub) * 1, 4)
#     ) %>% 
#     unique() %>% 
#     as.data.frame()
#   tmp <- tmp[, c(2, 1, 3, 4, 5)]
#   names(tmp)[2] <- "level"
#   tmp
# })
# cates <- do.call(rbind, cates) %>% 
#   mutate_if(is.character, as.factor) 
# 
# ate_grf <- mean(model_hat$predictions)
# 
# # visualize these
# library(ggplot2)
# ggplot(cates, aes(x = level, y = ate, color = variable)) +
#   theme_light() +
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     strip.text.y = element_text(colour = "black"),
#     strip.background = element_rect(colour = NA, fill = NA),
#     legend.position = "none"
#   ) +
#   geom_point() +
#   geom_errorbar(aes(ymin = ate_lb, ymax = ate_ub), width = .2) +
#   geom_hline(yintercept = ate_grf, linetype = 3) +
#   facet_grid(variable ~ ., scales = "free_y") +
#   coord_flip() +
#   labs(y = "CATE")
# 
# average_treatment_effect(mod, subset=(df_W$statin_use == TRUE))
# average_treatment_effect(mod, subset=(df_W$statin_use == F))
# 
# 
# varimp <- variable_importance(c_forest)
# selected_idx = which ( varimp > mean ( varimp ))
# colnames(W)[selected_idx]
