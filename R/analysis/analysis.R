library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(ggplot2)

rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))
source(paste0(repo_path, "R/est_function/vte.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))


### ------------  Part 1. import data  ------------ ###
outcome = 'diab'; t = 24
df <- get_data(outcome, t)

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data


### ------------  Part 2. Estimation ------------ ###
set.seed(123)
if (outcome == "a1c"){
  sl_Q <- Lrnr_sl$new(
    learners = lrnr_stack_Q,
    metalearner = ls_metalearner
  )
}


# ATE
# Diff in mean
mean(df$Y[which(df$A == 1)]) - mean(df$Y[which(df$A == 0)])


y_l <- min(df$Y)
y_u <- max(df$Y)
df$Y <- scale01(df$Y, y_l, y_u)

# CATE
ws = c('statin_use')
cv = F
dr = T
max.it = 1e4
Q_bounds = c(1e-4, 1-1e-4)
g_bounds = c(0.025, 0.975)
tau_bounds = c(-1+1e-4, 1-1e-4)
tau_s_bounds = c(-1+1e-4, 1-1e-4)
gamma_s_bounds = c(1e-8, 1-1e-8)
add_tau_sc = T

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
                   gamma_s_bounds = gamma_s_bounds,
                   add_tau_sc = add_tau_sc)
toc()

# visualize CATE
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

p_cate <- plot_cate(cm_names, cbind(df, "tau" = df_fit$tau))

# ggsave("tnp/plot/p_cate.png", p_cate)


# VIM
# tmle vim
theta_TMLE_a <- VIM(df_fit, method = "TMLE_a", y_l = 0, y_u = 1, max.it)
theta_TMLE_b <- VIM(df_fit, method = "TMLE_b", y_l = 0, y_u = 1, max.it)

# ee vim
df_fit$Y <- rescale(df_fit$Y, y_l, y_u)
df_fit$tau <- df_fit$tau*(y_u - y_l)
df_fit$tau_s <- df_fit$tau_s*(y_u - y_l)
df_fit$po <- df_fit$po*(y_u - y_l)

theta_EE <- VIM(df_fit, method = "AIPW")


# VTE
aipw_vte <- VTE(df_fit, method = "AIPW")
tmle_vte <- VTE(df_fit, method = "TMLE")


# GRF cate
library(grf)

df_W <- df %>% select(setdiff(names(df), c("Y", "A")))
W <- model.matrix(~. , data=df_W)
Y <- df$Y
A <- df$A

c_forest <- causal_forest(X = W , Y = Y, W = A)

p_cate_grf <- plot_cate(cm_names, 
                        cbind(df, "tau" = predict(c_forest)$predictions))

# ggsave("tnp/plot/p_cate_grf.png", p_cate_grf)

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
