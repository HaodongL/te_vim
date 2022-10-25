library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)

rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/example_helpers.R")) #Used for the current examples
source(paste0(repo_path, "R/sl3_config.R"))
source(paste0(repo_path, "R/fit_para.R"))
source(paste0(repo_path, "R/vim.R"))
source(paste0(repo_path, "R/tmle_v1.R"))
source(paste0(repo_path, "R/vte.R"))
source(paste0(repo_path, "R/aipw.R"))
source(paste0(repo_path, "R/tmle.R"))


df <- read_csv("data/df_all.csv")

cutoff_t <- 24

y_names <- c("AVAL_FRINSLTM", "EVENT_FRINSLTM")
df <- df %>% mutate_at(y_names, ~ifelse(is.na(.), 0, .))


# primary Y
df <- df %>% mutate(Y = ifelse((AVAL_FRINSLTM <= cutoff_t & EVENT_FRINSLTM == 1 )|
                               (AVAL_OADTM <= cutoff_t & EVENT_OADTM == 1 ),
                                1, 0)) %>% 
            select(-c("USUBJID", "AVAL_MACE",
                      "AVAL_NFMI", "AVAL_CVDEATH",
                      "AVAL_NFSTROKE", "AVAL_NONCVDEATH",
                      "EVENT_MACE", "EVENT_NFMI",
                      "EVENT_CVDEATH", "EVENT_NFSTROKE", "EVENT_NONCVDEATH",
                      "AVAL_FRINSLTM", "EVENT_FRINSLTM", 
                      "AVAL_OADTM", "EVENT_OADTM"))

# secondary Y
# df <- df %>% mutate(Y = ifelse((AVAL_MACE <= cutoff_t & EVENT_MACE == 1 )|
#                                (AVAL_NFMI <= cutoff_t & EVENT_NFMI == 1 )|
#                                (AVAL_CVDEATH <= cutoff_t & EVENT_CVDEATH == 1 )|
#                                (AVAL_NFSTROKE <= cutoff_t & EVENT_NFSTROKE == 1 )|
#                                (AVAL_NONCVDEATH <= cutoff_t & EVENT_NONCVDEATH == 1 ), 
#                                1, 0)) %>% 
#              select(-c("USUBJID", "AVAL_MACE",
#                        "AVAL_NFMI", "AVAL_CVDEATH",
#                        "AVAL_NFSTROKE", "AVAL_NONCVDEATH",
#                        "EVENT_MACE", "EVENT_NFMI",
#                        "EVENT_CVDEATH", "EVENT_NFSTROKE", "EVENT_NONCVDEATH",
#                        "AVAL_FRINSLTM", "EVENT_FRINSLTM", 
#                        "AVAL_OADTM", "EVENT_OADTM"))

df <- df %>% mutate(across(where(is.character), ~ as.factor(.)))
df <- df %>% mutate(A = ifelse(ACTARM == "Placebo", 0, 1)) %>% select(-"ACTARM")

cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "minera_cm", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

df <- df %>% mutate_at(cm_names, ~ifelse(is.na(.), FALSE, TRUE))

nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

df <- process_missing(df, nodes)$data
# sapply(df, class)
# sum(df$Y)


ws = c('statin_use')
cv = F
dr = TRUE
lfm_linear = FALSE
max.it = 1e4
Q_bounds = c(0.001, 0.999)
g_bounds = c(0.025, 0.975)
tau_bounds = c(-1+1e-3, 1-1e-3)
tau_s_bounds = c(-1+1e-3, 1-1e-3)
gamma_s_bounds = c(1e-6, 1-1e-6)

# ATE
# Diff in mean
mean(df$Y[which(df$A == 1)]) - mean(df$Y[which(df$A == 0)])

# CATE
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


df_plot <- cbind(df_fit, "statin_use" = df$statin_use)
ggplot(df_plot, 
       aes(x=statin_use, y=tau)) + 
  geom_boxplot()


df_plot <- cbind(df_fit, "SEX" = df$SEX)
ggplot(df_plot, 
       aes(x=SEX, y=tau)) + 
  geom_boxplot()


df_plot <- cbind(df_fit, "antihypertensives" = df$antihypertensives)
ggplot(df_plot, 
       aes(x=antihypertensives, y=tau)) + 
  geom_boxplot()

# VTE
aipw_vte <- VTE(df_fit,method="AIPW")
tmle_vte <- VTE(df_fit,method="TMLE")

# VIM
tic()
res <- run_VIM_Theta(df = df,
                     sl_Q = sl_Q, 
                     sl_g = sl_g,
                     sl_x = sl_x,
                     ws = c('betab'), 
                     cv = F,
                     dr = TRUE,
                     lfm_linear = FALSE, 
                     max.it = 1e4, 
                     Q_bounds = c(0.001, 0.999), 
                     g_bounds = c(0.025, 0.975),
                     tau_bounds = c(-1+1e-3, 1-1e-3),
                     tau_s_bounds = c(-1+1e-3, 1-1e-3),
                     gamma_s_bounds = c(1e-6, 1-1e-6)
)
toc()
res_ee <- res$resEE
res_tmle <- res$resTMLE


# GRF
library(grf)

df_W <- df %>% select(setdiff(names(df), c("Y", "A")))
W <- model.matrix(~. , data=df_W)
Y <- df$Y
A <- df$A

c_forest <- causal_forest(X = W , Y = Y, W = A)



# aggregating up ---------------------------------------------------------------
mod = c_forest
model_hat <- predict(mod, estimate.variance = TRUE)
model_sigma <- sqrt(model_hat$variance.estimates)

# saving estimated treatment effect and upper/lower bounds to data frame
# save these back to the original data frame and calculate confidence intervals
dat <- df_W
dat$pred_est <- c(model_hat$predictions)
dat$pred_var <- c(model_sigma)
dat$pred_est_lb <- dat$pred_est - 1.96 * dat$pred_var
dat$pred_est_ub <- dat$pred_est + 1.96 * dat$pred_var
dat <- dat %>% mutate_if(is.logical, as.factor)

# for each factor, get the average at each level
# get results for every level of every variable by aggregating up
library(dplyr)
cates <- lapply(names(dat[, c(1:2, 34:43)]), function(x) {
  tmp <- dat %>% 
    group_by(.dots = x) %>% 
    transmute(
      variable = x,
      ate = round(mean(pred_est) * 1, 2),
      ate_lb = round(mean(pred_est_lb) * 1, 2),
      ate_ub = round(mean(pred_est_ub) * 1, 2)
    ) %>% 
    unique() %>% 
    as.data.frame()
  tmp <- tmp[, c(2, 1, 3, 4, 5)]
  names(tmp)[2] <- "level"
  tmp
})
cates <- do.call(rbind, cates) %>% 
  mutate_if(is.character, as.factor) 
  

# visualize these
library(ggplot2)
ggplot(cates, aes(x = level, y = ate, color = variable)) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y = element_text(colour = "black"),
    strip.background = element_rect(colour = NA, fill = NA),
    legend.position = "none"
  ) +
  geom_point() +
  geom_errorbar(aes(ymin = ate_lb, ymax = ate_ub), width = .2) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_grid(variable ~ ., scales = "free_y") +
  coord_flip()

average_treatment_effect(mod, subset=(df_W$statin_use == TRUE))
average_treatment_effect(mod, subset=(df_W$statin_use == F))


varimp <- variable_importance(c_forest)
selected_idx = which ( varimp > mean ( varimp ))
colnames(W)[selected_idx]
