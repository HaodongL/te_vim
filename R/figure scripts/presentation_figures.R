library(kableExtra)
library(readr)
library(dplyr)
library(sl3)
library(tictoc)
library(tmle3)
library(ggplot2)
library(table1)
library(ggpubr)
library(R6)



rm(list = ls())
source(paste0(here::here(),"/R/0_config.R"))

source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))
source(paste0(repo_path, "R/est_function/fit_paraloop.R"))
df_w <- read_csv(file = paste0(repo_path, "data/supp/df_w.csv"))
# df_fit <- readRDS("~/Repo/te_vim/data/df_fit.RDS")
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

df_strat <- readRDS(paste0(repo_path, "analy_res/df_strat.RDS"))



# p_cate_strat <- plot_tmle_strat(cm_names, df_strat)
# p_cate_strat

varname=cm_names
df_plot=df_strat
buffer=0.01
  
  ate <- mean(df_plot$tau)
  ylim_l <- min(df_plot$tau - 1.96*df_plot$se)
  ylim_u <- max(df_plot$tau + 1.96*df_plot$se)
  
  df_plot$varname <- rep(varname, each=2)
  df_plot$group[grepl("T_",df_plot$group)] <- "Yes"
  df_plot$group[grepl("F_",df_plot$group)] <- "No"
  
  p_cate <- ggplot(df_plot, aes(x = group, y = tau, color = varname)) +
    theme_light() +
    theme(
      text=element_text(size=8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text.y.left = element_text(size=8,angle = 0, face = "bold", colour = "black"),
      #= element_text(size=10, angle = 15, face = "bold", colour = "black"),
      strip.background = element_rect(colour = NA, fill = NA),
      axis.text.y = element_text(size=8, hjust = 1),
      #strip.text.y = element_text(size=10, angle = 180, face = "bold"),
      strip.placement = "outside",
      axis.text.x = element_text(size=10, vjust = 0.5),
      panel.spacing = unit(0, "lines"),
      legend.position = "none"
    ) +
    geom_point() +
    geom_errorbar(aes(ymin = pmax(tau - 1.96*se, ylim_l), 
                      ymax = pmin(tau + 1.96*se, ylim_u)), width = .2) +
    geom_hline(yintercept = 0, linetype = 3) +
    facet_grid(varname ~ ., scales = "free_y", switch="both") +
    coord_flip() +
    labs(y = "CATE", x="Subgroup variable") +
    ylim(c(ylim_l, ylim_u)) 

  p_cate
  
  #XXXXX TO DO:
  #Fix angle and labels of strips




## 2. TMLE Effect Modification
df_em <-readRDS(paste0(repo_path, "analy_res/df_em.RDS"))
p_cate_strat <- plot_tmle_em(cm_names, df_em)
p_cate_strat


## 3. VTE Estimation
res_vte <-readRDS(paste0(repo_path, "analy_res/res_vte_t.RDS"))

print(paste0("VTE TMLE:"))
print(paste0("coef:", res_vte$res_tmle$coef))
print(paste0("ci low:", res_vte$res_tmle$ci_l))
print(paste0("ci up:", res_vte$res_tmle$ci_u))


print(paste0("VTE EE:"))
print(paste0("coef:", res_vte$res_ee$coef))
print(paste0("ci low:", res_vte$res_ee$ci_l))

df_theta <-readRDS(paste0(repo_path, "analy_res/df_vim_t.RDS"))


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

p_theta_ee <- plot_theta(df_theta, estimator = "EE")
p_theta_tmle <- plot_theta(df_theta, estimator = "TMLE")
p_theta_ss <- plot_theta(df_theta, estimator = "SS")

p_theta_all <-
  ggarrange(p_theta_ee + ggtitle("EE VIM estimates (T-learner)"),
            p_theta_tmle + ggtitle("TMLE VIM estimates (T-learner)"),
            p_theta_ss + ggtitle("SS VIM estimates (T-learner)"))

p_theta_all















