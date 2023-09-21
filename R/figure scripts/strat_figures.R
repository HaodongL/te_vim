



rm(list = ls())
source(paste0(here::here(),"/R/0_config.R"))


df_strat_diab <- readRDS(paste0(here::here(),"/analy_res/df_strat_diab.RDS"))
df_strat_diab2 <- readRDS(paste0(here::here(),"/analy_res/df_strat_diab2.RDS"))
df_strat_cv <- readRDS(paste0(here::here(),"/analy_res/df_strat_cv.RDS"))
df_strat_a1c <- readRDS(paste0(here::here(),"/analy_res/df_strat_a1c.RDS"))



#df_w <- read_csv(file = paste0(repo_path, "data/supp/df_w.csv"))
# cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
#               "vkantag", "caantag", "thiazide", "loopdiur")
covar_levels

#------------------------------------------------------------------------------
# Plotting function
#------------------------------------------------------------------------------

df_plot=df_strat_diab
varname=cm_names
unique(df_strat_diab$group)

str_split(df_plot$group, pattern ="_", simplify = T)[,2]

plot_strat_ate <- function(df_plot, covar_levels=drug_levels){
 
  ate <- mean(df_plot$tau)
  ylim_l <- min(df_plot$tau - 1.96*df_plot$se)
  ylim_u <- max(df_plot$tau + 1.96*df_plot$se)
  
  
  df_plot$level <- str_split(df_plot$group, pattern ="_", simplify = T)[,1]
  df_plot <- df_plot %>% mutate(level = case_when(
    level=="T" ~ "Yes", level=="F" ~ "No"
  ))
  df_plot$varname <- str_split(df_plot$group, pattern ="_", simplify = T)[,2]
  df_plot$varname <- fct_recode(df_plot$varname, !!!covar_levels)
  
  #Format label widths
  df_plot$varname <- str_wrap(df_plot$varname, width = 10)
  
  p_cate_strat <- ggplot(df_plot, aes(x = level, y = tau, color = varname)) +
    geom_point() +
    geom_errorbar(aes(ymin = tau - 1.96*se, 
                      ymax = tau + 1.96*se), width = .2) +
    geom_hline(yintercept = 0, linetype = 3) +
    facet_grid(varname ~ ., scales = "free_y", switch="both") +
    coord_flip() +
    labs(x="Baseline drug usage", y = "CATE") +
    ylim(c(ylim_l, ylim_u)) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.background = element_rect(colour = NA, fill = NA),
      strip.text.y.left = element_text(angle = 0, color = 'black'), 
      legend.position = "none",
      text=element_text(size=8),
      strip.placement = "outside"
    ) 
  
  return(p_cate_strat)
}



p_strat_diab <- plot_strat_ate(df_strat_diab)
p_strat_diab

p_strat_diab2 <- plot_strat_ate(df_strat_diab2)
p_strat_diab2

p_strat_cv <- plot_strat_ate(df_strat_cv)
p_strat_cv

p_strat_a1c <- plot_strat_ate(df_strat_a1c)
p_strat_a1c


#To do: add percentages/N's in each level, add a panel for an overall unstratified level
#finish fixing labels

save(p_strat_diab,p_strat_diab2,p_strat_cv,p_strat_a1c, file=here::here("figures","strat_figures.Rdata"))








