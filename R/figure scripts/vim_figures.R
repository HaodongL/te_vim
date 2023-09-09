

rm(list = ls())
source(paste0(here::here(),"/R/0_config.R"))

df_vim_t_a1c <- readRDS(paste0(here::here(),"/analy_res/df_vim_t_a1c.RDS"))
df_vim_t_cv <- readRDS(paste0(here::here(),"/analy_res/df_vim_t_cv.RDS"))
df_vim_t_diab <- readRDS(paste0(here::here(),"/analy_res/df_vim_t_diab.RDS"))
df_vim_t_diab2 <- readRDS(paste0(here::here(),"/analy_res/df_vim_t_diab2.RDS"))



df_plot=df_vim_t_diab
covar_levels=covar_levels

plot_vim <- function(df_plot, covar_names=covar_levels){

  df_plot$varname <- fct_recode(df_plot$varname, !!!covar_names)
  
  #Format label widths
  #df_plot$varname <- str_wrap(df_plot$varname, width = 20)
  df_plot$sig <- factor(ifelse(df_plot$ci_l>0 & df_plot$ci_u>0 | df_plot$ci_l<0 & df_plot$ci_u<0, 1, 0), levels=c("1","0"))
  
  plot_theta <- ggplot(data=df_plot %>% filter(method!="SS") %>%
                         arrange(importance) %>% 
                         mutate(varname = factor(varname, levels = unique(varname))),
                       aes(x=varname, y=importance, ymin=ci_l, ymax=ci_u, color=sig, fill=sig, alpha=sig)) +
    geom_pointrange() + 
    facet_grid(~method, scales = "free") +
    geom_hline(yintercept=0, lty=2) + 
    coord_flip() +  
    xlab("Variable") + ylab("Importance") +
    scale_color_manual(values=c("black","grey30")) + 
    scale_fill_manual(values=c("black","grey30")) + 
    scale_alpha_manual(values=c(1,0.5)) + 
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.background = element_rect(colour = NA, fill = NA),
      strip.text.y.left = element_text(angle = 0, color = 'black'), 
      legend.position = "none",
      text=element_text(size=8),
      strip.placement = "outside") 
  
  return(plot_theta)
}

p_vim_diab <- plot_vim(df_vim_t_diab)
p_vim_diab

p_vim_diab2 <- plot_vim(df_vim_t_diab2)
p_vim_diab2

p_vim_cv <- plot_vim(df_vim_t_cv)
p_vim_cv

p_vim_a1c <- plot_vim(df_vim_t_a1c)
p_vim_a1c

# To do: make option to subset to top ten or so
# fix rest of labels like RC

save(p_vim_diab,p_vim_diab2,p_vim_cv,p_vim_a1c, file=here::here("figures","vim_figures.Rdata"))

