
rm(list = ls())
source(paste0(here::here(),"/R/0_config.R"))

df_em_su_a1c <- readRDS(paste0(here::here(),"/analy_res/fit_em_su_a1c.RDS")) %>% mutate(Outcome = "A1C")
df_em_su_cv <- readRDS(paste0(here::here(),"/analy_res/fit_em_su_cv.RDS")) %>% mutate(Outcome = "CV event")
df_em_su_diab <- readRDS(paste0(here::here(),"/analy_res/fit_em_su_diab.RDS")) %>% mutate(Outcome = "Diabetes intensification")
df_em_su_diab2 <- readRDS(paste0(here::here(),"/analy_res/fit_em_su_diab2.RDS")) %>% mutate(Outcome = "Diabetes intensification - alt definition")


df_plot = bind_rows(df_em_su_diab, df_em_su_diab2, df_em_su_cv, df_em_su_a1c)

  #subset to ate
  df_plot <- df_plot %>% filter(type=="ATE")
  
  #rename parameters
  df_plot <- df_plot %>% mutate(
    param = case_when(param=="E[Y_{A=11}] - E[Y_{A=01}]" ~ "ATE in statin users",
                      param=="E[Y_{A=10}] - E[Y_{A=00}]" ~ "ATE in statin non-users",
                      param=="E[Y_{A=11}] - E[Y_{A=01}] - E[Y_{A=10}] - E[Y_{A=00}]" ~ "Effect modification\n(difference in ATE's)")
  )
  
  
  #Format label widths
  df_plot$Outcome <- str_wrap(df_plot$Outcome, width = 20)

  p_em <- ggplot(df_plot, aes(x = param, y =  tmle_est)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper), width = .2) +
    geom_hline(yintercept = 0, linetype = 3) +
    facet_grid(. ~ Outcome, scales = "free") +
    coord_flip() +
    labs(x="Parameter", y = "Effect size") +
    #ylim(c(ylim_l, ylim_u)) +
    theme_light() +
    theme(
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text = element_text(color = 'black'), 
      strip.background = element_rect(colour = NA, fill = NA),
      legend.position = "none",
      text=element_text(size=8)) 
  
  p_em

  
#Check correct outcome names  
  
save(p_em, file=here::here("figures","em_figures.Rdata"))
  