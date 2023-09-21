buffer <- 0.2
# helper function, plot cate est by subgroups
plot_cate <- function(varname, 
                      df_plot){
  
  cates <- lapply(varname, function(x) {
    tmp <- df_plot %>% 
      group_by_at(x) %>% 
      transmute(
        variable = x,
        cate = round(mean(tau) * 1, 4)
      ) %>% 
      unique() %>% 
      as.data.frame()
    tmp <- tmp[, c(2, 1, 3)]
    names(tmp)[2] <- "level"
    tmp
  })
  
  cates <- do.call(rbind, cates) %>%  mutate_if(is.character, as.factor) 
  
  ate <- mean(df_plot$tau)
  
  p_cate <- ggplot(cates, aes(x = level, y = cate, color = variable)) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text.y = element_text(colour = "black"),
      strip.background = element_rect(colour = NA, fill = NA),
      legend.position = "none"
    ) +
    geom_point() +
    geom_hline(yintercept = ate, linetype = 3) +
    facet_grid(variable ~ ., scales = "free_y") +
    coord_flip() +
    labs(y = "CATE") +
    ylim(c(ate - buffer, ate + buffer)) +
    theme(text=element_text(size=8))
  
  return(p_cate)
}


# varname <- cm_names
# df_plot <- df_strat

plot_tmle_strat <- function(varname, 
                      df_plot){
  
  ate <- mean(df_plot$tau)
  ylim_l <- ate - buffer
  ylim_u <- ate + buffer
  
  df_plot$varname <- rep(varname, each=2)
  
  p_cate <- ggplot(df_plot, aes(x = group, y = tau, color = varname)) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text.y = element_text(colour = "black"),
      strip.background = element_rect(colour = NA, fill = NA),
      legend.position = "none"
    ) +
    geom_point() +
    geom_errorbar(aes(ymin = pmax(tau - 1.96*se, ylim_l), 
                      ymax = pmin(tau + 1.96*se, ylim_u)), width = .2) +
    geom_hline(yintercept = 0, linetype = 3) +
    facet_grid(varname ~ ., scales = "free_y") +
    coord_flip() +
    labs(y = "CATE") +
    ylim(c(ylim_l, ylim_u)) +
    theme(text=element_text(size=8))
  
  return(p_cate)
}


plot_tmle_em <- function(varname, 
                            df_plot){
  
  ate <- mean(df_plot$tau)
  ylim_l <- ate - buffer
  ylim_u <- ate + buffer
  
  df_plot$varname <- rep(varname, each=3)
  
  p_cate <- ggplot(df_plot, aes(x = group, y = tau, color = varname)) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text.y = element_text(colour = "black"),
      strip.background = element_rect(colour = NA, fill = NA),
      legend.position = "none"
    ) +
    geom_point() +
    geom_errorbar(aes(ymin = pmax(tau - 1.96*se, ylim_l), 
                      ymax = pmin(tau + 1.96*se, ylim_u)), width = .2) +
    geom_hline(yintercept = 0, linetype = 3) +
    facet_grid(varname ~ ., scales = "free_y") +
    coord_flip() +
    labs(y = "CATE") +
    ylim(c(ylim_l, ylim_u)) +
    theme(text=element_text(size=8))
  
  return(p_cate)
}

# varname <- cm_names
# df_plot <- cbind(df, 
#       "tau" = c_forest_hat$predictions,
#       "var_tau" = c_forest_hat$variance.estimates)

plot_cate_grf <- function(varname, 
                      df_plot){
  
  cates <- lapply(varname, function(x) {
    tmp <- df_plot %>% 
      group_by_at(x) %>% 
      transmute(
        variable = x,
        cate = round(mean(tau), 4),
        cate_lb = round(mean(tau) - 1.96*sqrt(mean(var_tau)/n()), 4),
        cate_ub = round(mean(tau) + 1.96*sqrt(mean(var_tau)/n()), 4),
        # cate_lb = round(mean(tau - 1.96*sqrt(var_tau)), 4),
        # cate_ub = round(mean(tau + 1.96*sqrt(var_tau)), 4),
        group_n = n()
      ) %>% 
      unique() %>% 
      as.data.frame()
    tmp <- tmp[, c(2, 1, 3, 4, 5,6)]
    names(tmp)[2] <- "level"
    tmp
  })
  
  cates <- do.call(rbind, cates) %>%  mutate_if(is.character, as.factor) 
  
  ate <- mean(df_plot$tau)
  
  p_cate <- ggplot(cates, aes(x = level, y = cate, color = variable)) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text.y = element_text(colour = "black"),
      strip.background = element_rect(colour = NA, fill = NA),
      legend.position = "none"
    ) +
    geom_point() +
    geom_errorbar(aes(ymin = cate_lb, ymax = cate_ub), width = .2) +
    geom_hline(yintercept = ate, linetype = 3) +
    facet_grid(variable ~ ., scales = "free_y") +
    coord_flip() +
    labs(y = "CATE") +
    ylim(c(ate - buffer, ate + buffer)) +
    theme(text=element_text(size=8))
  
  return(p_cate)
}


plot_theta <- function(df_theta, estimator = "AIPW"){
  p <- ggplot(data=df_theta %>% 
                filter(method == estimator) %>% 
                arrange(importance) %>% 
                mutate(varname = factor(varname, levels = varname)),
              aes(x=varname, y=importance, ymin=ci_l, ymax=ci_u)) +
    geom_pointrange() + 
    # geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
    # geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    geom_hline(yintercept=0, colour = 'red') +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("") + ylab("Importance") + # ylim(c(-350,350)) +
    theme_bw()  # use a white background
  return(p)
}
