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
    ylim(c(-0.1, 0)) +
    theme(text=element_text(size=8))
  
  return(p_cate)
}