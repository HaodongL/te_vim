library(ggpubr)
library(ggh4x)
library(stringr)
library(forcats)
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


plot_tmle_em <- function(varname, df_plot, ate_hat, ate_l, ate_u){
  ate <- ate_hat
  
  ylim_l <- ate - buffer
  ylim_u <- ate + buffer
  
  df_plot$varname <- rep(varname, each=3)
  df_plot$level <- str_split(df_plot$group, pattern ="_", simplify = T)[,2]
  df_plot <- df_plot %>% mutate(level = case_when(
    level=="T" ~ "Yes", level=="F" ~ "No", level=="EM" ~"EM"
  ))
  df_plot$varname <- str_split(df_plot$group, pattern ="_", simplify = T)[,1]
  df_plot$varname <- fct_recode(df_plot$varname, !!!drug_levels)

  p_cate <- ggplot(df_plot, aes(x = interaction(level, varname), y = tau, color = varname, group = 1)) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text.y = element_text(colour = "black"),
      strip.background = element_rect(colour = NA, fill = NA),
      legend.position = "none"
    ) +
    geom_point() +
    geom_errorbar(aes(ymin = tau - 1.96*se, 
                      ymax = tau + 1.96*se), width = .2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_hline(yintercept = ate) +
    annotate('ribbon', x = c(-Inf, Inf), ymin = ate_l, ymax = ate_u, 
             alpha = 0.4, fill = 'grey') +
    # facet_grid(varname ~ ., scales = "free_y") +
    coord_flip() +
    labs(y = "Treatment Effect", x = "") +
    scale_x_discrete(NULL, guide = "axis_nested") +
    # ylim(c(ylim_l, ylim_u)) +
    theme(text=element_text(size=9)) +
    scale_color_manual(values = c("#c52233","#f9844a","#ffa62b","#90be6d",
                                  "#43aa8b","#1282a2","#003459","#3d348b", "#a44200"))
  
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


plot_vim <- function(df_theta, estimator, vte, vte_l, vte_u){
  p <- ggplot(data=df_theta %>% 
                filter(method == estimator) %>% 
                arrange(importance) %>% 
                mutate(varname = factor(varname, levels = varname)),
              aes(x=varname, y=importance, ymin=ci_l, ymax=ci_u)) +
    geom_pointrange() + 
    # geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
    # geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    geom_hline(yintercept=0, linetype=2) +
    # geom_hline(yintercept = vte) +
    # annotate('ribbon', x = c(-Inf, Inf), ymin = vte_l, ymax = vte_u, 
    #          alpha = 0.4, fill = 'grey') +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("") + ylab("Importance") + # ylim(c(-350,350)) +
    theme_bw()  # use a white background
  return(p)
}


covar_levels = c(`Country` = "COUNTRY", `Age` = "AGE", `Sex` = "SEX", `Race` = "RACE", 
                 `Smoker status` = "SMOKER", `Diabetes duration` = "DIABDUR", 
                 `Antidiabetic therapy at baseline` = "ANTDBFL", `NYHA class (I - IV)` = "NYHACLAS", 
                 `Serum creatinine at baseline` = "CREATBL", `eGFR (Baseline) using MDRD (Calc)` = "EGFMDRBC", 
                 `Diabetic retinopathy severity` = "RETINSEV", `Myocardial infarction flag` = "MIFL", 
                 `Revascularization flag` = "REVASFL", `Carotid >50% stenosis on angiography` = "STENFL", 
                 `Coronary heart disease flag` = "CHDFL", `Ischaemic heart disease flag` = "IHDFL", 
                 `Chronic heart failure NYHA II-III flag` = "CHFFL", `Chronic kidney failure flag` = "KIDFL", 
                 `Microalbuminuria or proteinuria flag` = "MICFL", `Hypertension and LVH flag` = "HYPFL", 
                 `Left ventricular systolic and diastolic dysfunction flag` = "LVSDFL", 
                 `Ankle/Brachial index <0.9 flag` = "PADFL", `CV risk category` = "CVRISK", 
                 `HbA1c group (N)` = "HBA1CGRN", `Diabetes duration group (N)` = "DDURGRN", 
                 `Antihypertensive therapy` = "AHYPERFL", `Number of inclusion criterion passed` = "INCPASSN", 
                 `BMI at baseline` = "BMIBL", `Pulse at baseline` = "PULSEBL", 
                 `Systolic BP at baseline` = "SYSBPBL", `Diastolic BP at baseline` = "DIABPBL", 
                 `HbA1c at baseline (SI)` = "HBA1CBL", `HDL cholesterol at baseline (SI)` = "HDL1BL", 
                 `Calc. LDL cholesterol at baseline (SI)` = "LDL1BL", `Total cholesterol at baseline (SI)` = "CHOL1BL", 
                 `Triglycerides at baseline (SI)` = "TRIG1BL", `eGFR (Baseline) using MDRD (Calculated)` = "EGFMDRBC", 
                 `Peptic Ulcer and GERD Meds baseline flag` = "GERDBLFL", `Proton pump inhibitors flag` = "PPIFL", 
                 `H2 blockers flag` = "H2BLFL", `Stroke flag` = "STROKEFL",
                 `Statin use` = "statin_use", 
                 `Antihypertensives` = "antihypertensives",
                 `Beta blockers` = "betab",
                 `Remant cholesterol` ="RC",
                 `Mineralocorticoid receptor antagonists` = "minera", 
                 `ADP receptor inhibitors`= "adp",
                 `Vitamin K antagonists` = "vkantag",
                 `Ca antagonists` = "caantag", 
                 `Thiazide` = "thiazide",
                 `Loop diuretic` = "loopdiur")

drug_levels = c(`Statin use` = "statin", 
                `Anti-hypertensives` = "antihy",
                `Beta blockers` = "betab",
                `Mineralocorticoid receptor antagonists` = "minera", 
                `ADP receptor inhibitors`= "adp",
                `Vitamin K antagonists` = "vkanta",
                `Ca antagonists` = "caanta", 
                Thiazide = "thiazi",
                `Loop diuretic` = "loopdi")

list_ws <- list(
  "Country" = c("COUNTRY"),
  "Age" = c("AGE"), 
  "Sex" = c("SEX"), 
  "Race" = c("RACE"),
  "Smoker status" = c("SMOKER"),
  "Diabetes duration" = c("DIABDUR"),
  "Antidiabetic therapy at baseline" = c("ANTDBFL"),
  "NYHA class (I - IV)" = c("NYHACLAS"),
  "Kidney function" = c("CREATBL", "EGFMDRBC"),
  "Diabetic retinopathy severity" = c("RETINSEV"),
  "Previous CVD" = c("MIFL", "REVASFL", "STENFL", "CHDFL", "IHDFL", "HYPFL", "PADFL", "STROKEFL"),
  "Heart failure" = c("CHFFL", "LVSDFL", "CVRISK"),
  "Kidney disease" = c("KIDFL", "MICFL"),
  "BP medication" = c("AHYPERFL", "antihypertensives", "betab", "minera", "caantag", "thiazide", "loopdiur"),
  "Number of inclusion criterion passed" = c("INCPASSN"),
  "BMI at baseline" = c("BMIBL"),
  "Pulse at baseline" = c("PULSEBL"),
  "Blood pressure" = c("SYSBPBL", "DIABPBL"),
  "Lipids" = c("HDL1BL", "LDL1BL", "CHOL1BL", "TRIG1BL", "statin_use", "RC", "RCoverHDL"),
  "Peptic Ulcer and GERD Meds" = c("GERDBLFL", "PPIFL", "H2BLFL"),
  "HbA1c at Baseline (SI)" = c("HBA1CBL"),
  "ADP" = c("adp"),
  "Vitamin K antagonists" = c("vkantag"),
  "Insulin naive" = c("INSNVFL")
)