
rm(list = ls())
source("R/data_process/data_helper.R")
library(survival)
library(ggsurvfit)

df <- get_data_tte(outcome = 'diab')


# KM plot
p_KM <- survfit2(Surv(Y, delta) ~ 1, data = df) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Survival Probability"
  ) + 
  add_confidence_interval() + 
  add_risktable()
p_KM


# survival plot (Liraglutide v.s. Placebo)
p_A <- survfit2(Surv(Y, delta) ~ A, data = df %>% 
                   mutate(A = ifelse(A == 1, 'Liraglutide', 'Placebo'))
                 ) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Survival Probability"
  ) + 
  add_confidence_interval() +
  add_risktable()
p_A


# survival plot (Liraglutide v.s. Placebo)
p_statin <- survfit2(Surv(Y, delta) ~ statin_use, data = df %>% 
                        mutate(statin_use = ifelse(statin_use == TRUE, 'Statin', 'Non-statin'))
                      ) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Survival Probability"
  ) + 
  add_confidence_interval() +
  add_risktable()
p_statin