library(speff2trial)
library(dplyr)
library(sl3)
library(tmle3)
library(foreach)
library(doParallel)


rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))

### ------------  Part 1. Data ------------ ###
data(ACTG175)
df <- ACTG175

df <- df %>% 
  filter(arms == 1 | arms == 3) %>% 
  select(c("arms", "cd420", "age", "wtkg", "karnof", "cd40", "cd80",
           "homo", "gender", "race", "symptom", "drugs", "hemo", "str2")) %>% 
  rename(A = "arms",
         Y = "cd420") %>% 
  mutate(A = as.numeric(A == 3))

### ------------  Part 2. Estimation ------------ ###
# vim loop over all covariates
set.seed(1)
ws = c("age", "wtkg", "karnof", "cd40", "cd80", "homo", 
       "gender", "race", "symptom", "drugs", "hemo", "str2")

n_ws <- length(ws)


registerDoParallel(3)
df_vim <- foreach(i = 1:n_ws, .combine = 'rbind') %dopar% {
  print(paste0("Covar name: ", ws[i]))
  res <- run_VIM(df = df,
                 sl_Q = sl_Q, 
                 sl_g = sl_g,
                 sl_x = sl_x,
                 ws = ws[i], 
                 cv = F,
                 dr = F,
                 max_it = 1e4)
  res_ee <- res$resEE
  res_tmle <- res$resTMLE
  res_ss <- res$resSS
  
  df_vim_i <- data.frame("varname" = c(ws[i], ws[i], ws[i]), 
                         "importance" = c(res_ee$coef, res_tmle$coef, res_ss$coef), 
                         "ci_l" = c(res_ee$ci_l, res_tmle$ci_l, res_ss$ci_l), 
                         "ci_u" = c(res_ee$ci_u, res_tmle$ci_u, res_ss$ci_u), 
                         "method" = c('EE', 'TMLE', 'SS'))
  return(df_vim_i)
}

saveRDS(df_vim, file = "~/Repo/te_vim/demo/demo_res_nocvt.RDS")