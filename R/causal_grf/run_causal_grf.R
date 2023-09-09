rm(list = ls())
source(paste0(here::here(),"/R/0_config.R"))

source(paste0(repo_path, "R/data_process/data_helper.R"))
source(paste0(repo_path, "R/simu/simu_dgd.R")) 
source(paste0(repo_path, "R/simu/simu_config.R"))
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/analysis/analy_helper.R"))

library(causalTree)
library(grf)
library(glmnet)
library(splines)
library(MASS)
library(lmtest)
library(sandwich)
library(ggplot2)
library(stringr)
library(tidyverse)
library(haven)
library(flextable)

#grf vignette
#https://grf-labs.github.io/grf/articles/grf_guide.html

### ------------  Part 1. import data  ------------ ###
df <- get_data(outcome="diab2", t=24, rm_baseIns=T)
nodes <- list(W = setdiff(names(df), c("Y", "A")), A = 'A', Y = 'Y')
df <- process_missing(df, nodes)$data

#Note! Drop imputed data or using in the causal forest estimation

ws = c('statin_use')


ws = c("AGE", "SEX", "RACE", "COUNTRY", "SMOKER", "NYHACLAS",
       "DIABDUR", "ANTDBFL", "AHYPERFL", "INCPASSN", 
       "BMIBL", "PULSEBL", "SYSBPBL", "DIABPBL", "HBA1CBL", "HDL1BL", "LDL1BL",
       "CHOL1BL", "RC", "RCoverHDL","TRIG1BL", "CREATBL", "EGFMDRBC", 
       "RETINSEV", "GERDBLFL", "PPIFL", "H2BLFL", 
       "MIFL", "STROKEFL", "REVASFL", "STENFL", "CHDFL", "IHDFL", "CHFFL",
       "KIDFL", "MICFL", "HYPFL", "LVSDFL", "PADFL", "CVRISK", "HBA1CGRN", "DDURGRN", 
       "statin_use", "antihypertensives", "betab", "minera", "adp",
       "vkantag", "caantag", "thiazide", "loopdiur")
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")

df <- as.data.frame(df)
#X <- df %>% subset(., select=cm_names)
X <- df %>% subset(., select=ws)
X <- washb::design_matrix(X)
colnames(X) <- gsub("TRUE","",colnames(X))
W <- df$A 
Y <- df$Y

set.seed(12345)
res_cf <- causal_forest(
  X = X,                       #covariates
  Y = Y,                       #outcomes  
  W = W                       #treatment
)

#average treatment effect
ate <- as.data.frame(t(average_treatment_effect(res_cf)))
ate

varimp = variable_importance(res_cf)
ranked.vars <- order(varimp, decreasing = TRUE)
ranked.vars

#Best linear projections
#grf_res <- as.matrix(best_linear_projection(res_cf, A=X[ranked.vars[1:10]])) #top 10
grf_blp <- as.matrix(best_linear_projection(res_cf, A=X[colnames(X) %in% cm_names])) #drugs to examine
grf_res <- as.data.frame(grf_blp[2:nrow(grf_blp),c(1,2,4)])
grf_res$Estimate <- grf_res$Estimate+grf_blp[1,1] # add intercept to get CATE. NOTE! Need to get the linear combo of SE
grf_res$cate=rownames(grf_res)
grf_res=bind_rows(ate %>% rename(Estimate =estimate,`Std. Error`=std.err) %>% mutate(cate="overall"),grf_res)
grf_res <- grf_res %>% mutate(cate=factor(cate, levels=unique(cate)))
grf_res$cate <- fct_recode(grf_res$cate, !!!subgroup_levels)

rownames(grf_res)=NULL
grf_res$ci.lb <- grf_res$Estimate - 1.96*grf_res$`Std. Error`
grf_res$ci.ub <- grf_res$Estimate + 1.96*grf_res$`Std. Error`

grf_plot <- ggplot(grf_res, aes(x=cate, y=Estimate)) + geom_point() + coord_flip() +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub)) + theme_bw() + geom_hline(yintercept=0) +
  xlab("Drug usage") + ylab("CATE")
grf_plot

ggsave(paste0(here(),"/tnp/plot/p_grf_cate.png"), grf_plot, width = 5, height = 5)
saveRDS(grf_res, paste0(here(),"/analy_res/res_grf_cate_diab2.rds"))





#NOTE! I don't think this is right... overall is way too different. I think need to covert the BLP predictions based on the regression equation
#look at Caitlin's manuscript
#need to include the non-usage CATES using the intercept
#Convert variable names to to plot labels
#Convert the

# estimand=test_calibration(res_cf)
# res_blp <- best_linear_projection(res_cf)
# plot(res_blp)
# 
# tau.hat <- predict(res_cf)$predictions
# 
# plot(tree <- get_tree(res_cf, 50))


# #Tree based policy learning
# 
# samples.by.school <- split(seq_along(school), school)
# num.schools <- length(samples.by.school)
# train <- unlist(samples.by.school[sample(1:num.schools, num.schools / 2)])
# 
# train.forest <- causal_forest(X[train, ], Y[train], W[train], W.hat = 0.5, clusters = school[train])
# tau.hat.eval <- predict(train.forest, X[-train, ])$predictions
# 
# eval.forest <- causal_forest(X[-train, ], Y[-train], W[-train], W.hat = 0.5, clusters = school[-train])
# 
# rate.cate <- rank_average_treatment_effect(eval.forest, tau.hat.eval)
# plot(rate.cate, main = "TOC: By decreasing estimated CATE")
# 
# library(policytree)
# # Use train/evaluation school split from above, but use non-missing units for policy_tree
# eval <- (1:nrow(X))[-train]
# not.missing <- which(complete.cases(X))
# train <- train[which(train %in% not.missing)]
# eval <- eval[which(eval %in% not.missing)]
# 
# # Compute doubly robust scores
# dr.scores <- get_scores(cf)
# # Use as the ATE as a "cost" of program treatment to find something non-trivial
# cost <- ate[["estimate"]]
# dr.rewards <- cbind(control=-dr.scores, treat=dr.scores - cost)
# 
# # Fit depth 2 tree on training subset
# tree <- policy_tree(X[train, ], dr.rewards[train, ], min.node.size = 100)
# plot(tree, leaf.labels = c("dont treat", "treat"))

