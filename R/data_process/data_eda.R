library(dplyr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggcorrplot)
library(RColorBrewer)
library(tmle3)

# read baseline variables
df_w <- read_csv("data/supp/df_w.csv")
df_w <- df_w %>% 
  mutate(across(where(is.character), ~ as.factor(.))) %>% 
  mutate(across(where(is.logical), ~ as.factor(.))) 
df_w <- df_w %>% select(-c('A', 'USUBJID'))
# sapply(df_w, class)

# remove INSNVFL == TRUE
df_w <- df_w %>% filter(INSNVFL == FALSE) %>% select(-INSNVFL)

df_w_numeric <- df_w %>% select(where(is.numeric))
df_w_factor <- df_w %>% select(where(is.factor))

# correlation (numeric W)
color_corr <- c("#023e8a", "white", "#ae2012")
model.matrix(~0+., data=df_w_numeric) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=2, colors = color_corr)

# correlation (concomitant med (factor))
cm_names <- c("statin_use", "antihypertensives", "betab", "minera", "adp",
              "vkantag", "caantag", "thiazide", "loopdiur")
df_cm <- df_w_factor %>% select(all_of(cm_names))

model.matrix(~., data=df_cm)[,-1] %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=2, colors = color_corr)

# remove highly correlated cols
var_corr <- c("EGFRMDRC", "EGFREPB", "EGFMDRBC", "EGFRMDR")
df_w <- df_w %>% select(-all_of(var_corr))

# deal with missingness
# nodes <- list(W = setdiff(names(df_w), 'SEX'),
#               Y = "SEX")
# df_w_impute <- process_missing(df_w, nodes)$data

# drop missing rows and constant cols
df_w_model <- model.matrix(~., data=df_w)[,-1]
df_w_model <- df_w_model[, -c(16,19)]
# for (i in c(1:ncol(df_w_model))){
#   if (var(df_w_model[,i]) == 0){
#     print(i)
#   }
# }


# saveRDS(df_w_model, file = "data/supp/df_w_model.RDS")

# PCA 
pca_c <- prcomp(x = df_w_model, center = T, scale = T)

ling_scores <- as.data.frame(pca_c$x)

a = 2
b = 3
p <-
  ggplot(ling_scores) + 
  geom_point(aes(x = ling_scores[, paste0("PC",a)], y = ling_scores[, paste0("PC",b)]),
             size = 0.8, alpha=0.7) +
  labs(x=paste0("PC",a), y=paste0("PC",b))
p






