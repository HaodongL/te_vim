library(here)
library(tidyverse)
library(dplyr)
library(mice)
library(sl3)
library(tmle3)
library(caret)

rm(list = ls())


### ------------  Part 0. Import Dataset  ------------ ###

  # df <- get_data(outcome, t, rm_baseIns=T)
  df <- read_csv(here("data/supp/df_analy1.csv"))
  

# impute missingness
nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

final <- process_missing(df, nodes)$data
no_variance<-nzv(final)
colnames(final)[no_variance]


head(final)

# library(DataExplorer)
# library(rmarkdown)
# 
# 
# final %>%
#   DataExplorer::create_report(
#     output_format = html_document(toc = F,  theme = "yeti"),
#     output_file = here::here("reports",paste("analysis_data_EDA", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), sep="-")),
#     report_title = "EDA Report"
#     #y = "CVS"
#   )



library(SmartEDA)

# similarly, with dplyr syntax: df %>% ExpReport(...)
ExpReport(
  final,
  #Target="cardio",
  label=NULL,
  op_file="EDA_report.html",
  op_dir=here::here("reports"))

