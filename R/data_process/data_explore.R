library(here)
library(tidyverse)
library(dplyr)
library(mice)
library(sl3)
library(tmle3)
library(caret)

rm(list = ls())


### ------------  Part 0. Import Dataset  ------------ ###
if(here::here()=="C:/Users/andre/Documents/jici/te_vim"){
  
  repo_path = "C:/Users/andre/Documents/jici/te_vim/"
  source(paste0(repo_path, "R/data_process/data_helper.R"))
  
  #cleaned data
  df <- read_csv(paste0(here(),"/data/supp/df_analy1.csv"))
  dim(df)
  colnames(df)
  
  head(df)  
  
  
}else{
  
  repo_path = "~/Repo/te_vim/"
  source(paste0(repo_path, "R/data_process/data_helper.R"))
  # read in analysis data
  ### ------------  Part 1. import data  ------------ ###
  # outcome = 'diab'; t = 24
  # df <- get_data(outcome, t, rm_baseIns=T)
  df <- read_csv("~/Repo/te_vim/data/supp/df_analy1.csv")
}
  

# impute missingness
nodes <- list(W = setdiff(names(df), c("Y", "A")),
              A = 'A',
              Y = 'Y')

final <- process_missing(df, nodes)$data
no_variance<-nzv(final)
colnames(final)[no_variance]


head(final)

library(DataExplorer)


final %>%
  create_report(
    output_format = html_document(toc = F,  theme = "yeti"),
    output_file = paste("analysis_data_EDA", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), sep=" - "),
    report_title = "EDA Report"
    #y = "CVS"
  )



library(SmartEDA)

# similarly, with dplyr syntax: df %>% ExpReport(...)
ExpReport(
  final,
  #Target="cardio",
  label=NULL,
  op_file="Report.html",
  op_dir=getwd())

