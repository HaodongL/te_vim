library(here)
library(tidyverse)
library(dplyr)
library(mice)
library(sl3)
library(tmle3)

rm(list = ls())
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/data_process/data_helper.R"))

### ------------  Part 0. Import Dataset  ------------ ###
if(here::here()=="C:/Users/andre/Documents/te_vim"){
  library(boxr)
  box_auth()
  
  
}else{
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

df <- process_missing(df, nodes)$data





