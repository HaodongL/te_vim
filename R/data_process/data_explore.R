

#ADCM 
library(boxr)
box_auth()
#concommitant medication usage
d <- box_read("946034255301")
head(d)
unique(d$CMTRT)


#hypoglycemia data
adhypo <- box_read("946034248101")
head(adhypo)
  


#cleaned data
final <- box_read("1048393668373")
dim(final)
colnames(final)

head(final)  

library(caret)
no_variance<-nzv(final)
colnames(final)[no_variance]


final#

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

