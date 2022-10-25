

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
  


