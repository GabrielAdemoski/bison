
source("src/correlation.R")

individual.corr(catagories=c("Seasons", "Years"), levels=c("l2", "l6"), xlsx.prefix='taxa-level')
delta.corr(catagories=c("Seasons", "Years"), levels=c("l2", "l6"), xlsx.prefix='taxa-level')
