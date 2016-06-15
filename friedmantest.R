/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(reshape2)
#options(scipen=999)
tt=read.table(args[1], header=T)
melt.tt=melt(tt)
colnames(melt.tt)=c("ID","TF","value")
#Friedman
fried.result=friedman.test(value~TF|ID,data=melt.tt)
print(fried.result)

#Friedman namenyi test
library(PMCMR)
friedman.result2=posthoc.friedman.nemenyi.test(value~TF|ID,data=melt.tt)
friedman.result2
print(fried.result2)
