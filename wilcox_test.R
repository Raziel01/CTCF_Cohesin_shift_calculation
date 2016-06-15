#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(reshape2)
#options(scipen=999)
tt=read.table(args[1], header=T)
melt.tt=melt(tt)
colnames(melt.tt)=c("ID","TF","value")

x=tt$CTCF                                                                    #Renaming  the variables
y=tt$Factor1
wilcox.test(x,y, paired=TRUE)
