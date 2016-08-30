setwd("/scratch/zojer/projects/varcap_data2/coverage/dirR/")
dataset <- read.table("cov_final.txt", header=FALSE, sep="\t")
dataset
dataset1 <- dataset[c(1:3),c(1:5)]
# dataset1
library('ggplot2')
library('reshape2')
colnames(dataset1, do.NULL = FALSE, prefix = "col")
dataset1
# melted = melt(dataset1, id.vars=dataset1[,c(2)])
# melted
# colnames(melted) <- c("index","pos","cov")
# ggplot(data=melted, aes(x=pos, y=cov, group=variable)) + geom_line()
# write.table(dataset, "dataset.csv", sep="\t", row.names=FALSE, col.names=FALSE)


