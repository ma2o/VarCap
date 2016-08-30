# samtools view -f66 file.bam|cut -f 9|sed 's/^-//' > InsertSizeMetrics.txt
# grep -v '[0-9]\{5,\}' InsertSizeMetrics_samtools.txt >InsertSizeMetrics_samtools_1k.txt
mydata = scan("InsertSizeMetrics_samtools_1k.txt")
summary(mydata)
# hist(mydata, breaks=20, main="Breaks=20")
# hist(mydata, main="Insert_Size _Metrics", xlab="IS length", ylab="Counts", breaks=c(0,1000,by=25),xlim = c(0,1000))
