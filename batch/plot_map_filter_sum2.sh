#!/bin/bash

echo "subtot <- read.delim(\"filter_mapped_info_sort.txt\",header=TRUE)" >sum_var.R
# echo "str(subtot)" >>sum_var.R
echo "subtot[\"READS\"] <- subtot[\"FILCOV\"]/1000000 " >>sum_var.R
echo "par(las=2)" >>sum_var.R
echo "par(mar=c(14,4,4,3))" >>sum_var.R
echo "plot(subtot\$RL,type=\"p\",col=\"red\",main=\"Filter and mapping results\",ylab=\"abundance\",xlab=NULL,ylim=c(0,100),xaxt=\"n\",pch=0,cex.axis=0.6)" >>sum_var.R
# echo "grid()" >>sum_var.R
echo "abline(v=c(1,6,11,16,21,26,31,36,41,45),lty=3,col=\"gray50\")" >>sum_var.R
echo "axis(1,at=c(1:24),labels=rownames(subtot), pos=-5, lty=3, col=\"gray50\", las=2, tck=1, cex.axis=0.6)" >>sum_var.R
echo "points(subtot\$FILPCT,col=\"lightblue\",main=\"Filter and mapping results\",ylab=\"pct_filter\")" >>sum_var.R
echo "points(subtot\$READS,col=\"blue\",main=\"Filter and mapping results\",ylab=\"tot_filter\",pch=ifelse(subtot\$RL > 40,3,4),cex=0.7)" >>sum_var.R
echo "points(subtot\$MAPPCT,col=\"green\",main=\"Filter and mapping results\",ylab=\"pct_mapped\")" >>sum_var.R
Rscript sum_var.R
mv Rplots.pdf filter_map_info_sum_Rplots.pdf

