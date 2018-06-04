setwd('E://Rdata')
CNV1 <- read.table('CNV1',stringsAsFactors = FALSE)
RNASeq1 <- read.table('RNASeq1',stringsAsFactors = FALSE)
pvalue1 <- c('H vs N')
pvalue2 <- c('N vs L')
js <- c('js')
for (i in 2:nrow(CNV1)){
  print(i)
  high <- c()
  normal <- c()
  low <- c()

for (j in 2:ncol(CNV1)){

    if ((CNV1[i,j] == 1) | (CNV1[i,j] == 2)) {high <- c(high,RNASeq1[i,j]) }
      else if ((CNV1[i,j] == -1) | (CNV1[i,j] == -2)){low <- c(low,RNASeq1[i,j])}
      else {normal <- c(normal,RNASeq1[i,j])}
  }

{
if ((length(high)>1) & (length(normal)>1) & (length(low)>1)){
  ttest1 <- t.test(as.numeric(high),as.numeric(normal),alternative = "greater")
  ttest2 <- t.test(as.numeric(normal),as.numeric(low),alternative = "greater")
  pvalue1 <- c(pvalue1,ttest1$p.value)
  pvalue2 <- c(pvalue2,ttest2$p.value)
  js1 <- 0
  
  
  if  (is.nan(ttest1$p.value)|(is.nan(ttest2$p.value))) {js1 <- -2}
  if  (!(is.nan(ttest1$p.value)))
         {if (ttest1$p.value < 0.05) {js1 = js1 + 1}}
  if  (!(is.nan(ttest2$p.value)))
         {if (ttest2$p.value < 0.05) {js1 = js1 + 1}}
  }
  else {js1 <- -1
        pvalue1 <- c(pvalue1,-1)
        pvalue2 <- c(pvalue2,-1)}
  js <- c(js,js1)

}
}
CNV1 <- cbind(CNV1,pvalue1,pvalue2,js)
RNASeq1 <- cbind(RNASeq1,pvalue1,pvalue2,js)
write.table(CNV1,"CNV-final")
write.table(RNASeq1,"RNAseq-final")
