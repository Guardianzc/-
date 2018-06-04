#1.SEARTH AND CAPTURE #抓取sample中的样本
setwd('E://Rdata')
sample <- read.table('sample.txt',stringsAsFactors = FALSE)    #不加这个参数后面提取一行出来的数据会带level
CNV <- read.table('Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes')
RNASeq <- read.table('HiSeqV2')

CNV_after <- data.frame(matrix(NA,24777,2))
RNASeq_after <- data.frame(matrix(NA,20531,2))  #这里是可以优化的一个部分，因为一行的dataframe会变成vector所以在初始化的时候先加了两个空行

CNV_after <- cbind(CNV_after,CNV[,1])
RNASeq_after <- cbind(RNASeq_after,RNASeq[,1])

for (i in 2:nrow(sample)){
  if ((sample[i,3] == 'YES') & (sample[i,4] == 'YES')){    #这里是CNA行和SEQ行
    sampleID = sample[i,1]
    set1 <- which(CNV[1,]==sampleID)
    if (length(set1)!=0){
    CNV_after <- cbind(CNV_after,CNV[,set1])
    }
    set2 <- which(RNASeq[1,]==sampleID)
    if (length(set2)!=0){
    RNASeq_after <- cbind(RNASeq_after,RNASeq[,set2])
    }
  }
}
CNV_after <- CNV_after[,-1]
CNV_after <- CNV_after[,-1]
RNASeq_after <- RNASeq_after[,-1]                    
RNASeq_after <- RNASeq_after[,-1]                        #这里把一开始的空行删掉
write.table(CNV_after,'CNV_after',row.names = FALSE,col.names = FALSE,sep = '\t')
write.table(RNASeq_after,'RNASeq_after',row.names = FALSE,col.names = FALSE,sep = '\t')

#2.col #调整两组数据具有相同的行
CNV_after <- read.table('CNV_after',stringsAsFactors = FALSE)
RNASeq_after <- read.table('RNASeq_after',stringsAsFactors = FALSE)
CNV1 <- data.frame(matrix(NA,2,168))
RNASeq1 <- data.frame(matrix(NA,2,168))   #老样子这里又多了两个空行，项南师兄说可以用drop= FALSE 把它去掉
cnames=paste("V",1:168,sep="")
colnames(CNV1)=cnames
colnames(RNASeq1)=cnames
CNV1 <- rbind(CNV1,CNV_after[1,])
RNASeq1 <- rbind(RNASeq1,RNASeq_after[1,])

for (i in 2:nrow(RNASeq_after)){
  print(i)
  GeneID = RNASeq_after[i,1]
  set1 <- which(CNV_after[,1]==GeneID)
  if (length(set1)!=0){
    CNV1 <- rbind(CNV1,CNV_after[set1,])
    RNASeq1 <- rbind(RNASeq1,RNASeq_after[i,])
    }
  
}
CNV1 <- CNV1[-1,]
CNV1 <- CNV1[-1,]
RNASeq1 <- RNASeq1[-1,]
RNASeq1 <- RNASeq1[-1,]                                                         #去掉空行
write.table(CNV1,'CNV1',row.names = FALSE,col.names = FALSE,sep = '\t')
write.table(RNASeq1,'RNASeq1',row.names = FALSE,col.names = FALSE,sep = '\t')

#3.Pearson & p-value  #按照pearson算相关性和pvalue
setwd('E://Rdata')
CNV2 <- read.table('CNV1')
RNASeq2 <- read.table('RNASeq1')
CNV3 <- t(CNV2)
CNV3 <- CNV3[-1,]

RNASeq3 <- t(RNASeq2)
RNASeq3 <- RNASeq3[-1,]
                                                            #这里把矩阵转置取行，因为算pearson时要用vector，取列好算
                                                            #这里RNASeq为0的值算出来相关是NA，注意这一点
cvalue=paste("cv",1:nrow(CNV2),sep="")
pvalue=paste("pv",1:nrow(CNV2),sep="")
for (i in 2:nrow(CNV2)){
  print(i)
  a <- as.double(CNV3[,i])
  b <- as.double(RNASeq3[,i])
  cor1 <- cor.test(a,b,method="pearson")
  cvalue[i] <- cor1$estimate
  pvalue[i] <- cor1$p.value
}