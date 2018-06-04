#3.Merge.R
setwd("d:/Rdata/CNA/Delly/ARD")
#merge 1
files <- list.files(pattern = "\\.dup.pass.precise.vcf.INFO")
num <- length(files) #the number of files
for (i in 1:num) {
  CNA<- read.table(files[i],sep="\t",header=F)
  CNA<-CNA[which((CNA[,3]-CNA[,2])<1000000),]
  dat<-data.frame()
  n <- nrow(CNA)
  j <- 1
  while(j <= n){
    k <- 0
    new_dat <- data.frame()
    while(k < 11){
      k <- k+1;if(is.na(CNA[j+k,2])==TRUE|CNA[j+k,1]!=CNA[j,1]|as.numeric(CNA[j+k,2])< as.numeric(CNA[j,2])|as.numeric(CNA[j+k,2]) > (as.numeric(CNA[j,3])+10000)) break;
    }
    new_dat <- data.frame(CNA[j,1],CNA[j,2],max(CNA[j:(j+k-1),3])) 
    dat <- rbind(dat,new_dat)
    j <- j+k
  }
  colnames(dat) <- c("Chr","Start","End")
  write.table(dat,paste(files[i],"merge1",sep="."),sep="\t",col.names=TRUE,row.names=FALSE)
}
#merge 2
files <- list.files(pattern = "\\.dup.pass.precise.vcf.INFO.merge1")
num <- length(files) #the number of files
for (i in 1:num) {
  #i <- 1
  CNA<- read.table(files[i],sep="\t",header=T)
  CNA<-CNA[which((CNA[,3]-CNA[,2])<1000000),]
  dat<-data.frame()
  n <- nrow(CNA)
  j <- 1
  while(j <= n){
    k <- 0
    new_dat <- data.frame()
    while(k < 11){
      k <- k+1;if(is.na(CNA[j+k,2])==TRUE|CNA[j+k,1]!=CNA[j,1]|as.numeric(CNA[j+k,2])< as.numeric(CNA[j,2])|as.numeric(CNA[j+k,2]) > (as.numeric(CNA[j,3])+10000)) break;
    }
    new_dat <- data.frame(CNA[j,1],CNA[j,2],max(CNA[j:(j+k-1),3])) 
    dat <- rbind(dat,new_dat)
    j <- j+k
  }
  colnames(dat) <- c("Chr","Start","End")
  write.table(dat,paste(files[i],"merge2",sep="."),sep="\t",col.names=TRUE,row.names=FALSE)
}
#merge 3
files <- list.files(pattern = "\\.dup.pass.precise.vcf.INFO.merge1.merge2")
num <- length(files) #the number of files
for (i in 1:num) {
  #i <- 1
  CNA<- read.table(files[i],sep="\t",header=T)
  CNA<-CNA[which((CNA[,3]-CNA[,2])<1000000),]
  dat<-data.frame()
  n <- nrow(CNA)
  j <- 1
  while(j <= n){
    k <- 0
    new_dat <- data.frame()
    while(k < 11){
      k <- k+1;if(is.na(CNA[j+k,2])==TRUE|CNA[j+k,1]!=CNA[j,1]|as.numeric(CNA[j+k,2])< as.numeric(CNA[j,2])|as.numeric(CNA[j+k,2]) > (as.numeric(CNA[j,3])+10000)) break;
    }
    new_dat <- data.frame(CNA[j,1],CNA[j,2],max(CNA[j:(j+k-1),3])) 
    dat <- rbind(dat,new_dat)
    j <- j+k
  }
  colnames(dat) <- c("Chr","Start","End")
  write.table(dat,paste(files[i],"merge3",sep="."),sep="\t",col.names=TRUE,row.names=FALSE)
}
#merge 4
files <- list.files(pattern = "\\.dup.pass.precise.vcf.INFO.merge1.merge2.merge3")
num <- length(files) #the number of files
for (i in 1:num) {
  #i <- 1
  CNA<- read.table(files[i],sep="\t",header=T)
  CNA<-CNA[which((CNA[,3]-CNA[,2])<1000000),]
  dat<-data.frame()
  n <- nrow(CNA)
  j <- 1
  while(j <= n){
    k <- 0
    new_dat <- data.frame()
    while(k < 11){
      k <- k+1;if(is.na(CNA[j+k,2])==TRUE|CNA[j+k,1]!=CNA[j,1]|as.numeric(CNA[j+k,2])< as.numeric(CNA[j,2])|as.numeric(CNA[j+k,2]) > (as.numeric(CNA[j,3])+10000)) break;
    }
    new_dat <- data.frame(CNA[j,1],CNA[j,2],max(CNA[j:(j+k-1),3])) 
    dat <- rbind(dat,new_dat)
    j <- j+k
  }
  colnames(dat) <- c("Chr","Start","End")
  write.table(dat,paste(files[i],"merge4",sep="."),sep="\t",col.names=TRUE,row.names=FALSE)
}

#4.Classify.R
setwd("d:/Rdata/CNA/Delly/ARD")
library(data.table)
##segment level 
#L1 0-0.1Kb
#L2 0.1-1kb
#L3 1-10Kb
#L4 10-100Kb
#L5 >100Kb
#Before Merge
files <- list.files(pattern = "dup.pass.precise.vcf.INFO")
n <- length(files)
for(i in 1:n) {
dat<-fread(files[i],sep="\t",header=TRUE,data.table=FALSE)
name <- sub(pattern="dup.pass.precise.vcf.INFO",replacement="dup.pass.precise.vcf.INFO.LEVEL",files[i])
sampleID <- sub(pattern=".dup.pass.precise.vcf.INFO",replacement="",files[i])
colnames(dat) <- c("Chr","Start","End")
Level <- apply(dat, 1, function(x){
len <- as.numeric(x[3])- as.numeric(x[2]) 
level <- ifelse(len <100,'L1',(ifelse(len<1000,'L2',(ifelse(len<10000,'L3',(ifelse(len<100000,'L4','L5')))))))
})
dat$LEVEL <- Level
dat$SampleID <- sampleID
write.table(dat,name,sep="\t",col.names=TRUE,row.names=FALSE)
}
#After Merge
files <- list.files(pattern = "dup.pass.precise.vcf.INFO.merge1.merge2.merge3.merge4")
n <- length(files)
for(i in 1:n) {
dat<-fread(files[i],sep="\t",header=TRUE,data.table=FALSE)
name <- sub(pattern="dup.pass.precise.vcf.INFO.merge1.merge2.merge3.merge4",replacement="dup.pass.precise.vcf.INFO.merge1.merge2.merge3.merge4.LEVEL",files[i])
sampleID <- sub(pattern=".dup.pass.precise.vcf.INFO.merge1.merge2.merge3.merge4",replacement="",files[i])
colnames(dat) <- c("Chr","Start","End")
Level <- apply(dat, 1, function(x){
len <- as.numeric(x[3])- as.numeric(x[2]) 
level <- ifelse(len <100,'L1',(ifelse(len<1000,'L2',(ifelse(len<10000,'L3',(ifelse(len<100000,'L4','L5')))))))
})
dat$LEVEL <- Level
dat$SampleID <- sampleID
write.table(dat,name,sep="\t",col.names=TRUE,row.names=FALSE)
}

#5.MergeChange.R
##Merge before and after the change.
#Merge before 
files <- list.files(pattern = ".dup.pass.precise.vcf.INFO.LEVEL$")
num <- length(files) #the number of files
before <- data.frame() #Number of deltetion segment that merge before
for(i in 1:num){
	merge_before<- read.table(files[i],sep="\t",header=T) #Merge before
	before_EachLevel_Number <- t(cbind(table(merge_before$LEVEL)))
	before_SegmentNumber <- nrow(merge_before)
	SampleID <- sub(pattern=".dup.pass.precise.vcf.INFO.LEVEL",replacement="",files[i])
	Process <- c("Before")
	before_message <- data.frame(SampleID,Process,before_EachLevel_Number,before_SegmentNumber)
	before <- rbind(before,before_message)
}
colnames(before) <- c("SampleID","Process","L2","L3","L4","L5","Total")

#Merge after
files <- list.files(pattern = ".dup.pass.precise.vcf.INFO.merge1.merge2.merge3.merge4.LEVEL$")
num <- length(files) #the number of files
after <- data.frame() #Number of deltetion segment that merge after
for(i in 1:num){
	merge_after<- read.table(files[i],sep="\t",header=T) #Merge after
	after_EachLevel_Number <- t(cbind(table(merge_after$LEVEL)))
	after_SegmentNumber <- nrow(merge_after)
	SampleID <- sub(pattern=".dup.pass.precise.vcf.INFO.merge1.merge2.merge3.merge4.LEVEL",replacement="",files[i])
	Process <- c("After")
	after_message <- data.frame(SampleID,Process,after_EachLevel_Number,after_SegmentNumber)
	after <- rbind(after,after_message)
}
colnames(after) <- c("SampleID","Process","L2","L3","L4","L5","Total")
merge_change <- rbind(before,after)
write.table(merge_change,"ARD_merge_change.txt",sep="\t",col.names=TRUE,row.names=FALSE)

#The change of segment number
#boxplot
library("ggplot2")
library("data.table")
merge_change <- fread("ARD_merge_change.txt",header=T,sep="\t",data.table=FALSE)
merge_change$Process <- ordered(merge_change$Process, levels = c("Before","After"))
p <- ggplot(merge_change,aes(Process,Total))+ geom_boxplot(aes(fill=Process),show.legend = FALSE)
p <- p+theme(panel.grid=element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1.3))+theme_classic()
pdf("1.ARD_MergeChange_SegmentsNumber.pdf")
p + xlab("Merge")+ylab("Number of duplication segment")+labs(title ="")+
		theme(axis.title.x = element_text(face='bold',size=15,hjust=0.5),
        axis.title.y = element_text(face='bold',size=15,vjust=1),
        axis.text.x = element_text(face='bold',size=15,color='black'),
        axis.text.y = element_text(face='bold',size=15,color='black'),
        legend.title =  element_blank(),
        legend.text = element_text(face="bold", size=15),
        plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5)) + 
		scale_y_continuous(limits=c(100,400),breaks=c(100,200,300,400))
dev.off()
#boxplot
library(dplyr) 
library(tidyr)
merge_change_Level <- merge_change[,-7]
dat <- merge_change_Level %>% gather(Level,SegmentNumber,-1,-2)
p <- ggplot(dat, aes(x=Level, y=SegmentNumber, fill=Process))+geom_boxplot()
p <- p+theme(panel.grid=element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1.3))+theme_classic()
pdf("2.ARD_MergeChange_Level.pdf")
p + xlab("Segment Level")+ylab("Number of duplication segment")+labs(title ="")+
  theme(axis.title.x = element_text(face='bold',size=15,hjust=0.5),
        axis.title.y = element_text(face='bold',size=15,vjust=1),
        axis.text.x = element_text(face='bold',size=15,color='black'),
        axis.text.y = element_text(face='bold',size=15,color='black'),
        legend.title =  element_blank(),
        legend.text = element_text(face="bold", size=15),
        plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5))+theme(legend.position = c(0.8, 0.8))+ 
		scale_y_continuous(limits=c(100,200),breaks=c(100,150,200))
dev.off()

#6.Concordance.R
setwd("d:/Rdata/CNA/Delly/ARD")
library(data.table)
library(ggplot2)
library(dplyr) 
library(tidyr)
files <- fread("sample.txt",header=FALSE,sep="\t",data.table=FALSE)
#The Concordance of Replicate1 and Replicate2 
sample <- data.frame(rep(1:3,each=3),rep(1:3))
colnames(sample) <- c("Sample1","Sample2")
for(sam in 1:nrow(sample)){
  a <- sample[sam,1]
  b <- sample[sam,2]
  Concordance23 <- data.frame()
  if (a!=b){
    for (i in 1:nrow(files)) {
      #i <- 1
      data1 <- fread(paste0(files[i,],paste0("_",a,"."),"dup.pass.precise.vcf.INFO.merge1.merge2.merge3.merge4.LEVEL"),sep="\t",header=TRUE,data.table=FALSE)
      data2 <- fread(paste0(files[i,],paste0("_",b,"."),"dup.pass.precise.vcf.INFO.merge1.merge2.merge3.merge4.LEVEL"),sep="\t",header=TRUE,data.table=FALSE)
      ##LEVEL
      L <- 5
      ##OneSample Concordance
      One_Con <- data.frame()
      for(j in 2:L){
        #j <- 4
        dat1_LEVEL <- data1[data1$LEVEL==paste0("L",j),]
        dat2_LEVEL <- data2[data2$LEVEL==paste0("L",j),]
        ##dat1_dat2_match
        dat1_dat2_match<-data.frame() 
        for(k in unique(dat1_LEVEL$Chr)){
          dat1<- dat1_LEVEL[dat1_LEVEL$Chr==k,]
          dat2<- dat2_LEVEL[dat2_LEVEL$Chr==k,]
          if (nrow(dat2) > 0){
            ##INDEX
            dat2$INDEX <- c(1:nrow(dat2))
            #match1
            #x and y are the two endpoints of the dat1.
            match1 <- apply(dat1,1,function(Dat1,Dat2){
              x <- as.numeric(Dat1[2])
              y <- as.numeric(Dat1[3])
              Dat2_INDEX<-union(intersect(Dat2$INDEX[Dat2$End>=y],Dat2$INDEX[Dat2$Start<=0.5*(x+y)]),intersect(Dat2$INDEX[Dat2$End>=y],Dat2$INDEX[Dat2$Start<=(2*y-as.numeric(Dat2$End))]))
              n <- length(Dat2_INDEX)
              return(n)
            },dat2
            )
            #match2
            match2 <- apply(dat1,1,function(Dat1,Dat2){
              x <- as.numeric(Dat1[2])
              y <- as.numeric(Dat1[3])
              Dat2_INDEX<-union(intersect(Dat2$INDEX[Dat2$Start<x],Dat2$INDEX[Dat2$End>=0.5*(x+y)]),intersect(Dat2$INDEX[Dat2$Start<x],Dat2$INDEX[Dat2$End>=(2*x-as.numeric(Dat2$Start))]))
              n <- length(Dat2_INDEX)
              return(n)
            },dat2
            )
            #match3
            match3 <- apply(dat1,1,function(Dat1,Dat2){
              x <- as.numeric(Dat1[2])
              y <- as.numeric(Dat1[3])
              Dat2_INDEX<-intersect(Dat2$INDEX[Dat2$Start>=x],Dat2$INDEX[Dat2$End<y])
              n <- length(Dat2_INDEX)
              return(n)
            },dat2
            )
            dat1_dat2_new <- data.frame(match1<-c(match1),match2<-c(match2),match3<-c(match3))
            colnames(dat1_dat2_new) <- c("match1","match2","match3")
            dat1_dat2_match <- rbind(dat1_dat2_match,dat1_dat2_new)
          }
        }
        dat1_dat2_match$Sum <- rowSums(dat1_dat2_match[1:3])
        len <- nrow(dat1_dat2_match[dat1_dat2_match$Sum>0,])
        #dat1_dat2_concordance <- data.frame(sum(dat1_dat2_match$Sum)/((nrow(dat1_LEVEL)+nrow(dat2_LEVEL)+sum(dat1_dat2_match$Sum)-len)/2))
        #dat1_dat2_concordance <- data.frame(sum(dat1_dat2_match$Sum)/((nrow(dat1_LEVEL)+nrow(dat2_LEVEL))/2))
        dat1_dat2_concordance <- data.frame(len/((nrow(dat1_LEVEL)+nrow(dat2_LEVEL))/2))
        One_Con <- rbind(One_Con,dat1_dat2_concordance)
      }
      colnames(One_Con) <- c(files[i,])
      Concordance23 <- rbind(Concordance23,t(One_Con))
    }
    write.table(Concordance23,paste0(a,"Vs",b,".txt"),sep="\t",col.names=TRUE,row.names=TRUE)		
  }
}		

#7.Concordance_Boxplot.R  
files <- list.files(pattern="Vs")
n <- length(files)
dat <- data.frame()
for (i in 1:n){
  data <- read.table(files[i],header=TRUE,sep="\t",check.names=FALSE)
  dat <- rbind(dat,data)
}
colnames(dat) <- c("0.1~1","1~10","10~100",">100")
dat_new <- dat %>% gather(Level,Concordance)
dat_new$Level <- ordered(dat_new$Level, levels = c("0.1~1","1~10","10~100",">100"))
p <- ggplot(dat_new,aes(Level,Concordance,fill=Level))+geom_boxplot(show.legend = FALSE)
p <- p + geom_point(show.legend=FALSE,size=1.3)
p <- p +scale_fill_brewer(palette="Blues") +theme(panel.grid=element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1.3))+ theme_classic()
pdf("3.NVG_Concordance.pdf")
q <- p+xlab("Segment length (Kb)")+ylab("Concordance (%)")+theme(axis.title.x = element_text(face='bold',size=15,hjust=0.5),
                                                                 axis.title.y = element_text(face='bold',size=15,vjust=1),
                                                                 axis.text.x = element_text(face='bold',size=15,color='black'),
                                                                 axis.text.y = element_text(face='bold',size=15,color='black'),
                                                                 legend.title = element_text(face="bold", size=15),
                                                                 legend.text = element_text(face="bold", size=15),
                                                                 plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5))+
  scale_y_continuous(limits=c(0.4,1),breaks=c(0.4,0.6,0.8,1.0))

q
dev.off()										  
