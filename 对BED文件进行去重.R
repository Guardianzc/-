setwd("E:\\Rdata")
dat <- read.table("test.bed")
dat1 <- data.frame()
dat1 <- dat[1,]
n <- nrow(dat)
js <- 0
for (i in 2:n) {
  if (dat[i,1]!=dat[i-1,1]||dat[i,2]!=dat[i-1,2]||dat[i,3]!=dat[i-1,3]) 
  {dat1 <- rbind(dat1,dat[i,])}
}

write.table(dat1,'test2.bed',col.names=FALSE,row.names=FALSE,quote = F)

