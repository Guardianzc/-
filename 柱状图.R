library(ggplot2)
setwd('E://Rdata')
fc <- read.csv('forecast.csv')
data1=na.omit(fc[c("years","Earning")])
label <- fc$label

library(showtext)
showtext.auto(enable = TRUE)
font.add('SimSun', 'simsun.ttc')
pdf(file="forecast.pdf")
  ggplot(data=data1,mapping=aes(x=years,y=Earning,fill=label,group=factor(1)))+
  geom_bar(stat="identity",width = 0.5)+
  geom_text(aes(label = Earning, vjust = -0.8, hjust = 0.5, color = label), show.legend = TRUE)+
  labs(x="年份",y="全国餐饮收入（亿元）",title="全国餐饮收入预测")
dev.off()

