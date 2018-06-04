install.packages('forecast')
library(forecast)
setwd('E://Rdata')
data <- read.csv("forecast.csv")

plot(data$Earning,type='l')

freq <- spec.pgram(data$Earning, taper=0, log='no', plot=FALSE);

start <- which(freq$spec==max(freq$spec))
frequency <- 1/freq$freq[which(freq$spec==max(freq$spec))]

meanTS <- ts(
  data$Earning[start:length(data$Earning)], 
  frequency=frequency
)

meanARIMA = auto.arima(meanTS)
meanARIMAForecast = forecast(meanARIMA, h=5);
meanARIMAForecast$mean
plot(meanARIMAForecast)
