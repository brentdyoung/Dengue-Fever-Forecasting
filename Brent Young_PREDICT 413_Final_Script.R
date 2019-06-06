#Final Project Code
#R Code
#Load the data & create time series
library(fpp)
library(tseries)
library(forecast)
library(ggplot2)
library(forecastHybrid)
library(opera)
library(dplyr)
library(corrplot)
library(zoo)
library(MASS)
library(reshape2)
library(mice)
library(VIM)
library(Hmisc)
library(forecastxgb)
library(fGarch)
library(rugarch)

memory.size(max=TRUE)

#Clear Global Environment

#Load Data 
setwd("~/R/PREDICT 413/Final")
submissions <- read.csv("submission_format.csv")

#Full Training Data (Unclean)
de <- read.csv("dengue_features_train_combined.csv")
View(de)
summary(de)
str(de)
describe(de)

#Data Preparation
library(mice)
sju= de[1:936,]
summary(sju)
str(sju)
describe(sju)

iqu= de[937:1456,]
summary(iqu)
str(iqu)
describe(iqu)

#Missing Data EDA
#Check for missing values
sapply(de, function(x) sum(is.na(x)))
sum(is.na(de))

sapply(sju, function(x) sum(is.na(x)))
sum(is.na(sju))

sapply(iqu, function(x) sum(is.na(x)))
sum(is.na(iqu))

#Check missing data percentage
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(de,2,pMiss)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(sju,2,pMiss)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(iqu,2,pMiss)

library(VIM)
aggr_plot <- aggr(de, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(de), cex.axis=.5, gap=2, ylab=c("Histogram of missing data","Pattern"))

aggr_plot <- aggr(sju, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(de), cex.axis=.5, gap=2, ylab=c("Histogram of missing data","Pattern"))

aggr_plot <- aggr(iqu, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(de), cex.axis=.5, gap=2, ylab=c("Histogram of missing data","Pattern"))

#Split datasets into numerical and categorical

#Numeric
subdatnum <- subset(de, select=c(
  "INDEX",
  "total_cases",
  "ndvi_ne",
  "ndvi_nw",
  "ndvi_se",
  "ndvi_sw",
  "precipitation_amt_mm",
  "reanalysis_air_temp_k",
  "reanalysis_avg_temp_k",
  "reanalysis_dew_point_temp_k",
  "reanalysis_max_air_temp_k",
  "reanalysis_min_air_temp_k",
  "reanalysis_precip_amt_kg_per_m2",
  "reanalysis_relative_humidity_percent",
  "reanalysis_sat_precip_amt_mm",
  "reanalysis_specific_humidity_g_per_kg",
  "reanalysis_tdtr_k",
  "station_avg_temp_c",
  "station_diur_temp_rng_c",
  "station_max_temp_c",
  "station_min_temp_c",
  "station_precip_mm"))

subdatnum.df <- data.frame(subdatnum)

#Categorical
subdatcat <- subset(de, select=c(
  "INDEX",
  "city",
  "year",
  "weekofyear",
  "week_start_date"))

subdatcat.df <- data.frame(subdatcat)

#Run imputation
tempData <- mice(subdatnum.df,m=5,maxit=50,meth='pmm',seed=500)
summary(tempData)

#Check N/A values have been removed
subdatnumimp <- complete(tempData,1)
apply(subdatnumimp,2,pMiss)
summary(subdatnumimp)
sapply(subdatnumimp, function(x) sum(is.na(x)))

#Manually filled reanalysis_sat_precip_amt_mm and precipitation_amt_mm in Excel to 0. Imputation was not picking it up. 

#Merge Numeric and Categorical datasets back
DE <- merge(subdatcat.df,subdatnumimp, by=c("INDEX"))

#Check data
str (DE)
summary (DE)

#Full Training Data for SJ 
sj= DE[1:936,]
summary(sj)

#Full Training Data for IQ
iq= DE[937:1456,]
summary(iq)

sjtimeseries <- ts(sj$total_cases)
iqtimeseries <- ts(iq$total_cases)
iqtimeseries

#Test Data for SJ 
sj_test= DE[1457:1716,]
summary(sj_test)

#Test Data for IQ
iq_test= DE[1717:1872,]
summary(iq_test)

#EDA
#Timeplot on for SJ and IQ
par(mfrow=c(2,1)) 
plot(sjtimeseries, main='SJ Total Cases', xlab = "Week", ylab= "Cases")
plot(iqtimeseries, main='IQ Total Cases', xlab = "Week", ylab= "Cases")
par(mfrow=c(1,1)) 

#Barplot of SJ and IQ of Total Cases
p1<- ggplot(sj, aes(x = total_cases)) +
  geom_vline(xintercept = 29, color = "red") +
  theme_bw() +
  ggtitle('Cases of Dengue in San Juan') +
  geom_histogram(fill = "dark blue", bins=30) +
  theme(plot.title = element_text(hjust = 0.5))
p1 + annotate("text", x = 77, y = 250, label = "Median Cases is 19")

p2<- ggplot(iq, aes(x = total_cases)) +
  geom_vline(xintercept = 5, color = "red") +
  theme_bw() +
  ggtitle('Cases of Dengue in Iquitos') +
  geom_histogram(fill = "light blue", bins = 30) +
  theme(plot.title = element_text(hjust = 0.5))
p2 + annotate("text", x = 17, y = 190, label = "Median Cases is 5")

#Correlation Plots 
m_sj_clean <- data.matrix(sj)
m_sj_clean <- cor(x = m_sj_clean [,c(6:26)], use = 'complete.obs', method = 'pearson')

m_iq_clean <- data.matrix(iq)
m_iq_clean <- cor(x = m_iq_clean [,c(6:26)], use = 'everything', method = 'pearson')

# Correlation Heatmap
corrplot(m_sj_clean, type = 'full', tl.col = 'black', method="shade", shade.col=NA, tl.cex=0.5)
corrplot(m_iq_clean, type = 'full', tl.col = 'black', method="shade", shade.col=NA,tl.cex=0.5)

# Correlation Bar plot
df_m_sj_clean <- data.frame(m_sj_clean)[2:21,] 
df_m_sj_clean <- dplyr::select(df_m_sj_clean, total_cases) 

df_m_iq_clean <- data.frame(m_iq_clean)[2:21,]
df_m_iq_clean <- dplyr::select(df_m_iq_clean, total_cases) 

ggplot(df_m_sj_clean, aes(x= reorder(rownames(df_m_sj_clean), -total_cases), y = total_cases)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  ggtitle('Correlation of variables in San Juan') +
  ylab('Correlation') +
  xlab('Variables') +
  coord_flip()

ggplot(df_m_iq_clean, aes(x= reorder(rownames(df_m_iq_clean), -total_cases), y = total_cases)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  ggtitle('Correlation of variables in Iquitos') +
  ylab('Correlation') +
  xlab('Variables') +
  coord_flip()

#Split Data into Train/Test for both SJ and IQ using train (clean) dataset

#SJ - 260 for test to match test dataset
SJtimeseries_train = DE[1:675,]
SJtimeseries_test = DE[676:936,]

sjtimeseries_train <- ts(SJtimeseries_train$total_cases)
sjtimeseries_test <-ts(SJtimeseries_test$total_cases)

#IQ - 156 for test to match test dataset
IQtimeseries_train = DE[937:1299,]
IQtimeseries_test = DE[1300:1456,]

iqtimeseries_train <- ts(IQtimeseries_train$total_cases)
iqtimeseries_test <-ts(IQtimeseries_test$total_cases)

#Model Selection for SJ
#ETS Model
fit1_ets <- ets(sjtimeseries_train)
summary(fit1_ets)
#Auto.Arima Pre-work
tsdisplay(sjtimeseries_train)

#Unit Root Tests
adf.test(sjtimeseries_train, alternative = "stationary")
kpss.test(sjtimeseries_train)

#Differencing 
kpss.test(diff(sjtimeseries_train))
tsdisplay(diff(sjtimeseries_train))

#Auto.Arima Model 
fit2_arima <- auto.arima(sjtimeseries_train, stepwise=FALSE, approximation=FALSE)
summary(fit2_arima)

tsdisplay(residuals(fit2_arima))
Box.test(residuals(fit2_arima), fitdf=3, lag=10, type="Ljung")
tsdiag(fit2_arima)

#Auto.Arima Model w/Regressors
fit2_arimaR <- auto.arima(SJtimeseries_train [,6], xreg= SJtimeseries_train [,c(12:16,20,22,24:25)], stepwise=FALSE, approximation=FALSE)
summary(fit2_arimaR)

tsdisplay(residuals(fit2_arimaR))
Box.test(residuals(fit2_arimaR), fitdf=3, lag=10, type="Ljung")
tsdiag(fit2_arimaR)

#Neutral Net 
fit3_nn <- nnetar(sjtimeseries_train)
summary(fit3_nn)

#Neutral Net w/Regressors
fit3_nnR <- nnetar(SJtimeseries_train [,6], xreg= SJtimeseries_train [,c(4,16,20)])
summary(fit3_nnR)

#Benchmark
fit6_naive<-naive(sjtimeseries_train)

#Training Set Accuracy - Summary
summary(fit1_ets) # training set
summary(fit2_arima) # training set
summary(fit2_arimaR) # training set
summary(fit3_nn) # training set
summary(fit3_nnR) # training set
summary(fit6_naive) # training set

#Training Set Accuracy - Goodness-of-fit
accuracy (fit1_ets) # training set
accuracy (fit2_arima) # training set
accuracy (fit2_arimaR) # training set
accuracy (fit3_nn) # training set
accuracy (fit3_nnR) # training set
accuracy (fit6_naive) # training set

#Forecast on Test Set
par(mfrow=c(3,2)) 
ETS <-forecast(fit1_ets, h=length(sjtimeseries_test))
plot(ETS, ylab="Total Cases")
lines(sjtimeseries, col="red",ylab="Actual")
ETS

Auto.ARIMA <-forecast(fit2_arima, h=length(sjtimeseries_test))
plot(Auto.ARIMA, ylab="Total Cases")
lines(sjtimeseries, col="red",ylab="Actual")
Auto.ARIMA

Auto.ARIMAR <- forecast(object=fit2_arimaR, xreg =SJtimeseries_train[,c(12:16,20,22,24:25)],h=260)
plot(Auto.ARIMAR,
     main="Forecasts from regression ", ylab="Total Cases")
lines(sjtimeseries, col="red",ylab="Actual")
Auto.ARIMAR

NN <-forecast(fit3_nn, h=length(sjtimeseries_test))
plot(NN, ylab="Total Cases")
lines(sjtimeseries, col="red",ylab="Actual")
NN

NNR <-forecast(object=fit3_nnR, xreg =SJtimeseries_train[,c(4,16,20)],h=260)
plot(NNR, ylab="Total Cases")
lines(sjtimeseries, col="red",ylab="Actual")
NNR

NAIVE <-naive(sjtimeseries_train, h=length(sjtimeseries_test))
plot(NAIVE)
lines(sjtimeseries, col="red",ylab="Actual")

par(mfrow=c(1,1)) 

print(accuracy(ETS, sjtimeseries))
print(accuracy(Auto.ARIMA, sjtimeseries))
print(accuracy(Auto.ARIMAR, sjtimeseries))
print(accuracy(NN, sjtimeseries))
print(accuracy(NNR, sjtimeseries))
print(accuracy(NAIVE, sjtimeseries))

#Model Selection for IQ
#ETS Model
fit1_ets <- ets(iqtimeseries_train)
summary(fit1_ets)

#Auto.Arima Pre-work
tsdisplay(iqtimeseries_train)

#Unit Root Tests
adf.test(iqtimeseries_train, alternative = "stationary")
kpss.test(iqtimeseries_train)

#Differencing 
kpss.test(diff(iqtimeseries_train))
tsdisplay(diff(iqtimeseries_train))

#Auto.Arima Model 
fit2_arima <- auto.arima(iqtimeseries_train, stepwise=FALSE, approximation=FALSE)
summary(fit2_arima)

tsdisplay(residuals(fit2_arima))
Box.test(residuals(fit2_arima), fitdf=3, lag=10, type="Ljung")
tsdiag(fit2_arima)

#Auto.Arima Model w/Regressors
fit2_arimaR <- auto.arima(IQtimeseries_train [,6], xreg= IQtimeseries_train [,c(14,16,20,25)], stepwise=FALSE, approximation=FALSE)
summary(fit2_arimaR)

tsdisplay(residuals(fit2_arimaR))
Box.test(residuals(fit2_arimaR), fitdf=3, lag=10, type="Ljung")
tsdiag(fit2_arimaR)

#Neutral Net 
fit3_nn <- nnetar(iqtimeseries_train)
summary(fit3_nn)

#Neutral Net w/Regressors
fit3_nnR <- nnetar(IQtimeseries_train [,6], xreg= IQtimeseries_train [,c(4,16,20)])
summary(fit3_nnR)

#Benchmark
fit6_naive<-naive(iqtimeseries_train)

#Training Set Accuracy - Summary
summary(fit1_ets) # training set
summary(fit2_arima) # training set
summary(fit2_arimaR) # training set
summary(fit3_nn) # training set
summary(fit3_nn) # training set
summary(fit6_naive) # training set

#Training Set Accuracy - Goodness-of-fit
accuracy(fit1_ets) # training set
accuracy (fit2_arima) # training set
accuracy (fit2_arimaR) # training set
accuracy (fit3_nn) # training set
accuracy (fit3_nnR) # training set
accuracy (fit6_naive) # training set

#Forecast on Test Set
par(mfrow=c(3,2)) 
ETS_ANN <-forecast(fit1_ets, h=length(iqtimeseries_test))
plot(ETS_ANN, ylab="Total Cases")
lines(iqtimeseries, col="red",ylab="Actual")
ETS_ANN

Auto.ARIMA <-forecast(fit2_arima, h=length(iqtimeseries_test))
plot(Auto.ARIMA, ylab="Total Cases")
lines(iqtimeseries, col="red",ylab="Actual")
Auto.ARIMA

Auto.ARIMAR <- forecast(object=fit2_arimaR, xreg =IQtimeseries_train[,c(14,16,20,25)],h=156)
plot(Auto.ARIMAR,
     main="Forecasts from regression ", ylab="Total Cases")
lines(iqtimeseries, col="red",ylab="Actual")
Auto.ARIMAR

NN <-forecast(fit3_nn, h=length(iqtimeseries_test))
plot(NN, ylab="Total Cases")
lines(iqtimeseries, col="red",ylab="Actual")
NN

NNR <-forecast(object=fit3_nnR, xreg = IQtimeseries_train [,c(4,16,20)],h=156)
plot(NNR, ylab="Total Cases")
lines(iqtimeseries, col="red",ylab="Actual")
NNR

NAIVE <-naive(iqtimeseries_train, h=length(iqtimeseries_test))
plot(NAIVE)
lines(iqtimeseries, col="red",ylab="Actual")

par(mfrow=c(1,1))

print(accuracy(ETS_ANN, iqtimeseries))
print(accuracy(Auto.ARIMA, iqtimeseries))
print(accuracy(Auto.ARIMAR, iqtimeseries))
print(accuracy(NN, iqtimeseries))
print(accuracy(NNR, iqtimeseries))
print(accuracy(NAIVE, iqtimeseries))

#Final Submission: Auto.ARIMA
fit_sj <- auto.arima(sjtimeseries)
fit_Arima_sj=forecast(fit_sj, h = 260)
plot(forecast(fit_sj))
tsdisplay(residuals(fit_sj))
Box.test(residuals(fit_Arima_sj),fitdf=2, lag=151, type="Ljung")  

fit_iq <- auto.arima(iqtimeseries)
fit_Arima_iq=forecast(fit_iq, h = 156)
plot(forecast(fit_Arima_iq))
tsdisplay(residuals(fit_Arima_iq))
Box.test(residuals(fit_Arima_iq), fitdf=3, lag=151, type="Ljung")

arima_sj_sol <- data.frame(submissions[1:260,-4], total_cases = round(fit_Arima_sj$mean))
arima_iq_sol <- data.frame(submissions[261:416,-4], total_cases =round(fit_Arima_iq$mean))
arima_solution <- bind_rows(arima_sj_sol,arima_iq_sol)

write.csv(arima_solution, file = 'submission.csv', row.names = F)

#Final Submission: ETS
fit_sj <- ets(sjtimeseries)
fit_ETS_sj=forecast(fit_sj, h = 260)
plot(forecast(fit_sj))
tsdisplay(residuals(fit_sj))
Box.test(residuals(fit_ETS_sj),fitdf=2, lag=151, type="Ljung")  

fit_iq <- ets(iqtimeseries)
fit_ETS_iq=forecast(fit_iq, h = 156)
plot(forecast(fit_ETS_iq))
tsdisplay(residuals(fit_ETS_iq))
Box.test(residuals(fit_ETS_iq), fitdf=3, lag=151, type="Ljung")

ETS_sj_sol <- data.frame(submissions[1:260,-4], total_cases = round(fit_ETS_sj$mean))
ETS_iq_sol <- data.frame(submissions[261:416,-4], total_cases =round(fit_ETS_iq$mean))
ETS_solution <- bind_rows(ETS_sj_sol,ETS_iq_sol)

write.csv(ETS_solution, file = 'submission.csv', row.names = F)

#Final Submission: NNETAR 10,6 and 5,3
set.seed(101)
fit_sj <- nnetar(sjtimeseries)
fit_Nnetar_sj=forecast(fit_sj, h = 260)
plot(forecast(fit_sj))
tsdisplay(residuals(fit_Nnetar_sj))
Box.test(residuals(fit_Nnetar_sj),fitdf=2, lag=151, type="Ljung")  

fit_iq <- nnetar(iqtimeseries)
fit_Nnetar_iq=forecast(fit_iq, h = 156)
plot(forecast(fit_Nnetar_iq))
tsdisplay(residuals(fit_Nnetar_iq))
Box.test(residuals(fit_Nnetar_iq), fitdf=3, lag=151, type="Ljung")

nnetar_sj_sol <- data.frame(submissions[1:260,-4], total_cases = round(fit_Nnetar_sj$mean))
nnetar_iq_sol <- data.frame(submissions[261:416,-4], total_cases =round(fit_Nnetar_iq$mean))
nnetar_solution <- bind_rows(nnetar_sj_sol,nnetar_iq_sol)

write.csv(nnetar_solution, file = 'submission.csv', row.names = F)

#Final Submission: ARIMAR 
fit_sj <- auto.arima(sj[,6], xreg= sj[,c(12:16,20,22,24:25)], stepwise=FALSE, approximation=FALSE)
fit_Arima_sj=forecast(object=fit_sj, xreg =sj_test[,c(12:16,20,22,24:25)],h=260)
plot(forecast(fit_Arima_sj))
tsdisplay(residuals(fit_sj))
Box.test(residuals(fit_Arima_sj),fitdf=2, lag=151, type="Ljung")  

fit_iq <- auto.arima(iq[,6], xreg= iq[,c(14,16,20,25)], stepwise=FALSE, approximation=FALSE)
fit_Arima_iq=forecast(object=fit_iq, xreg =iq_test[,c(14,16,20,25)],h=156)
plot(forecast(fit_Arima_iq))
tsdisplay(residuals(fit_Arima_iq))
Box.test(residuals(fit_Arima_iq), fitdf=3, lag=151, type="Ljung")

arima_sj_sol <- data.frame(submissions[1:260,-4], total_cases = round(fit_Arima_sj$mean))
arima_iq_sol <- data.frame(submissions[261:416,-4], total_cases =round(fit_Arima_iq$mean))
arima_solution <- bind_rows(arima_sj_sol,arima_iq_sol)

write.csv(arima_solution, file = 'submission.csv', row.names = F)

#Final Submission: NNETAR with Regressors (BEST)
set.seed(101)
fit_sj <- nnetar(sj[,6], xreg= sj[,c(4,16,20)])
fit_Nnetar_sj=forecast(object=fit_sj, xreg =sj_test[,c(4,16,20)],h=260)
plot(forecast(fit_Nnetar_sj))
tsdisplay(residuals(fit_sj))
Box.test(residuals(fit_Nnetar_sj),fitdf=2, lag=151, type="Ljung")  

fit_iq <- nnetar(iq[,6], xreg= iq[,c(4,16,20)])
fit_Nnetar_iq= forecast(object=fit_iq, xreg =iq_test[,c(4,16,20)],h=156)
plot(forecast(fit_Nnetar_iq))
tsdisplay(residuals(fit_Nnetar_iq))
Box.test(residuals(fit_Nnetar_iq), fitdf=3, lag=151, type="Ljung")

nnetar_sj_sol <- data.frame(submissions[1:260,-4], total_cases = round(fit_Nnetar_sj$mean))
nnetar_iq_sol <- data.frame(submissions[261:416,-4], total_cases =round(fit_Nnetar_iq$mean))
nnetar_solution <- bind_rows(nnetar_sj_sol,nnetar_iq_sol)

write.csv(nnetar_solution, file ='submission.csv', row.names = F)
