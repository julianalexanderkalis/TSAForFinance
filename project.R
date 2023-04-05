# 1 ----- loading required libraries

library(tseries)
library(forecast)
library(lmtest)
library(quantmod)
library(tidyverse)
library(vars)
library(fGarch)


# 2 ----- Setup


# set working directory

setwd("...") # use own wd


# set seed for repeatable results
set.seed(420)


# read in csv data and remove columns with NA values
ecto <- read_csv("gw2_glob_of_ecto.csv",col_types = cols())
ecto <- ecto[ , colSums(is.na(ecto)) < nrow(ecto)]

coin <- read_csv("gw2_mystic_coin.csv",col_types = cols())
coin <- coin[ , colSums(is.na(coin)) < nrow(coin)]


# create ts object from buy_price_max
ts_coin <- ts(coin$buy_price_max, start = 2013, frequency = 365, end = 2021)
ts_ecto <- ts(ecto$buy_price_max, start = 2013, frequency = 365, end = 2021)


# calculating returns
coin_return <- na.omit(diff(log(ts_coin)))
ecto_return <- na.omit(diff(log(ts_ecto)))



# Plotting the price development over the years
par(mfrow=c(2,1))
plot(ts_coin, type = 'l', ylab = 'Price of Mystic Coin', xlab = 'Years')
plot(ts_ecto, type='l', ylab = 'Price of Ectoplasm', xlab = 'Years')


# testing for stationarity
adf.test(ts_coin)
adf.test(ts_ecto)

adf.test(coin_return)
adf.test(ecto_return)


# 3 ----- Analysis of Seasonal Trends


# coin
coin_decomposed <- decompose(ts_coin)
plot(coin_decomposed, xlab="Year")
seasonality_coin <- coin_decomposed$seasonal
trend_coin <- coin_decomposed$trend


# ecto
ecto_decomposed <- decompose(ts_ecto)
plot(ecto_decomposed, xlab="Year")
seasonality_ecto <- ecto_decomposed$seasonal
trend_ecto <- ecto_decomposed$trend


# --- plots


# seasonality
par(mfrow=c(2,1))
plot(seasonality_coin, main="Seasonality of Mystic Coin", ylab="In copper coins", xlab="Year") +
  grid()
plot(seasonality_ecto, main="Seasonality of Ectoplasm", ylab="In copper coins", xlab="Year") +
  grid()

# trend
par(mfrow=c(2,1))
plot(trend_coin, main="Trend of Mystic Coin", ylab="In copper coins", xlab="Year") +
  grid()
plot(trend_ecto, main="Trend of Ectoplasm", ylab="In copper coins", xlab="Year") +
  grid()


# 4 ----- correlation of the two resources


# in levels
par(mfrow=c(1,1))
plot(x=ts_ecto, y=ts_coin, xlab = "Ectoplasm", ylab = "Mystic Coin", pch = 20) + 
  grid()
cor(ts_ecto,ts_coin)
abline(lm(ts_coin~ts_ecto), col='blue')

# in returns
par(mfrow=c(1,1))
plot(x=ecto_return, y=coin_return, xlab = "Ectoplasm", ylab = "Mystic Coin") + 
  grid()
cor(ecto_return, coin_return)


# 5 ----- VAR analysis


VAR_ecto_coin <- VAR(cbind(coin_return, ecto_return), ic="AIC", lag.max = 10)

# testing for causality using granger criteria
causality(VAR_ecto_coin, cause="coin_return")["Granger"]
causality(VAR_ecto_coin, cause="ecto_return")["Granger"]

# testing using the Impulse Response Function
ecto_impulse <- irf(VAR_ecto_coin, impulse="ecto_return", response="coin_return")
coin_impulse <- irf(VAR_ecto_coin, impulse="coin_return", response="ecto_return")

# plotting the irf
plot(ecto_impulse)
plot(coin_impulse)


# 6 ----- forecasting


# arima coin
arima_coin <- auto.arima(log(ts_coin), ic="aic") 
pred_arima_coin <- forecast(arima_coin,level=0.95,h=200)

# arima ecto
arima_ecto <- auto.arima(log(ts_ecto), ic="aic") 
pred_arima_ecto <- forecast(arima_ecto,level=0.95,h=250)

# plotting arima
par(mfrow=c(2,1))
plot(pred_arima_coin, ylab="Mystic Coin (log)") + grid()
plot(pred_arima_ecto, ylab="Ectoplasm (log)") + grid()

# garch model
par(mfrow=c(1,1))

gar_model <- garchFit(~arma(1,1) + garch(1,1),data=ecto_return, trace = F)
gar_pre_ecto <- predict(gar_model, n.ahead = 5, plot=TRUE, nx=50)

gar_model <- garchFit(~arma(1,1) + garch(1,1),data=coin_return, trace = F)
gar_pre_coin <- predict(gar_model, n.ahead = 5, plot=TRUE, nx=50)


# ---- best parameters for Ectoplasm, due to time for calculations I focused on one and not both items in the ...
# ... paper, the other one can be found as a picture in the appendix

# subsetting timeseries
best_param_ecto <- ecto_return
best_param_test <- ts(best_param_ecto[2914:2920], start = 2021, frequency = 365)
best_param_train <- ts(best_param_ecto[0:2913], start = 2013, frequency = 365, end = 2021)

# initializing values
best_value <- 1
best_values <- c(0,0)

# running for loops over model, second parameter returned an 'infinite' error for certain combinations when ...
# ... increased over 8, thus to achieve consistency I limited it to this value, but since the optimal ...
# ... solutions never reached the limit I deemed it acceptable
for(i in 1:30){
  for(j in 1:8){

    # passing parameters into fit function using substitution, otherwise it doesn't work for some reason
    gar_model <- garchFit(substitute(~ arma(1,u)+garch(v,1),list(u=i, v=j)),data = best_param_train,trace = F)
    condPre <- predict(gar_model, n.ahead = 7, plot=F, nx=50)
    
    # calculating metric
    diff_pred_test <- mean(abs(condPre$meanForecast - best_param_test))
    
    # print combination, can be ignored, was just for seeing if all combinations are really tested
    print(c(i,j))
    
    # if new calculation is better then the one previously saved, replace combination
    if(diff_pred_test < best_value) {
      best_value <- diff_pred_test
      best_values <- c(i, j)
      
      # print new best value, can be ignored, was just for seeing if it works
      print('----- New best values: -----')
      print(best_values)
      print('----- continuing -----')
    }
  }
  
}

# plotting the optimal model
gar_model <- garchFit(formula = ~arma(1,27)+garch(1,1), data = best_param_train, trace = F)
condPre <- predict(gar_model, n.ahead = 7, plot=TRUE, nx=50)





# ----- Appendix 





# subsetting timeseries
best_param_coin <- coin_return
best_param_test <- ts(best_param_coin[2914:2920], start = 2021, frequency = 365)
best_param_train <- ts(best_param_coin[0:2913], start = 2013, frequency = 365, end = 2021)

# initializing values
best_value <- 1
best_values <- c(0,0)

# running for loops over model, second parameter returned an 'infinite' error for certain combinations when ...
# ... increased over 8, thus to achieve consistency I limited it to this value, but since the optimal ...
# ... solutions never reached the limit I deemed it acceptable
for(i in 1:30){
  for(j in 1:8){
    
    # passing parameters into fit function using substitution, otherwise it doesn't work for some reason
    gar_model <- garchFit(substitute(~ arma(1,u)+garch(v,1),list(u=i, v=j)),data = best_param_train,trace = F)
    condPre <- predict(gar_model, n.ahead = 7, plot=F, nx=50)
    
    # calculating metric
    diff_pred_test <- mean(abs(condPre$meanForecast - best_param_test))
    
    # print combination, can be ignored, was just for seeing if all combinations are really tested
    print(c(i,j))
    
    # if new calculation is better then the one previously saved, replace combination
    if(diff_pred_test < best_value) {
      best_value <- diff_pred_test
      best_values <- c(i, j)
      
      # print new best value, can be ignored, was just for seeing if it works
      print('----- New best values: -----')
      print(best_values)
      print('----- continuing -----')
    }
  }
  
}

# -> opitimal combination: c(27, 3)

# plotting the optimal model
gar_model <- garchFit(formula = ~arma(1,27)+garch(3,1), data = best_param_train, trace = F)
condPre <- predict(gar_model, n.ahead = 7, plot=TRUE, nx=50)
