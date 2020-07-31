library(tidyverse)
library(forecast)
library(expsmooth)
library(TTR)
library(fpp)
library(fpp2)


# Time series demand example
len <- 100

y.err <- rnorm(len, 0 , 5)
y.const <- rep(50, len)
y.trend <- 1:len*0.1
x.seas <- 0:(len-1)%%4 + 1 
y.seas <- (x.seas-1)*2 # quarterly fluctuations

x <- 1:len

#png(file = "forecast_ts_const.png", bg = "transparent", width=600, height = 400)
win.metafile("forecast_ts_const.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x , y.const , type = "l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
dev.off()


#png(file = "forecast_ts_trend.png", bg = "transparent", width=600, height = 400)
win.metafile("forecast_ts_trend.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x , y.trend , type = "l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
dev.off()

#png(file = "forecast_ts_seas.png", bg = "transparent", width=600, height = 400)
win.metafile("forecast_ts_seas.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x , y.seas , type = "l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
dev.off()

win.metafile("forecast_ts_err.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x , y.err , type = "l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
dev.off()


#png(file = "forecast_ts_all.png", bg = "transparent", width=600, height = 400)
win.metafile("forecast_ts_all.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x , y.const + y.trend + y.seas , type = "l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
dev.off()

win.metafile("forecast_ts_all_err.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x ,y.err +  y.const + y.trend + y.seas , type = "l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
dev.off()


# ts values
y.com <- y.err +  y.const + y.trend + y.seas
y.com.ts <- ts(y.com, deltat = 1/4)

# TS decomosition by stl
stl.y <- stl(y.com.ts, s.window = 4)
plot(stl.y)

# TS decomosition by decompose
decompose.y <- decompose(y.com.ts)
plot(decompose.y)

# TS decomosition by lm
y.lm.decom <- lm(y.com ~ x + as.factor(x.seas) )
summary(y.lm.decom)
anova(y.lm.decom)

y.lm.pred <- predict(y.lm.decom)
y.lm.res <- residuals(y.lm.decom)

# MSE
sum(y.lm.res^2)/100
# MAD
median(abs(y.lm.res))
# MAPE
mean(abs(y.lm.res)/y.lm.pred)
# log q
mean(log(y.lm.pred/y.com)^2)

# function for calculating accuracy measures
fc.stats <- function(yhat, y){
  y.res <- y -  yhat
  # MSE
  mse <- mean(y.res^2, na.rm =T)
  # MAD
  mad <- median(abs(y.res), na.rm =T)
  # MAPE
  mape <- mean(abs(y.res)/yhat, na.rm =T)
  # log q
  msla <- mean(log(yhat/y)^2, na.rm =T)
  #   return
  res <- c(mse,mad,mape,msla)
  names(res) <- c("MSE","MAD","MAPE","MSLA")
  return(res)
  }

# test functions
fc.stats(yhat = y.lm.pred, y = y.com)


# Regression analysis of examplary time series  ##############################

y.lm.trend <- coefficients(y.lm.decom)[2] * x
lm.seas.coeff <- c(0,tail(coefficients(y.lm.decom),3))
y.lm.seas <- rep(lm.seas.coeff, ceiling(len/4))[1:len]

win.metafile("forecast_ts_lm_trend.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.lm.trend, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
lines(x, y.trend, type="l", col="red", lwd = 1.5)
legend("topleft", col=c("red","black"), legend = c("true model","estimate"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()

win.metafile("forecast_ts_lm_seas.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.lm.seas, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
lines(x, y.seas, type="l", col="red")
legend("topleft", col=c("red","black"), legend = c("true model","estimate"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()


win.metafile("forecast_ts_lm_pred.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.com, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.5 , cex.lab = 2, cex.axis = 1.5)
lines(x, y.lm.pred, type="l", col="red")
legend("topleft", col=c("red","black"), legend = c("estimate","observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()

win.metafile("forecast_ts_lm_err.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, residuals(y.lm.decom), type="h", xlab= "time/periods", ylab ="residual", lwd = 1.25 , cex.lab = 2, cex.axis = 1.5, pch =15)
points(x-.25, y.err, type="h", col="red", lwd=1.25)
legend("topleft", col=c("red","black"), legend = c("estimate","observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()

# moving average ##############################

y.sma.4 <- SMA(y.com.ts, n =4)

win.metafile("forecast_ts_sma_4.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.com.ts, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.75 , cex.lab = 2, cex.axis = 1.5, pch =15)
lines(x, c(NA, head(y.sma.4,99)), col="red", lwd=1.75)
legend("topleft", col=c("red","black"), legend = c("esti. SMA (n=4)","observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()

fc.stats(yhat = y.sma.4, y = y.com)

# find best moving average parameter
sim.sma <- function(x.n , y=y.com.ts){
  # calculates SMA statistics for a vector of time windows x.n
  sapply(x.n, function(x){
    tmp <- SMA(y, n =x)
    fc.stats(yhat = c(NA, head(tmp,99)), y = as.numeric(y))
  })
}
# Results
x.n <- 2:20
res.sma <- sim.sma(x.n = x.n)
# best n w.r.t. MSE
x.n[which.min(res.sma["MSE",])]
# optimal statistics
res.sma[,which.min(res.sma["MSE",])]
res.sma[,which.min(res.sma["MAPE",])]


win.metafile("forecast_ts_sma_best.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,5))
plot(x.n, res.sma["MSE",], type="h", xlab= expression(n), ylab ="MSE", lwd = 2 , cex.lab = 2, cex.axis = 1.5)
par(new = TRUE)
plot(x.n+.25, res.sma["MAPE",], type="h", col="red", lwd=2, xaxt = "n", yaxt = "n", xlab="", ylab="", xlim=c(2,20))
axis(side = 4, cex.axis = 1.5)
mtext("MAPE", side = 4, line = 3, cex = 2)
legend("topright", col=c("red","black"), legend = c("MAPE","MSE"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()


# 1st order exp. smoothing #####################################################################

#y.com <- y.com.ts <- c(50.46,49.78,58.39,61.97,54.33,55.67,60.07,50.33,60.45,52.21)

y.1st.exp.smoo <- HoltWinters(y.com.ts, alpha = .4, beta = F, gamma = F)

round(head(cbind(y.com, c(NA, y.1st.exp.smoo$fitted[,1])),10),2)


win.metafile("forecast_ts_1stem_04.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.com.ts, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.75 , cex.lab = 2, cex.axis = 1.5, pch =15)
lines(x, c(NA, y.1st.exp.smoo$fitted[,1]), col="red", lwd=1.75)
legend("topleft", col=c("red","black"), legend = c(expression(paste("1st ES with (", alpha==0.4,")")),"observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()


fc.stats(yhat = c(NA, y.1st.exp.smoo$fitted[,1]), y=y.com)

# find best alpha 
sim.1st.es <- function(x.alpha , y=y.com.ts){
  # calculates 1stES statistics for a vector of alphas x.alpha
  sapply(x.alpha, function(x){
    tmp <- HoltWinters(y, alpha = x, beta = F, gamma = F)
    fc.stats(yhat = c(NA, tmp$fitted[,1]), y = as.numeric(y))
  })
}
# results
x.alp <- seq(0.05,.95, length.out = 100)
res.sim.1stes <- sim.1st.es(x.alpha = x.alp)
# optimal alphas
x.alp[which.min(res.sim.1stes["MSE",])]
x.alp[which.min(res.sim.1stes["MAPE",])]
# optimal accuracy measures
res.sim.1stes[,which.min(res.sim.1stes["MSE",])]
res.sim.1stes[,which.min(res.sim.1stes["MAPE",])]


win.metafile("forecast_ts_1stem_best.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,5))
plot(x.alp, res.sim.1stes["MSE",], type="l", xlab= expression(alpha), ylab ="MSE", lwd = 2 , cex.lab = 2, cex.axis = 1.5)
par(new = TRUE)
plot(x.alp, res.sim.1stes["MAPE",], type="l", col="red", lwd=2, xaxt = "n", yaxt = "n", xlab="", ylab="")
axis(side = 4, cex.axis = 1.5)
mtext("MAPE", side = 4, line = 3, cex = 2)
legend("topleft", col=c("red","black"), legend = c("MAPE","MSE"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()


# 2nd order exponential smoothing #######################################

# function for calculating 2ndES forecasts
sec.es <- function(y, alpha = .4, beta = .6, initial = c(mean(y), 0 )){
  
  n <- length(y)
  res <- matrix(NA, ncol=4, nrow=n+2)
  rownames(res) <- 0:(n+1)
  colnames(res) <- c("y","a","b","y.hat")
  
  res["0", c("a","b")] <- initial
  res[2:(n+1),"y"] <- y
  
  for(i in 2:(nrow(res)-1) ){
    res[i, "y.hat"] <- res[i-1, "a"] + res[i-1, "b"]
    res[i, "a"] <- alpha * res[i, "y"] + (1 - alpha) * res[i, "y.hat"]
    res[i, "b"] <- beta * (res[i, "a"]-res[i-1, "a"]) + (1 - beta) * res[i-1, "b"]
  }
  res[n+2, "y.hat"] <- res[n+1, "a"] + res[n+1, "b"]
  return(res)
}

y.2nd.exp.smoo <- sec.es(y.com, alpha = .4, beta = 0.2, initial = c(50,0))
fc.stats(yhat = y.2nd.exp.smoo[2:101,"y.hat"] , y = y.com)

# alternative from stats package (initial values cannot be controlled)
# y.2nd.exp.smoo <- HoltWinters(y.com.ts, alpha = T, beta = T, gamma = F, l.start = 50, b.start = 0, start.periods = 0)

# alternative from forecast package 
y.2nd.exp.smoo <- holt(y.com.ts , h = 4 , initial = "simple")
fc.stats(yhat = y.2nd.exp.smoo$model$fitted , y = y.com)

#png(file = "forecast_ts_2ndem_best.png", bg = "transparent", width=600, height = 400)
win.metafile("forecast_ts_2ndem_best.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.com.ts, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.75 , cex.lab = 2, cex.axis = 1.5, pch =15)
lines(x, y.2nd.exp.smoo$model$fitted, col="red", lwd=1.75)
legend("topleft", col=c("red","black"), legend = c(expression(paste("2nd ES (", alpha==0.08,",",beta==0.18,")")),"observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()


# 3rd order exponential smoothing
#######################################

# function for calculating 3rdES forecasts
thi.es <- function(y, alpha = .4, beta = .2, gamma=.3, periods = 4, initial = list( a = mean(y), b= 0, c= rep(0,periods) )){
  n <- length(y)
  res <- matrix(NA, ncol=5, nrow=n+1+periods)
  rownames(res) <- (-periods+1):(n+1)
  colnames(res) <- c("y","a","b","c","y.hat")
  
  res["0", c("a","b")] <- c(initial$a, initial$b)
  res[1:periods, "c"] <- initial$c
  res[(periods+1):(n+periods),"y"] <- y
  
  for(i in (periods+1):(nrow(res)-1) ){
    res[i, "y.hat"] <- res[i-1, "a"] + res[i-1, "b"]+ res[i-periods, "c"]
    res[i, "a"] <- alpha * (res[i, "y"] - res[i-periods, "c"]) + (1 - alpha) * (res[i-1, "a"] + res[i-1, "b"])
    res[i, "b"] <- beta * (res[i, "a"] - res[i-1, "a"]) + (1 - beta) * res[i-1, "b"]
    res[i, "c"] <- gamma * (res[i, "y"] - res[i-1, "a"] - res[i-1, "b"]) + (1 - gamma) *res[i-periods, "c"]
  }
  res[nrow(res), "y.hat"] <- res[nrow(res)-1, "a"] + res[nrow(res)-1, "b"] + res[nrow(res)-periods, "c"]
  return(res)
}

# alternative from forecast package 
y.3rd.exp.smoo <- hw(y.com.ts , h = 1 , initial = "simple", seasonal ="additive", alpha = 0.4, beta = 0.2, gamma = 0.3)

# optimal 3rd ES
y.3rd.exp.smoo <- hw(y.com.ts , h = 1 , initial = "simple", seasonal ="additive")

fc.stats(yhat = y.3rd.exp.smoo$model$fitted , y = y.com)

#png(file = "forecast_ts_3rdem_best.png", bg = "transparent", width=600, height = 400)
win.metafile("forecast_ts_3rdem_best.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.com.ts, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.75 , cex.lab = 2, cex.axis = 1.5, pch =15)
lines(x, y.3rd.exp.smoo$model$fitted, col="red", lwd=1.75)
legend("topleft", col=c("red","black"), legend = c(expression(paste("3rd ES (", alpha==0.037,",",beta==0.031,",",gamma==0.089,")")),"observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()

##########################
# best smoothing model - AICc minimization

# determine optimal ES model 
best.smooth <- ets(y.com.ts)

fc.stats(yhat = best.smooth$fitted , y = y.com)

#png(file = "forecast_ts_ets_best.png", bg = "transparent", width=600, height = 400)
win.metafile("forecast_ts_ets_best.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.com.ts, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.75 , cex.lab = 2, cex.axis = 1.5, pch =15)
lines(x, best.smooth$fitted, col="red", lwd=1.75)
legend("topleft", col=c("red","black"), legend = c(expression(paste("opt. ES (", alpha==0,",",beta==0,",",gamma==0,")")),"observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()


##################
# Arima
##################

# optimal arima model
best.arima <- auto.arima(y.com.ts) 
fc.stats(yhat = best.arima$fitted , y = y.com)

checkresiduals(best.arima)
autoplot(forecast(best.arima))

win.metafile("forecast_ts_arima_best.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.com.ts, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.75 , cex.lab = 2, cex.axis = 1.5, pch =15)
lines(x, best.arima$fitted, col="red", lwd=1.75)
legend("topleft", col=c("red","black"), legend = c(expression(paste("opt. ARIMA (3,1,0)",group("(",list(0,0,2),")")[4])),"observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()


# optimal arima model with external regressors
best.arimax <- auto.arima(y.com.ts, xreg = 1:100) 

win.metafile("forecast_ts_arimax_best.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(x, y.com.ts, type="l", xlab= "time/periods", ylab ="demand", lwd = 1.75 , cex.lab = 2, cex.axis = 1.5, pch =15)
lines(x, best.arimax$fitted, col="red", lwd=1.75)
legend("topleft", col=c("red","black"), legend = c(expression(paste("opt. ARIMA (2,0,2)",group("(",list(1,0,0),")")[4])),"observation"), lwd=c(2,2), cex=1.75, bg ="white")
dev.off()

fc.stats(yhat = best.arimax$fitted , y = y.com)


# forecasting intervals ########################

y.3rd.exp.smoo <- hw(y.com.ts , h = 8 , initial = "simple", seasonal ="additive", alpha = 0.4, beta = 0.2, gamma = 0.3)

win.metafile("forecast_ts_for-err_3es.wmf", width=7, height = 7)
autoplot(forecast(y.3rd.exp.smoo))+ xlim(20,28) + ggtitle(expression(paste("Forecasts 3rd ES (", alpha==0.4,",",beta==0.2,",",gamma==0.3,")")))+ ylab("demand") +xlab("time/period")
dev.off()


y.ets <- ets(y.com.ts)

win.metafile("forecast_ts_for-err_optes.wmf", width=7, height = 7)
autoplot(forecast(y.ets, 8)) + xlim(20,28) + ggtitle(expression(paste("Forecasts optimal ES (", alpha==0,",",beta==0,",",gamma==0,")")))+ ylab("demand") +xlab("time/period") + ylim(15,120)
dev.off()

autoplot(forecast(ets(y.com.ts),28)) + xlim(20,33)

################################
# linear and quadratic trends
tmp.len <- 50
tmp.y.err <- rnorm(tmp.len, 0 , .5)
tmp.y.const <- rep(10, tmp.len)
tmp.y.trend <- 1:tmp.len*0.1 

tmp.y.lin <- tmp.y.err + tmp.y.const +tmp.y.trend

tmp.x <- (1:tmp.len)

tmp.y.trend <- (1:tmp.len)^2*0.0075 
tmp.y.quad <- tmp.y.err + tmp.y.const +tmp.y.trend 

win.metafile("forecast_ts_diff_comp.wmf", width=9, height = 6)
par(mfrow =c(2,3), family = "serif") 
par(mar = c(5,5,.1,.1))
plot.ts(tmp.y.lin, xlab= "time/periods", ylab ="demand", lwd = 2, cex.lab = 2, cex.axis = 1.5)
plot.ts(c(NA,diff(tmp.y.lin)), xlab= "time/periods", ylab =expression(Delta*y), lwd = 2, cex.lab = 2, cex.axis = 1.5)
abline(lm(c(NA,diff(tmp.y.lin)) ~ tmp.x), col="red", lwd=2)
plot.ts(c(NA,NA,diff(diff(tmp.y.lin))), xlab= "time/periods", ylab =expression(Delta^2*y), lwd = 2, cex.lab = 2, cex.axis = 1.5)
abline(lm(c(NA,NA,diff(diff(tmp.y.lin))) ~ tmp.x), col="red", lwd=2)

plot.ts(tmp.y.quad, xlab= "time/periods", ylab ="demand", lwd = 2, cex.lab = 2, cex.axis = 1.5)
plot.ts(c(NA,diff(tmp.y.quad)), xlab= "time/periods", ylab =expression(Delta*y), lwd = 2, cex.lab = 2, cex.axis = 1.5)
abline(lm(c(NA,diff(tmp.y.quad)) ~ tmp.x), col="red", lwd=2)
plot.ts(c(NA,NA,diff(diff(tmp.y.quad))), xlab= "time/periods", ylab =expression(Delta^2*y), lwd = 2, cex.lab = 2, cex.axis = 1.5)
abline(lm(c(NA,NA,diff(diff(tmp.y.quad))) ~ tmp.x), col="red", lwd=2)
dev.off()

#####################################################
# sporadic demand
y.spor.inter <- round(rexp(20, 6)*10)

y.spor.dem <- numeric(sum(y.spor.inter+1))

y.dem.set <- rpois(100,50)
y.spor.dem[1] <- y.dem.set[1]
tmp <- 1
for(i in 1:length(y.spor.inter)){
  tmp <- tmp + y.spor.inter[i] + 1
  y.spor.dem[tmp] <- y.dem.set[i]
  }

win.metafile("forecast_spor_dem.wmf", width=9, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
plot(1:length(y.spor.dem),y.spor.dem, type="h", xlab= "time/periods", ylab ="demand", lwd = 2 , cex.lab = 2, cex.axis = 1.5, pch =15)
dev.off()


win.metafile("forecast_spor_hist.wmf", width=6, height = 6)
par(family="serif", mar = c(5,5,.1,.1))
hist(y.spor.dem[which(y.spor.dem>0)], 10, main ="", xlab="demand",ylab = "", cex.lab = 2, cex.axis = 1.75, col="lightgrey")
dev.off()


dens.it <- table(factor(y.spor.inter, levels=0:7), useNA = "ifany" )/length(y.spor.inter)
which(y.spor.dem>0)

mu.dem <- mean(y.spor.dem[which(y.spor.dem>0)])
sd(y.spor.dem[which(y.spor.dem>0)])

cond.exp <- sapply(1:length(dens.it), function(x) sum(tail(dens.it, length(dens.it)-x+1)))

list.cond.exp <- lapply(1:length(cond.exp), function(x) as.numeric(tail(dens.it, length(dens.it)-x+1) /cond.exp[x] ))

win.metafile("forecast_spor_pred.wmf", width=7.5, height = 6.5)
par(mar = c(4.25,4.25,1.5,.1), family="serif")
plot(tail(1:length(y.spor.dem),10), tail(y.spor.dem,10), type="h", xlab= "time/periods", ylab ="demand", lwd = 2 , cex.lab = 1.75, cex.axis = 1.25, pch =15, xlim=c(56,length(y.spor.dem)+9), ylim=c(0,65), main =expression(paste("expected demand depending on ",j," and ",h," at ",t==66,",...,",73)), cex.main = 1.75)
for(i in 1:length(list.cond.exp) ){
  lines(0.1*i+(length(y.spor.dem)+i):(length(list.cond.exp[[i]])+length(y.spor.dem)+(i-1)), list.cond.exp[[i]]*mu.dem, type="h", lwd=3, col=i+1)
}
legend("top", col=2:(length(list.cond.exp)+1), lwd =3 , legend =  paste("j=",0:7, sep="") , horiz = T)
dev.off()

#####################################################
# real world example ########################
# plot.ts(hyndsight)
# 
# 
# fc <- hw(subset(hyndsight,end=length(hyndsight)-35),
#          damped = TRUE, seasonal="multiplicative", h=35)
# 
# autoplot(forecast(auto.arima(hyndsight))) 
# checkresiduals(auto.arima(log(hyndsight)))
# 
# tmp <- auto.arima(log(hyndsight))
# tmp.res <- tmp$residuals
# ub <- 1.97*sd(tmp.res)
# lb <- -1.97*sd(tmp.res)
# plot.ts(tmp.res)
# abline(h = c(ub,lb), col="red",lty =2, lwd=2)
