# standard normal distribution plot  with 3 sigma interval
lims <- c(-5,5)
par(family="serif", mar=c(3,5,1,1), bg="white")
curve(dnorm, lims[1], lims[2],  yaxs="i", ylim=c(0,.45), lwd=2, ylab="density", xaxt="n", xlab="x", cex.lab=1.5, cex.axis=1.25, n=1000)
axis(1, labels=c(expression(mu-4*sigma),expression(mu-3*sigma),expression(mu-2*sigma),expression(mu-1*sigma),expression(mu),expression(mu+1*sigma),expression(mu+2*sigma),expression(mu+3*sigma),expression(mu+4*sigma)), at=c(-4:4), cex.axis=1.25)

# right area
cord.x <- c(3,seq(3,lims[2], length.out = 100),lims[2]) 
cord.y <- c(0,dnorm(seq(3,lims[2], length.out = 100)),0) 
polygon(cord.x,cord.y,col='grey', border="black")
# central area
cord.x <- c(lims[1],seq(lims[1],-3, length.out = 100), -3) 
cord.y <- c(0,dnorm(seq(lims[1],-3, length.out = 100)), 0) 
polygon(cord.x,cord.y,col='grey', border="black")
# left area
cord.x <- c(-3,seq(-3,3, length.out = 1000),3) 
cord.y <- c(0,dnorm(seq(-3,3, length.out = 1000)),0) 
polygon(cord.x,cord.y,col='grey', border="black", density=12, angle=45)
# annotations
abline(v=c(-3,3), lwd=1.75, lty=2, col="black")
abline(v=0, lwd=1.75, lty=3, col="black")
text(0,.115,round(pnorm(3)-pnorm(-3),4), cex=1.25)
text(4.25,.1,round(1-pnorm(3),4), cex=1.25)
lines(rbind(c(4.25,.08),c(3.2,.003)), lwd=2, col="grey")
text(-4.25,.1,round(1-pnorm(3),4), cex=1.25)
lines(rbind(c(-4.25,.07),c(-3.2,.003)), lwd=2, col="grey")


############################
# Newsvendor model
############################

##########################
# profit optimization

# Profit function 
obj.news.prof <- function(x, mu, sigma, p, c.unit ){
  tmp1 <- function(y) y*dnorm(y, mean=mu, sd = sigma)
  p * integrate(tmp1, lower = 0, upper=x)$value + p * x * (1 - pnorm(x , mean=mu, sd = sigma) ) - c.unit * x
} 

# optimal q (profit)
q.opt.prof <- mu + sigma * qnorm( (p-c.unit)/p )

# Example
p <- 1
c.unit <- .25
mu <- 50
sigma <- 10

# profit values 
x.val <- seq(45,65, length.out = 100)
y.val <- sapply(x.val, function(x) obj.news.prof(x, mu = mu, sigma = sigma, p = p , c.unit = c.unit) )

# plot of profit
plot(x.val, y.val , type="l")

val.mat <- sapply(seq(0.00001, 0.99999, length.out = 100), function(x){
  q.opt.prof <- mu + sigma * qnorm( x )
  c(x, q.opt.prof, obj.news.prof(q.opt.prof, mu=mu, sigma = sigma, p=1, c.u = 1-x))
})

par(family="serif", mar = c(6,5,4,4.5), bg ="white")
plot(val.mat[1,], val.mat[2,] , type="l", lwd=3, xlab="", ylab = expression(paste("opt. order quantity ", q^opt)) , main = expression(paste("normal distribution with ", list(mu==50,sigma==10))), cex.main = 1.75, cex.axis = 1.25, cex.lab = 1.75)
mtext(side = 1, expression(frac(p-c^unit, p)), line = 4.5, cex = 1.75)
par(new = T)
plot(val.mat[1,], val.mat[3,] , type="l", xaxt = "n", yaxt = "n", xlab="", ylab ="", col="red", lwd=3)
axis(side=4, cex.axis = 1.25)
mtext(side = 4, "optimal exp. profit", line = 3, cex = 1.75)
legend("topleft", lwd=3, col=c("black","red"), legend = c("opt. q","opt. profit") , cex=1.75 )


############################################################
# cost optimization
############################################################

mu <- 100
sigma <- 20
co = .25
cu = 1

# the direct way
zf.news <- function(x, cu, co, mu, sigma){
  tmp1 <- function(y) co*(x-y)*dnorm(y, mean=mu, sd = sigma)
  tmp2 <- function(y) cu*(y-x)*dnorm(y, mean=mu, sd = sigma)
  integrate(tmp1, lower = 0, upper=x)$value + integrate(tmp2, lower = x, upper=Inf)$value
} 


# the elegant way 
zf.news.elegant <- function(x, cu, co, mu, sigma){
  co*(x - mu) + (cu + co)*sigma*(dnorm((x-mu)/sigma) - ((x-mu)/sigma) * (1- pnorm((x-mu)/sigma)))
}

# find optimum
opt.news <- optim(par = 100, fn=zf.news, lower = 0, upper = 200, method="L-BFGS-B", co = co, cu = cu, sigma = sigma, mu=mu)
opt.news.elegant <- optim(par = 100, fn=zf.news.elegant, lower = 0, upper = 200, method="L-BFGS-B", co = co, cu = cu, sigma = sigma, mu=mu)

# theoretical optimum
cr <- cu/(cu+co) # critical ratio

q.star <- qnorm(cr, mean=mu, sd = sigma) 	# direct
q.star <- mu + sigma * qnorm(cr) 			    # transformation

# difference numerical vs. theoretical optimum
opt.news$par - q.star

########################################################
# Plot cost function
########################################################

# Example 1 ############################

p <- 1
c.unit <- .25
mu <- 50
sigma <- 10
co <- c.unit
cu <- p - c.unit

q.opt.cost <- mu + sigma * qnorm(cu/(cu+co)) 			# Transformation
q.opt.prof <- mu + sigma * qnorm( (p-c.unit)/p )

n <- 1000 		# nb. points
x.low <- 45		# x-Interval start
x.high<- 70	  # x-Interval end
# x vector 
x.vals <- seq(x.low, x.high , length.out = n)
# y vector (cost)
y.vals.cost <- sapply(x.vals, function(x) zf.news(x, co = co, cu = cu, sigma = sigma, mu=mu)) 
# y vector (profit)
y.val.prof <- sapply(x.vals, function(x) obj.news.prof(x, mu = mu, sigma = sigma, p = p , c.unit = c.unit) )

par(family="serif", mar = c(6,5,4,4.5), bg ="white")
plot(x.vals, y.vals.cost, type="l", xlab="order quantity q", ylab="Total expected cost", lwd=3 , cex.main = 1.75, cex.axis = 1.25, cex.lab = 1.75, main=expression(paste(mu==50,", ",sigma==10,", ",p==1,", ",c^unit==0.25)), ylim=c(3,7.5))
# add opt. cost
abline(h=zf.news.elegant(q.opt.cost, co = co, cu = cu, sigma = sigma, mu=mu), lty=2, col="red", lwd=1.5)
abline(v=q.opt.cost, lty=2, col="red", lwd=1.5)
# add profit function
par(new =T)
plot(x.vals , y.val.prof ,type="l", xaxt="n", yaxt="n", col="blue", xlab="", ylab="", lwd=3)
# add opt. profit
abline(h=obj.news.prof(q.opt.cost, mu = mu, sigma = sigma, p = p , c.unit = c.unit), lty=2, col="blue")
axis(side=4, cex.axis = 1.25)
mtext(side = 4, "total expected profit", line = 3, cex = 1.75)
legend("topleft", lwd=3, col=c("black","blue"), legend = c("expected cost","expected profit") , cex=1.5 )


########################################################
# Fixing service levels
########################################################
mu=0
sigma=1

# alpha-SL ############################
alpha = .95
q.alpha <- mu + sigma * qnorm(alpha)


zf.news(q.alpha) 						        # cost for alpha=.95
mu*cu - zf.news.elegant(q.alpha)		# profit for alpha=.95


# beta-SL ############################
# first order loss function
v.fun <- function(x, mu = 0, sigma = 1) {
  sigma * (dnorm((x-mu)/sigma) - (x-mu)/sigma * (1-pnorm((x-mu)/sigma) ))
}
# Plot standard los function
par(family ="serif", mar=c(4,5,0.1,0.1))
curve(-3, 3, expr=v.fun, xlab = expression(paste("normalized order quantity ", q^"'")), ylab =expression(paste("expected loss", L(Z,q^"'"))), cex.lab=1.75, cex.axis=1.5, lwd=3)
abline(v=0, lty=5)
abline(h=0, lty=5)

# optimizer for Beta
v.fun.opt <- function(x) (beta - (1 - v.fun(x)/mu))^2 
# optimize q for given Beta
q.beta <- optim(fn = v.fun.opt, par = mu, lower = 0, upper = 200, method="L-BFGS-B")$par
# associated cost
zf.news(q.beta) 				
# associated profit 
mu*cu - zf.news(q.beta)


# Example 1 ############################
# code as before
p <- 1
c.unit <- .25
mu <- 50
sigma <- 10
co <- c.unit
cu <- p - c.unit

q.opt.cost <- mu + sigma * qnorm(cu/(cu+co)) 			# Transformation
q.opt.prof <- mu + sigma * qnorm( (p-c.unit)/p )

n <- 1000 		
x.low <- 45		
x.high<- 70	

x.vals <- seq(x.low, x.high , length.out = n) 
y.val.prof <- sapply(x.vals, function(x) obj.news.prof(x, mu = mu, sigma = sigma, p = p , c.unit = c.unit) )
y.val.alpha <- sapply(x.vals, function(x) pnorm(x, mean=mu, sd = sigma) )
y.val.beta <- sapply(x.vals, function(x) 1- v.fun(x, mu = mu, sigma = sigma)/mu )

par(family="serif", mar = c(6,5,4,4.5), bg ="white")
plot(x.vals, y.val.prof, type="l", xlab="order quantity q", ylab="Total expected profit", lwd=3 , cex.main = 1.75, cex.axis = 1.25, cex.lab = 1.75, main=expression(paste(mu==50,", ",sigma==10,", ",p==1,", ",c^unit==0.25)), ylim=c(31.5,35))

abline(h=mu*cu - zf.news.elegant(q.opt.cost, co = co, cu = cu, sigma = sigma, mu=mu), lty=2, col="red", lwd=1.5)
abline(v=q.opt.cost, lty=2, col="red", lwd=1.5)
# add service levels
par(new =T)
plot(x.vals , y.val.beta ,type="l", xaxt="n", yaxt="n", col="blue", xlab="", ylab="", lwd=3, ylim = c(.4,1))
lines(x.vals , y.val.alpha , col="red", lwd=3)

abline(h=obj.news.prof(q.opt.cost, mu = mu, sigma = sigma, p = p , c.unit = c.unit), lty=2, col="blue")
axis(side=4, cex.axis = 1.25)
mtext(side = 4,expression(paste(alpha/beta, "-service level")) , line = 3, cex = 1.75)
legend("bottomright", lwd=3, col=c("black","red","blue"), legend = c("exp. profit", "Alpha SL","Beta SL") , cex=1.5 )


###############################################
# Discrete Newsvendor 
###############################################

# PLot Poisson distributions

lambda.vec <- c(1,4,7,10)
x.vals <- 0:20

par(family="serif", mar = c(4.5,4.5,.1,.1), bg ="white")
plot(x.vals, dpois(x.vals, lambda.vec[1]), type="h", xlab="k", ylab="density", lwd=3 , cex.main = 1.75, cex.axis = 1.25, cex.lab = 1.75)
points(x.vals, dpois(x.vals, lambda.vec[1]), pch = 15, col="black")
for(i in 2:length(lambda.vec)){
  lines(x.vals+(i-1)*.15,  dpois(x.vals, lambda.vec[i]) , col=i, type="h", lwd = 3)
  points(x.vals+(i-1)*.15, dpois(x.vals, lambda.vec[i]), pch = 15+i-1, col=i)
}
legend("topright", lwd=3, col=1:length(lambda.vec), legend = parse(text=paste("lambda==", lambda.vec )) , cex=1.5 )


# Example 1 ####################################################

lambda <- 4
p <- 1
c.unit <- .25
co <- c.unit
cu <- p - c.unit

x.vals <- 0:13
exa.1.pois <- data.frame(p=x.vals, 
                         dens = dpois(x.vals, lambda),
                         os = NA, us = NA, tc = NA)
for(i in 1:nrow(exa.1.pois)){
  tmp.p <- exa.1.pois$p[i]
  tmp.os <- sum((tmp.p - exa.1.pois$p[1:(i-1)])*exa.1.pois$dens[1:(i-1)])
  tmp.us <- sum((exa.1.pois$p[i:nrow(exa.1.pois)] - tmp.p)*exa.1.pois$dens[i:nrow(exa.1.pois)])
  exa.1.pois$os[i] <- tmp.os
  exa.1.pois$us[i] <- tmp.us
  exa.1.pois$tc[i] <- tmp.os*(cu+co) + cu*(lambda - tmp.p  )  #tmp.us*cu + tmp.os*co
}

round(exa.1.pois,3)

exa.1.pois <- exa.1.pois %>% add_column(beta = 1-exa.1.pois$us/lambda , alpha = cumsum(exa.1.pois$dens) )

# Example 2 ##################################################

lambda <- 2
p <- 1
c.unit <- .75
co <- c.unit
cu <- p - c.unit

x.vals <- 0:13

exa.2.pois <- data.frame(p=x.vals, 
                         dens = dpois(x.vals, lambda),
                         os = NA, us = NA, tc = NA)
for(i in 1:nrow(exa.2.pois)){
  tmp.p <- exa.2.pois$p[i]
  tmp.os <- sum((tmp.p - exa.2.pois$p[1:(i-1)])*exa.2.pois$dens[1:(i-1)])
  tmp.us <- sum((exa.2.pois$p[i:nrow(exa.2.pois)] - tmp.p)*exa.2.pois$dens[i:nrow(exa.2.pois)])
  exa.2.pois$os[i] <- tmp.os
  exa.2.pois$us[i] <- tmp.us
  exa.2.pois$tc[i] <- tmp.os*(cu+co) + cu*(lambda - tmp.p  )  #tmp.us*cu + tmp.os*co
}

round(exa.2.pois,3)

exa.2.pois <- exa.2.pois %>% add_column(beta = 1-exa.2.pois$us/lambda , alpha = cumsum(exa.2.pois$dens) )


# Plot results ###########################
library(tidyverse)
library(reshape2)

ggplot(exa.1.pois, aes(x = p)) +
  geom_col(aes( y = tc*10, fill="redfill")) +
  ylab("total expected cost") + xlab("order quantity q")+
  geom_text(aes(y = tc*10, label = round(tc, 1)), fontface = "bold", vjust = 1.4, color = "black", size = 4) +
  geom_line(aes(y = alpha*100, group = 1, color = 'black'), size = 1.25) +
  geom_text(aes(y = alpha * 100, label = round(alpha * 100,1)), vjust = 2, color = "black", size = 4) +
  geom_line(aes(y = beta*100, group = 1, color = 'blue'), size = 1.25) +
  geom_text(aes(y = beta * 100+4, label = round(beta * 100,1)), vjust = 1.4, color = "blue", size = 4) +
  scale_y_continuous(labels = function(x) x/10, sec.axis = sec_axis(name = "service level (%)", trans = ~ .  )) +
  scale_fill_manual('', labels = 'Total expected cost', values = "#C00000") +
  scale_color_manual('', labels = c("Alpha SL","Beta SL"), values = c('black',"blue") ) +
  theme_minimal() + theme(legend.position = c(0.1, .8)) #
dev.off()


ggplot(exa.2.pois, aes(x = p)) +
  geom_col(aes( y = tc*10, fill="redfill")) +
  xlim(-.5,8.5) + ylab("total expected cost") + xlab("order quantity q")+
  geom_text(aes(y = tc*10, label = round(tc, 1)), fontface = "bold", vjust = 1.4, color = "black", size = 4) +
  geom_line(aes(y = alpha*100, group = 1, color = 'black'), size = 1.25) +
  geom_text(aes(y = alpha * 100, label = round(alpha * 100,1)), vjust = 2, color = "black", size = 4) +
  geom_line(aes(y = beta*100, group = 1, color = 'blue'), size = 1.25) +
  geom_text(aes(y = beta * 100+4, label = round(beta * 100,1)), vjust = 1.4, color = "blue", size = 4) +
  scale_y_continuous(labels = function(x) x/10, sec.axis = sec_axis(name = "service level (%)", trans = ~ .  )) +
  scale_fill_manual('', labels = 'Total expected cost', values = "#C00000") +
  scale_color_manual('', labels = c("Alpha SL","Beta SL"), values = c('black',"blue") ) +
  theme_minimal() + theme(legend.position = c(0.1, .8)) #



