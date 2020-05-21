

ts.lager.opt <- function(par, demand, wwbz, llini, obj = "cost", co = 0, cu = 0, B = 0, alpha.z = 95, beta.z = 99, tt){
  
  sim.R <- tS.lager(d = demand, t = tt , S = par, l.ini = llini, n = length(demand), wbz = wwbz)
  
  tmp <- sim.R[-c(1),]
  tmp.l <- tmp[,"Lagerbestand"]
  tmp.id <- 1:length(tmp.l)
  tmp.id.sub <- tmp[,"Lieferung"] > 0
  tmp.id.sub.w <- which(tmp.id.sub) - 1
  tmp.id.sub.w <-tmp.id.sub.w[tmp.id.sub.w >0]
  
  # Alpha-Servicegrad
  sim.alpha <- (100*(1 - sum( tmp.l[tmp.id.sub.w] < 0)/length(tmp.id.sub.w)))
  sim.alpha.diff <- (alpha.z - sim.alpha)^2
  # Beta-Servicegrad
  sim.beta <- (100*(1 - sum(-tmp.l[ tmp.l < 0])/sum(tmp[,"Bedarf"])))
  sim.beta.diff <- (beta.z - sim.beta)^2
  
  # Kosten
  sim.cost <- ( sum(tmp[,"Bestellung"] > 0)*B +  sum(tmp.l[ tmp.l > 0]) * co + sum(-tmp.l[ tmp.l < 0]) * cu)/length(tmp.l)
  
  #list(sim.res = sim.R, cost = sim.cost, alpha = sim.alpha, beta = sim.beta, res.dat = tmp, res.lag = tmp.l)
  return(ifelse(obj == "cost", sim.cost, ifelse(obj == "alpha", sim.alpha.diff, sim.beta.diff)) )
  
}


v.fun <- function(x, mu3, sigma3) {
  sigma3 * (dnorm((x-mu3)/sigma3) - (x-mu3)/sigma3 * (1-pnorm((x-mu3)/sigma3) ))
}

v.fun.opt <- function(x, mu2, sigma2) {
  (97/100 - (1 - v.fun(x, mu3 = mu2, sigma3 = sigma2)/mu2))^2
} # Abweichung von Zielwert

# beta-Servicegrad
v.fun <- function(x, mu = 0, sigma = 1) {
  sigma * (dnorm((x-mu)/sigma) - (x-mu)/sigma * (1-pnorm((x-mu)/sigma) ))
}

v.fun.opt.loss <- function(x, loss = loss.norm, ...) (loss - v.fun(x, ...))^2 # Abweichung von Zielwert

v.fun.opt.loss(x=0.11, loss=.03*20*8/sqrt(10)/4)

optim(fn = v.fun.opt.loss, par = 1, lower = -5, upper = 5, method="L-BFGS-B", loss = .03*20*8/sqrt(10)/4)$par


######################################
# (t,S) inventory with lead time     #
######################################

tS.inv.L <- function(t = 1, S = 100, l.ini, dem, L = 2){
  # n ... number of periods
  n <- length(dem) 
  # t ... order interval
  # S ... order level
  # l.ini ... initial stock
  # dem ... demand
  
  tmp.l.mat <- matrix(0, ncol = 6, nrow = n+1+L) # initialize result matrix
  colnames(tmp.l.mat) <- c("d","i","o", "ip", "x","del")
  rownames(tmp.l.mat) <- -L:n # set period names
  names(dem) <- 1:n
  tmp.l.mat["0",c("i", "ip")] <- l.ini		# set initial stock
  tmp.l.mat[as.character(1:n),"d"] <- dem # initialize demand
  
  # calculate stocks
  for(i in 1:n){
    tmp.l.mat[as.character(i), "i"] <- tmp.l.mat[as.character(i-1), "i"] - dem[as.character(i)] + tmp.l.mat[as.character(i), "del"]
    tmp.l.mat[as.character(i), "ip"] <- tmp.l.mat[as.character(i-1),"ip"] - dem[as.character(i)] 
    tmp.l.mat[as.character(i), "o"] <- tmp.l.mat[as.character(i-1), "o"] - tmp.l.mat[as.character(i), "del"]
    
    if(i %% t == 0){
      tmp.l.mat[as.character(i), "x"] <- S - tmp.l.mat[as.character(i) ,"ip"]
      tmp.l.mat[as.character(i), "ip"] <- tmp.l.mat[as.character(i), "ip"] + tmp.l.mat[as.character(i), "x"]
      tmp.l.mat[as.character(i), "o"] <- tmp.l.mat[as.character(i), "o"] + tmp.l.mat[as.character(i), "x"]
      if(i + L + 1 <= n){ 
        tmp.l.mat[as.character(i + L + 1), "del"] <- tmp.l.mat[as.character(i + L + 1), "del"] + tmp.l.mat[as.character(i), "x"]
      }
    }
  }
  
  return(as.data.frame(tmp.l.mat))
}

# L < T ####################################
# number periods
n.1 <- 20
# demand
dem.1 <- rpois(n.1, 20)
# simulate inventory
sim.1 <- tS.inv.L(l.ini = 120, dem = dem.1, t=4, L = 2, S = 120 )
# prepare x axis
x.vec <- 1:(n.1+1)
o.vec <- x.vec[ x.vec %% 4 == 0]

par(family="serif", mar = c(4.5,4.5,3,.1), bg ="white")
plot(tail(sim.1$i,n.1+1), type ="s", xaxt="n", xlab="periods",ylab="stock", xlim=c(4, n.1), lwd=2, cex.lab=1.5, cex.axis=1.75, cex.main = 1.75, main = parse(text = "list(T==4,L==2,S==120,mu==20)"), ylim=c(min(sim.1$i),120))
lines(1:(n.1+1) + .2 , tail(sim.1$ip, n.1+1), type="s", lty=2, lwd=2)
abline(h=0)
abline(h=120)
axis(side = 1, at = 1+0:(n.1/4)*4, labels = parse(text = paste(0:(n.1/4),"%.%T", sep="")), cex.axis=1.75)
axis(side = 1, at = 1+0:(n.1/4)*4+2, labels = parse(text=paste(0:(n.1/4),"%.%T+L", sep="")) , tck = -.005, cex=.75, mgp=c(3,.5,0), cex.axis=1.5)
for(i in 1:length(o.vec)){
  segments(x0= o.vec[i] + 1.2, y0 = sim.1$i[which(sim.1$x > 0)[i]], x1 = o.vec[i] + 1.2, y1 =sim.1$i[which(sim.1$x > 0)[i]] + sim.1$x[which(sim.1$x > 0)[i]], col="blue", lwd = 3)
}

legend("bottom", horiz = T, col=c("black","black","blue"), lwd=c(2,2,3), lty=c(1,2,1), legend = c("phys. stock", "inv. pos.","q"), bty="n", bg=NULL, cex=1.5)

# L < T ####################################
# number perios
n.1 <- 40
# demand
dem.1 <- rpois(n.1, 20)
wbzz <- 6
tt <- 3
# simulate inventory
sim.1 <- tS.inv.L( l.ini = 170, dem = dem.1, t=tt, L = wbzz, S = 180 )
x.vec <- 1:(n.1+1)
o.vec <- x.vec[ x.vec %% 4 == 0]


par(family="serif", mar = c(4.5,4.5,3,.1), bg ="white")
S <- 180
sim.sub <- tail(sim.1, tt*6)
per.vec <- as.numeric(rownames(sim.sub))
o.id.vec <- per.vec %% 3 == 0

plot(per.vec, sim.sub$i, type ="s", xaxt="n", xlab="periods",ylab="stock", lwd=2, cex.lab=1.5, cex.axis=1.75, ylim = c(-40,180), main = parse(text = "list(T==3,L==6,S==180,mu==20)"), cex.main = 1.75, xlim= range(per.vec))
lines(per.vec+.2, sim.sub$ip, type="s", lty=2, lwd=2)
abline(h=0)
abline(h=180)
axis(side = 1, at = per.vec[o.id.vec], labels = parse(text = paste(1:sum(o.id.vec),"%.%T", sep="")), cex.axis=1.75)
legend("bottom", horiz = T, lwd=c(2,2), lty=c(1,2), legend = c("on-hand stock", "inventory position"), bty="n", bg=NULL, cex=1.75)

o.id <- per.vec[o.id.vec]+.2
ord.vec <- sim.sub$x[o.id.vec]
del.vec <- sim.sub$del[sim.sub$del > 0]
ord.vec <- c(head(del.vec, 2),ord.vec)
for(j in 1:length(o.id)){
  segments(x0=o.id[j], x1=o.id[j], y0=S, y1 = S - ord.vec[j+2], col=3+j, lwd =3  )
  text(x=o.id[j]+.5 , y= S - (ord.vec[j+2])/2, label=parse(text=paste("q[",j,"*T]",sep="")))
  segments(x0=o.id[j], x1=o.id[j], y0 = S - ord.vec[j+2], y1 = S - ord.vec[j+2] - ord.vec[j+1] , col=3+j-1, lwd =3  )
  text(x=o.id[j]+.5 , y= S - ord.vec[j+2] - ord.vec[j+1]/2, label=parse(text=paste("q[",j-1,"*T]",sep="")))
  segments(x0=o.id[j], x1=o.id[j], y0 = S - ord.vec[j+2] - ord.vec[j+1], y1 = S - ord.vec[j+2] - ord.vec[j+1] - ord.vec[j] , col=3+j-2, lwd =3  )
  text(x=o.id[j]+.5 , y=  S - ord.vec[j+2] - ord.vec[j+1] - ord.vec[j]/2, label=parse(text=paste("q[",j-2,"*T]",sep="")))
}

###############################
# Service level definition
###############################
mu <- 20
sd <- 4

oi <- 8
lead <- 4

n.1 <- 100+oi+lead
set.seed(7654321)
dem.1 <- round(rnorm(n.1, mu, sd))

rp <- lead + oi
S <- round(rp*mu + sqrt(rp)*sd*qnorm(.30))

sim.1 <- tS.inv.L( l.ini = 50, dem = dem.1, t = oi, L = lead, S = S )

par(family="serif", mar = c(4.5,4.5,3,.1), bg ="white")
tmp.plot <- tail(sim.1, n.1-lead-oi)
plot(tmp.plot$i, type ="s",  xlab="periods",ylab="inventory level", lwd=2, cex.lab=1.5, cex.axis=1.75, cex.main = 1.75, main = parse(text = "list(T==8,L==4,S==233,n==100,y%~%N(20,4^2))"))
abline(h=0)

# inventory measures
# discard starting periods
tmp.plot <- tail(sim.1, n.1-lead-oi)

# number of periods with stock-outs
sum(tmp.plot$i < 0)
# number order intervals with stock outs
oi.end <- tmp.plot$i[ which(tmp.plot$del > 0) - 1 ]
sum(oi.end < 0, na.rm =T)
# sum of backorders
sum(oi.end[oi.end < 0], na.rm =T)
# sum of backlogs
sum(tmp.plot$i[tmp.plot$i < 0])


# Cost-optimal order interval

c.sh <- 1
c.bo <- 3
c.or <- 50
mu <- 20
sd <- 4

k.star <- qnorm(c.bo/(c.bo+c.sh))
S.star <- mu + sd * k.star

tot.cost.func <- function(t) c.sh*mu*t^2 - (c.sh +c.bo)*dnorm(k.star)*sd*t^.5 - 2*c.or

cost.tS <- function(t, L=4, sd, mu , c.or, c.sh, c.bo ){
  
  k.star <- qnorm(c.bo/(c.bo+c.sh))
  S.star <- mu + sd * k.star
  
  tmp.end <- (c.sh+c.bo)*dnorm(k.star)*sd*sqrt(L+t)/t
  tmp.stock <- t*mu*c.sh/2
  tmp.order <- c.or/t
  
  return(list(S.star = S.star,  k.star =k.star, stock.cost = tmp.stock, order.cost = tmp.order, end.cost = tmp.end, total.cost = tmp.order + tmp.end + tmp.stock))
} 

t.seq <- seq(1,20,length.out = 100)
res.TS.cost <- cost.tS(t=t.seq, sd =4, mu =20, c.or = 300 , c.sh = .5, c.bo = 1)

par(family="serif", mar = c(4.5,4.5,3,.1), bg ="white")
plot(t.seq,res.TS.cost$total.cost, type="l", xlab="order interval T", ylab="cost per period", ylim=c(0,max(res.TS.cost$total.cost)), lwd=2, cex.lab = 1.75, cex.axis=1.5 , cex.main =1.75, main = parse(text = "list(c^or==300,c^sh==0.5,c^bo==1,mu==20,sigma==4,L==4)"))
lines(t.seq,res.TS.cost$stock.cost, col="green", lwd=2)
lines(t.seq, res.TS.cost$order.cost, col="blue", lwd=2)
lines(t.seq, res.TS.cost$end.cost, col="red", lwd=2)
legend("top" , horiz = T, lwd=2, col=c("black","green","blue","red"), legend=c("total","stock","order","end-of-oi"), title = "cost types", bty="n", bg=NULL, cex=1.5)

# simulate inventory systems

TS.sim.func <- function(n = 500, reps = 70, mu = 20, sd = 4, c.sh = .5, c.bo = 1, c.or = 300, oi.range = 1:12, lead.range = 0:6 ){
  k.star <- qnorm(c.bo/(c.bo + c.sh))
  exp.des <- expand.grid(oi.range, lead.range, 1:reps)
  colnames(exp.des) <- c("order.intervall","lead.time","replication")
  res.vec <- numeric(nrow(exp.des))
  
  for(i in 1:nrow(exp.des)){
    print(paste(i,"|",nrow(exp.des)))
    
    tmp.oi <- exp.des[i,1]
    tmp.lead <- exp.des[i,2]
    tmp.rp <- tmp.lead + tmp.oi
    tmp.n <- n + tmp.oi + tmp.lead
    tmp.dem <- rnorm(tmp.n, mu, sd)
    
    S <- tmp.rp * mu + sqrt(tmp.rp) * sd * qnorm(k.star)
    
    tmp.sim <- tS.inv.L(l.ini = 50, dem = tmp.dem, t = tmp.oi, L = tmp.lead, S = S )
    tmp.plot <- tail(tmp.sim, tmp.n - tmp.lead - tmp.oi)
    
    res.vec[i] <- (sum(tmp.plot$i[tmp.plot$i > 0]*c.sh) - sum(tmp.plot$i[tmp.plot$i < 0]*c.bo) + sum(tmp.plot$x > 0)*c.or)/nrow(tmp.plot)
    
  }
  return(data.frame(exp.des, total.cost = res.vec ))
}


res.sim <- TS.sim.func()


library(tidyverse)

par(family="serif", mar = c(4.5,4.5,3,.1), bg ="white")
ggplot(res.sim, aes(x = as.factor(lead.time), y = total.cost)) + 
  geom_boxplot() + ylab("Total cost") + xlab("Lead time") + #title("order interval") +
  facet_wrap(~ as.factor(order.intervall), scale="free")


###################################################
# stochastic lead times
###################################################

# calculate performance measures
stock.perf <- function(res, c.bo = 1 , c.sh = .5, c.or = 0, cutoff = 10){
  n <- nrow(res)
  tmp.plot <- tail(res, n - cutoff)
  
  tot.cost <- (sum(tmp.plot$i[tmp.plot$i > 0]*c.sh) - sum(tmp.plot$i[tmp.plot$i < 0]*c.bo) + sum(tmp.plot$x > 0)*c.or)/nrow(tmp.plot)
  nb.oi <- sum(tmp.plot$x > 0)
  id.oi <- which(tmp.plot$del > 0) - 1
  alpha <- sum(tmp.plot$i[id.oi] >= 0)/nb.oi
  id.oi.beta <- id.oi[tmp.plot$i[id.oi] < 0]
  #beta <- 1 + (sum(tmp.plot$i[id.oi.beta])/nb.oi)/mean(tmp.plot$d)
  beta <- 1 + sum(tmp.plot$i[id.oi.beta])/sum(tmp.plot$d)
  return(list(cost = tot.cost, alpha = alpha, beta = beta))
}

tS.inv.L.pois <- function(n = 100, t = 1, S = 100, l.ini, dem, L.lambda = 2){
  # n ... nb. periods
  # t ... order interval
  # S ... order level
  # l.ini ... initial stock
  # dem ... demand vector
  # L.lambda ... Poisson parameter/exp. lead time
  n <- length(dem)
  tmp.l.mat <- matrix(0, ncol = 6, nrow = n+1+L.lambda) # initialize inventory
  colnames(tmp.l.mat) <- c("d","i","o", "ip", "x","del")
  # set period names
  rownames(tmp.l.mat) <- -L.lambda:n 
  # initial stocks
  tmp.l.mat["0",c("i", "ip")] <- l.ini		
  names(dem) <- 1:n
  tmp.l.mat[as.character(1:n),"d"] <- dem
  
  
  # calculate stocks
  for(i in 1:n){
    tmp.l.mat[as.character(i), "i"] <- tmp.l.mat[as.character(i-1), "i"] - dem[as.character(i)] + tmp.l.mat[as.character(i), "del"]
    tmp.l.mat[as.character(i), "ip"] <- tmp.l.mat[as.character(i-1),"ip"] - dem[as.character(i)] 
    tmp.l.mat[as.character(i), "o"] <- tmp.l.mat[as.character(i-1), "o"] - tmp.l.mat[as.character(i), "del"]
    
    if(i %% t == 0){
      tmp.l.mat[as.character(i), "x"] <- S - tmp.l.mat[as.character(i) ,"ip"]
      tmp.l.mat[as.character(i), "ip"] <- tmp.l.mat[as.character(i), "ip"] + tmp.l.mat[as.character(i), "x"]
      tmp.l.mat[as.character(i), "o"] <- tmp.l.mat[as.character(i), "o"] + tmp.l.mat[as.character(i), "x"]
      tmp.lead <- rpois(1, L.lambda)
      
      if(i + tmp.lead + 1 <= n){ 
        tmp.l.mat[as.character(i + tmp.lead + 1), "del"] <- tmp.l.mat[as.character(i + tmp.lead + 1), "del"] + tmp.l.mat[as.character(i), "x"]
      }
    }
  }
  
  return(as.data.frame(tmp.l.mat))
}


mu <- 20
sd <- 4
lambda <- 2
c.sh <- .5
c.bo <- 1
oi <- 8
rp <- lambda + oi

mu.rp <- rp * mu
sd.rp <- sqrt(lambda*sd^2 + mu^2*lambda^2+oi*sd^2)

S.cost <- mu.rp + sd.rp * qnorm(c.bo/(c.bo+c.sh))
S.cost.fix <- mu.rp + sd * sqrt(rp) * qnorm(c.bo/(c.bo+c.sh))
S.alpha <- mu.rp + sd.rp * qnorm(0.95)
S.alpha.fix <- mu.rp + sd * sqrt(rp) * qnorm(0.95)



# beta-Servicelevel
v.fun <- function(x, mu = 0, sigma = 1) {
  sigma * (dnorm((x-mu)/sigma) - (x-mu)/sigma * (1-pnorm((x-mu)/sigma) ))
}
# optimizer beta SL
v.fun.opt.loss <- function(x, loss = loss.norm, ...) (loss - v.fun(x, ...))^2 


S.beta <- mu.rp + sd.rp * optim(fn = v.fun.opt.loss, par = 1, lower = -5, upper = 5, method="L-BFGS-B", loss = .03*20*8/sd.rp)$par
S.beta.fix <- mu.rp + sqrt(rp)*sd * optim(fn = v.fun.opt.loss, par = 1, lower = -5, upper = 5, method="L-BFGS-B", loss = .03*20*8/sqrt(10)/4)$par

# simulation results ############
# stochastic vs fixed lead time

n <- 10000 + rp
dem <-  rnorm(n, mu, sd)
# S cost optimal 
res.cost.stoch <- tS.inv.L.pois( t = 8, S = S.cost, l.ini = 200 , dem = dem, L.lambda = lambda)
res.cost.stoch.fix <- tS.inv.L.pois( t = 8, S = S.cost.fix, l.ini = 200 , dem = dem, L.lambda = lambda)
res.cost.fix <- tS.inv.L(t = 8, S = S.cost.fix, l.ini = 200 , dem = dem, L = lambda)

stock.perf(res.cost.stoch)
stock.perf(res.cost.stoch.fix)
stock.perf(res.cost.fix)

# alpha optimal 
res.alpha.stoch <- tS.inv.L.pois( t = 8, S = S.alpha, l.ini = 200 , dem = dem, L.lambda = lambda)
res.alpha.stoch.fix <- tS.inv.L.pois( t = 8, S = S.alpha.fix, l.ini = 200 , dem = dem, L.lambda = lambda)
res.alpha.fix <- tS.inv.L( t = 8, S = S.alpha.fix, l.ini = 200 , dem = dem, L = lambda)

stock.perf(res.alpha.stoch)
stock.perf(res.alpha.stoch.fix)
stock.perf(res.alpha.fix)


# beta optimal 
res.beta.stoch <- tS.inv.L.pois( t = 8, S = S.beta, l.ini = 200 , dem = dem, L.lambda = lambda)
res.beta.stoch.fix <- tS.inv.L.pois( t = 8, S = S.beta.fix, l.ini = 200 , dem = dem, L.lambda = lambda)
res.beta.fix <- tS.inv.L( t = 8, S = S.beta.fix, l.ini = 200 , dem = dem, L = lambda)

stock.perf(res.beta.stoch)
stock.perf(res.beta.stoch.fix)
stock.perf(res.beta.fix)
