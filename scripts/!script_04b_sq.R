# sq inventory simulation
sq.inv.L <- function(d, s, q, l.ini, n, L){
  
  tmp.mat <- matrix(0, nrow= n+1+L, ncol= 5 )
  colnames(tmp.mat) <- c("demand","inv.level","inv.pos", "order", "delivery")
  
  tmp.mat[,"demand"] <- c(0,d,rep(0,L))
  tmp.mat[1,"Lagerbestand"] <- l.ini
  tmp.mat[1,"inv.pos"] <- l.ini
  
  for(i in 2:(n+1)){
    
    if(tmp.mat[i-1, "inv.pos"] <= s){
      tmp.mat[i, "order"] <- q
      tmp.mat[i+L, "delivery"] <- q
    }
    
    tmp.mat[i, "inv.pos"] <- tmp.mat[i-1, "inv.pos"] - tmp.mat[i, "demand"] + tmp.mat[i, "order"]
    tmp.mat[i, "inv.level"] <- tmp.mat[i-1, "inv.level"] - tmp.mat[i, "demand"] + tmp.mat[i, "delivery"]
  }
  return(tmp.mat[2:(n+1),]) 
}

# sq inventory optimization
sq.inv.opt <- function(pars, demand, wwbz, llini, obj = "cost", co = 0, cu = 0, B = 0, alpha.z = 95, beta.z = 99){
  
  sim.R <- sq.lager(d = demand, s = pars[1], q = pars[2], l.ini = llini, n = length(demand), L = wwbz)
  
  tmp <- sim.R[-c(1),]
  tmp.l <- tmp[,"inv.level"]
  tmp.id <- 1:length(tmp.l)
  tmp.id.sub <- tmp[,"delivery"] > 0
  tmp.id.sub.w <- which(tmp.id.sub) - 1
  tmp.id.sub.w <-tmp.id.sub.w[tmp.id.sub.w >0]
  
  # Alpha-SL
  sim.alpha <- (100*(1 - sum( tmp.l[tmp.id.sub.w] < 0)/length(tmp.id.sub.w)))
  sim.alpha.diff <- (alpha.z - sim.alpha)^2
  # Beta-SL
  sim.beta <- (100*(1 - sum(-tmp.l[ tmp.l < 0])/sum(tmp[,"demand"])))
  sim.beta.diff <- (beta.z - sim.beta)^2
  
  # cost
  sim.cost <- ( sum(tmp[,"order"] > 0)*B +  sum(tmp.l[ tmp.l > 0]) * co + sum(-tmp.l[ tmp.l < 0]) * cu)/length(tmp.l)
  
  return(ifelse(obj == "cost", sim.cost, ifelse(obj == "alpha", sim.alpha.diff, sim.beta.diff)) )
  
}

# Normal loss functions - direct
v.fun.direct.norm <- function(s, mu = 0, sigma = 1){
  integrate(function(y) (y-s)*dnorm(y, mean=mu, sd = sigma), lower = s, upper=Inf)$value
}

# deviation Beta vs. exp. loss 
v.diff.norm <- function(x, beta, mu = 0, sigma = 1, L,  c.or, c.sh){
  (v.fun.direct.norm(x[1], mu = L*mu, sigma = sqrt(L)*sigma) - (1-beta)*x[2])
}



################################################
# Optimization under beta service level constraint

sq.cost <- function(x, c.or, c.sh,mu, L, beta, sigma){
  c.sh *(x[2]/2 + x[1] - mu*L) + c.or*mu/x[2]
}

sq.cost(x=c(772, 1000), mu = 100, L = 8, c.or = 120, c.sh = 0.024)


# solve with Rsolnp ###################

library(Rsolnp)

# general-purpose constrined NLP optimizer
x.opt.solnp <- solnp(pars = c(80,300), fun = sq.cost, eqfun = v.diff.norm, eqB = 0, LB = c(10,10), UB = c(2000,3000), mu = 20, sigma = 4, L = 4, c.or = 100, c.sh = 0.1, beta = .99 )$pars

# setting of initial values for search important
x.opt.solnp <- solnp(pars = c(700,1000), fun = sq.cost, eqfun = v.diff.norm, eqB = 0, LB = c(100,100), UB = c(2000,3000), mu = 100, sigma = 30, L = 8, c.or = 120, c.sh = 0.024, beta = .95 )$pars

# solve with optim ###################

sq.cost.lagrange.norm <- function(x, mu, c.or, c.sh, L, sigma, beta){
  # x...(s,q,lambda)
  c.sh *(x[2]/2 + x[1] - mu*L) + c.or*mu/x[2] + x[3]*( ( v.fun.direct.norm(x[1], mu = L*mu, sigma = sqrt(L)*sigma) -(1-beta)*x[2] )^2 )
}

sq.cost.lagrange.norm(x=c(765,1085,0), mu = 100, L = 8, c.or = 120, c.sh = 0.024, sigma = 30, beta = 0.95)

# test result
v.fun.direct.norm(765, mu = 8*100, sigma = sqrt(8)*30)
(1-.95)*1085

x.opt.optim <- head(optim(fn = sq.cost.lagrange.norm, par = c(700, 1000, 10), lower = c(500, 800, 10), upper = c(1000, 2000, 11), method="L-BFGS-B", mu = 100, sigma = 30, L = 8, c.or = 120, c.sh = 0.024, beta = .95)$par,2)

x.opt.optim <- head(optim(fn = sq.cost.lagrange.norm, par = c(80, 800, 1), lower = c(50, 100, 1), upper = c(1000, 2000, 2), method="L-BFGS-B", mu = 20, sigma = 4, L = 4, c.or = 100, c.sh = 0.1, beta = .99)$par,2)

# with iterative algorithm ###################

v.diff.norm.q <- function(s, q, beta, mu = 0, sigma = 1, L,  c.or, c.sh){
  (v.fun.direct.norm(s, mu = L*mu, sigma = sqrt(L)*sigma) - (1-beta)*q)^2
}


iter.sq <- function(mu = 100, sigma = 30, L = 8, c.or = 120, c.sh = 0.024, beta = .95){
  
  iter <- 1
  q.opt.old <- Inf
  lambda.opt <- 0
  ef <- 0

  repeat{
    q.opt <- sqrt(2*(mu*c.or+lambda.opt*ef) /(c.sh))
    if( abs(q.opt.old - q.opt) <= 1e-4 )  return(c(s.opt, q.opt))
    ef <- (1-beta)*q.opt
    s.opt <- optim(fn = v.diff.norm.q, par = mu*L, lower = 0, upper = mu*L*5, method="L-BFGS-B", mu = mu, sigma = sigma, q = q.opt, beta = beta, L = L)$par
    lambda.opt <- c.sh*q.opt/(1-pnorm(s.opt, mean = L*mu, sd = sqrt(L)*sigma ))
    res <- c(iter, q.opt, ef, s.opt, lambda.opt)
    names(res) <- c("iter","q","ef","s","lambda")
    print(res)
    q.opt.old <- q.opt
    iter <- iter + 1
    if(iter > 100) break
  }
  
}

x.opt.iter <- iter.sq()

x.opt.iter <- iter.sq(mu = 20, sigma = 4, L = 4, c.or = 100, c.sh = 0.1, beta = .99)
sq.cost(x = x.opt.iter, mu = 20, sigma = 4, L = 4, c.or = 100, c.sh = 0.1, beta = .99)

# comparisonof solutions
sq.cost(x = x.opt.iter, mu = 100, L = 8, c.or = 120, c.sh = 0.024)
sq.cost(x = x.opt.solnp, mu = 100, L = 8, c.or = 120, c.sh = 0.024)
sq.cost(x = x.opt.optim, mu = 100, L = 8, c.or = 120, c.sh = 0.024)

v.diff.norm.q(x.opt.iter[1], x.opt.iter[2], mu = 100, L = 8, beta=.95, sigma = 30)

# 3D plot cost function #############################

library(plotly)

sq.cost.plot <- function(s, q, c.or, c.sh,mu, L){
  c.sh *(q/2 + s - mu*L) + c.or*mu/q
}
beta.plot <- function(s,q, mu, L, sigma) 1 - vv.fun.direct.norm(s, mu = L*mu, sigma = sqrt(L)*sigma)/q

# vectorize loss diff function
vv.fun.direct.norm <- Vectorize(v.fun.direct.norm)

# create plotting data
s.seq <- seq(82, 83, length= 120)
q.seq <- seq(200, 210, length= 120)
sq.mat <- expand.grid(s.seq, q.seq)
z.mat <- apply(sq.mat, 1, function(x) c(
  sq.cost(x=x, mu = 20,  L = 4, c.or = 100, c.sh = 0.1),
  1-v.fun.direct.norm(x[1], mu = 4*20, sigma = sqrt(4)*4)/x[2]
) )

sq.df <- data.frame(s=sq.mat[,1], q = sq.mat[,2], cost = z.mat[1,], beta  = z.mat[2,], beta.diff = (z.mat[2,] - .99)^2)

# expand to matrices
z.mat <- outer(s.seq, q.seq, FUN=sq.cost.plot, mu = 20,  L = 4, c.or = 100, c.sh = 0.1)
beta.mat <- outer(s.seq, q.seq, FUN = beta.plot, mu = 20,  L = 4, sigma = 4)


library(rgl)
sub.df <- subset(sq.df , beta.diff  <= 1e-11)

plot3d(x=rep(x.opt.iter[1],2), y=rep(x.opt.iter[2],2), z=rep(sq.cost(x = x.opt.iter, mu = 20, sigma = 4, L = 4, c.or = 100, c.sh = 0.1, beta = .99),2), xlim = range(s.seq), ylim=range(q.seq), zlim=range(sq.df$cost), xlab="s", ylab="q", zlab ="cost", pch=14, size=8)
lines3d(sub.df$s, sub.df$q, sub.df$cost)

rbPal <- colorRampPalette(c('green','red'))
rb.col <- rbPal(20)[as.numeric(cut(z.mat,breaks = 20))]
surface3d(s.seq, q.seq, z.mat, col=rb.col, alpha=.5, back = "lines")
rgl.bringtotop()

# contour plot ###################

library(plotly)
#  colorscale = 'Jet',

fig <- plot_ly(sq.df, x=~q, y=~s, z =~cost, type = "contour", 
        contours = list(showlabels = TRUE),showlegend = FALSE,showscale =F,
        line = list(smoothing = 0.05),
        colorscale = 'Jet')
fig <- fig %>% add_trace( type = 'contour',
                          x = sq.df$q,
                          y = sq.df$s,
                          z = sq.df$beta,
                          showlegend = FALSE,
                          opacity = .99,
                          autocolorscale = F,
                          autocontour =F,
                          showscale =F,
                          line = list(width=2),
                          contours=list(
                            start = .988,
                            end = .998,
                            size = .001,
                            showlabels = T,
                            coloring = 'none'
                          ))
# add opt. point don't work after update
# fig <- fig %>% add_trace( type = 'scatter',showlegend = FALSE, size=3,showscale =F, x = x.opt.iter[2], y = x.opt.iter[1])

fig


##############################
# density plot gamma & normal distribution 
##############################

mu <- 10
sigma <- 3

x.seq <- seq(0, mu+3.5*sigma, length.out = 100)
y.norm.seq <- dnorm(x.seq, mean=mu, sd = sigma)
y.gamma.seq <- dgamma(x.seq, shape = mu^2/sigma^2, rate = mu/sigma^2)

par(family="serif", mar = c(4.5,4.5,3,.1), bg ="white")
plot(x.seq, y.norm.seq, type="l", col="black", lwd=2, ylim=range(y.norm.seq,y.gamma.seq), main = parse(text = "list(mu==10,sigma==3)"), ylab = "density", xlab="y", cex.axis = 1.25, cex.lab = 1.5, cex.main = 1.5)
lines(x.seq, y.gamma.seq,  col="red", lwd=2)
legend("topleft", lwd=2, col=c("black","red"), legend =c(parse(text = "y%~%N(10,3)"),parse(text = "y%~%Gamma(alpha==11.1,beta==1.1)")), bty="n", bg=NULL, cex=1.25)



mu <- 10
sigma <- 5

x.seq <- seq(0, mu+3.5*sigma, length.out = 100)
y.norm.seq <- dnorm(x.seq, mean=mu, sd = sigma)
y.gamma.seq <- dgamma(x.seq, shape = mu^2/sigma^2, rate = mu/sigma^2)

par(family="serif", mar = c(4.5,4.5,3,.1), bg ="white")
plot(x.seq, y.norm.seq, type="l", col="black", lwd=2, ylim=range(y.norm.seq,y.gamma.seq), main = parse(text = "list(mu==10,sigma==5)"), ylab = "density", xlab="y", cex.axis = 1.25, cex.lab = 1.5, cex.main = 1.5)
lines(x.seq, y.gamma.seq,  col="red", lwd=2)
legend("topleft", lwd=2, col=c("black","red"), legend =c(parse(text = "y%~%N(10,5)"),parse(text = "y%~%Gamma(alpha==4,beta==0.4)")), bty="n", bg=NULL, cex=1.25)



mu <- 10
sigma <- 7

x.seq <- seq(0, mu+3.5*sigma, length.out = 100)
y.norm.seq <- dnorm(x.seq, mean=mu, sd = sigma)
y.gamma.seq <- dgamma(x.seq, shape = mu^2/sigma^2, rate = mu/sigma^2)

par(family="serif", mar = c(4.5,4.5,3,.1), bg ="white")
plot(x.seq, y.norm.seq, type="l", col="black", lwd=2, ylim=range(y.norm.seq,y.gamma.seq), main = parse(text = "list(mu==10,sigma==7)"), ylab = "density", xlab="y", cex.axis = 1.25, cex.lab = 1.5, cex.main = 1.5)
lines(x.seq, y.gamma.seq,  col="red", lwd=2)
legend("topleft", lwd=2, col=c("black","red"), legend =c(parse(text = "y%~%N(10,7)"),parse(text = "y%~%Gamma(alpha==2.04,beta==0.204)")), bty="n", bg=NULL, cex=1.25)


###############################################
# Gamma distributions
###############################################
alpha.vec <- c(.5,1,2)
beta.vec <- c(.5,1,2)

sets.mat <- expand.grid(alpha.vec, beta.vec)
sets.mat <- sets.mat[-c(4,5,7,8),]
sets.mat <- rbind(sets.mat, c(4,1))

x.vals <- seq(0,10,length.out = 100)

leg.vals <- apply(sets.mat, 1, function(x){
  paste("list(alpha==", x[1],",beta==",x[2],")", sep="" )
})

par(family="serif", mar = c(4.5,4.5,.1,.1), bg ="white")
plot(x.vals, dgamma(x.vals, shape = sets.mat[1,1], rate = sets.mat[1,2]), type="l",  xlab="y", ylab="density", lwd=3 , cex.main = 1.75, cex.axis = 1.25, cex.lab = 1.75)
for(i in 2:nrow(sets.mat)){
  lines(x.vals, dgamma(x.vals, shape = sets.mat[i,1], rate = sets.mat[i,2]) , col=i, type="l", lwd = 3)
}
legend("topright", lwd=3, col=1:length(leg.vals), legend = parse(text=leg.vals) , cex=1.5 , bty="n", bg=NULL)

###################################
# stochastic nature of demand
###################################

# inter-arrival times
irr.times <- rexp(10, 2.5)
# order sizes
ord.sizes <- rpois(10, 2)+1

# arrival times
arr.times <- cumsum(irr.times)
# plot
par(family="serif", mar = c(4,4,3,.1), bg ="white")
plot(arr.times, ord.sizes, type="h", xlim=c(0,ceiling(max(arr.times))), cex.lab = 1.75, cex.axis = 1.25 , xlab="time", ylab="order size in units", lwd=1.75, main = parse(text = "list(I%~%Exp(2.5),X%~%Pois(2)+1)"), cex.main = 1.75, ylim=c(0,max(ord.sizes)+1),xaxs ="i", yaxs="i")
points(arr.times, ord.sizes, pch=15)
abline(v = 1:4 , lwd=2, col="darkgrey", lty =2)


# simulation #############
reps <- 100000
z <- numeric(reps)

for(i in 1:reps){
  tmp <- rexp(8, 2)
  if(tmp[1] > 1){
    z[i] <- 0
  }
  else{
    z[i] <- max(which(cumsum(tmp) < 1))
  }
}

par(family="serif", mar = c(4.25,4.25,.1,.1), bg ="white")
plot(0:8, table(z)/reps, type="h", cex.lab = 1.75, cex.axis = 1.25 , xlab="demand in units", ylab="frequency", lwd=2, cex.main = 1.75, ylim=c(0,.3), yaxs="i")
lines(0:8+.125, dpois(0:8,2), pch=15, lwd=2, type="h", col="red")
legend("topright", lwd =2, col=c("black","red"), legend =c("simulated","Pois(2)"), bty="n" , bg=NULL, cex=1.75)

#############################################
# sq with Gamma distribution
#############################################

v.fun.direct.gamma <- function(s, mu = 0, sigma = 1){
  integrate(function(y) (y-s)*dgamma(y, shape = mu^2/sigma^2, rate = mu/sigma^2), lower = s, upper=Inf)$value
}

v.fun.gamma <- function(s, mu = 0, sigma = 1){
  alpha <- mu^2/sigma^2
  beta <- mu/sigma^2
  alpha/beta * (1 - pgamma(s, shape = alpha+1, rate = beta)) - s*(1 - pgamma(s, shape = alpha, rate = beta))
}


# comparison normal and gamma loss functions ################
v.fun.direct.norm(s = 20, mu = 7, sigma = 7)
v.fun.direct.gamma(s= 20, mu = 7, sigma = 7)

par(mfrow = c(3,1))
x.seq <- 10:60
for(sig in c(5, 10, 20)){
  par(family="serif", mar = c(4.25,4.25,3,.1), bg ="white")
  plot(x.seq, v.fun.gamma(s= x.seq, mu = 20, sigma = sig), type="l", lwd=2, col="black", ylim=c(0,15), xlab="reorder point s", ylab="expected loss L(s)", cex.lab = 1.75, cex.axis = 1.25 , main = parse(text = paste("list(mu==20,sigma==",sig,")", sep="")), cex.main =1.75)
  abline(v=20, lty=2, lwd=2, col="grey")
  lines(x.seq, vv.fun.direct.norm(x.seq,  mu=20, sigma =sig), lwd=2, col="red")
  legend("topright", lwd =2, col=c("black","red"), legend =c(parse(text = paste("Gamma(list(alpha==",round(20^2/sig^2,2),",beta==",round(20/sig^2,2),"))")), parse(text=paste("N(20,",sig,"^2)"))), bty="n" , bg=NULL, cex=1.75)
}



iter.sq.mult <- function(mu = 100, sigma = 30, L = 8, c.or = 120, c.sh = 0.024, beta = .95, dist = "norm"){
  
  # beta vs. loss normal
  v.diff.norm.q <- function(s, q, beta, mu = 0, sigma = 1, L,  c.or, c.sh){
    (v.fun.direct.norm(s, mu = L*mu, sigma = sqrt(L)*sigma) - (1-beta)*q)^2
  }
  # beta vs. loss gamma
  v.diff.gamma.q <- function(s, q, beta, mu = 0, sigma = 1, L,  c.or, c.sh){
    (v.fun.gamma(s, mu = L*mu, sigma = sqrt(L)*sigma) - (1-beta)*q)^2
  }
  
  if(!(dist %in% c("gamma","norm"))) stop("Neither normal nor gamma distribution defined")
     
  iter <- 1
  q.opt.old <- Inf
  lambda.opt <- 0
  ef <- 0
  
  repeat{
    q.opt <- sqrt(2*(mu*c.or+lambda.opt*ef) /(c.sh))
    if( abs(q.opt.old - q.opt) <= 1e-4 )  return(c(s.opt, q.opt))
    ef <- (1-beta)*q.opt
    if(dist == "norm"){
      s.opt <- optim(fn = v.diff.norm.q, par = mu*L, lower = 0, upper = mu*L*5, method="L-BFGS-B", mu = mu, sigma = sigma, q = q.opt, beta = beta, L = L)$par
      lambda.opt <- c.sh*q.opt/(1-pnorm(s.opt, mean = L*mu, sd = sqrt(L)*sigma ))
    }
    if(dist == "gamma"){
      s.opt <- optim(fn = v.diff.gamma.q, par = mu*L, lower = 0, upper = mu*L*5, method="L-BFGS-B", mu = mu, sigma = sigma, q = q.opt, beta = beta, L = L)$par
      lambda.opt <- c.sh*q.opt/(1-pgamma(s.opt, shape = (L*mu)^2/L/sigma^2, rate = L*mu/L/sigma^2 ))
    }
      
    res <- c(iter, q.opt, ef, s.opt, lambda.opt)
    names(res) <- c("iter","q","ef","s","lambda")
    #print(res)
    q.opt.old <- q.opt
    iter <- iter + 1
    if(iter > 100) break
  }
  
}

mu <- 20
wbz <- 4
beta.vec <- c(0.95, 0.99)
sig.vec <- 4:20
exp.des <- expand.grid(beta.vec, sig.vec)

res.list <- NULL

for(i in 1:nrow(exp.des)){
  sig <- exp.des[i,2]
  bet <- exp.des[i,1]
  x.opt.norm <- iter.sq.mult(mu = mu, sigma = sig, L = wbz, c.or = 100, c.sh = 0.1, beta = bet, dist = "norm")
  
  x.opt.gamma <- iter.sq.mult(mu = mu, sigma = sig, L = wbz, c.or = 100, c.sh = 0.1, beta = bet, dist = "gamma")
  
  cost.norm <- sq.cost(x = x.opt.norm, mu = mu, sigma = sig, L = wbz, c.or = 100, c.sh = 0.1, beta = bet)
  cost.gamma <- sq.cost(x = x.opt.gamma, mu = mu, sigma = sig, L = wbz, c.or = 100, c.sh = 0.1, beta = bet)
  
  res.norm <- data.frame(beta = bet, sigma = sig, dist = "norm" , s= x.opt.norm[1] , q= x.opt.norm[2], cost = cost.norm)
  res.gamma <- data.frame(beta = bet, sigma = sig, dist = "gamma", s=x.opt.gamma[1], q = x.opt.gamma[2], cost = cost.gamma)
  
  res.list <- rbind(res.list, res.norm, res.gamma)
  
}

colnames(res.list) <- c("beta", "sigma","distribution","s","q","cost")
res.list <- data.frame(res.list)

library(plotly)

style <- list(family = "Times New Roman", size =20, color="black")

plot_ly(res.list, x= ~sigma, y= ~s , type="scatter", mode="markers+lines", linetype = ~beta, color = ~distribution) %>%
  layout(xaxis = list(title = "sigma", titlefont = style), yaxis=list(title="reorder point s", titlefont = style), plot_bgcolor ="#CBCBCB60")


plot_ly(res.list, x= ~sigma, y= ~cost , type="scatter", mode="markers+lines", linetype = ~beta, color = ~distribution) %>%
  layout(xaxis = list(title = "sigma", titlefont = style), yaxis=list(title="period cost", titlefont = style), plot_bgcolor ="#CBCBCB60")

