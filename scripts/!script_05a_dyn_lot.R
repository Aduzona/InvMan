len <- 30

y.err <- rnorm(len, 0 , .001)
y.const <- rep(50, len)
x.seas <- 0:(len-1)%%6 + 1 
y.seas <- (x.seas-3.5)*2 # quarterly fluctuations
# ts values
y.com <- y.const + y.seas #+ y.err 
y.com.ts <- ts(y.com, deltat = 1/6)

win.metafile("sto_dyn_ls_exa_dem.wmf", width=5, height = 3)
par(family="serif", mar = c(4.25,4.25,.1,.1), bg ="white")
plot(y.com, xlab="Time", type="l", lwd=2, ylab="demand", cex.lab = 1.5, cex.axis = 1.25)
abline(h=50, col="grey", lwd=2)
dev.off()

res <- tS.lager.wbz(dem = y.com.ts, t=3, wbz=1, l.ini = 150, S = 200)

par(family="serif", mar = c(4.25,4.25,3,.1), bg ="white")
plot(tail(res$i,24), type="s", lwd=2, ylab="on-hand stock", cex.lab = 1.75, cex.axis = 1.5, xlab="Time", main=parse(text="list(S==200,t==3,L==1)"), cex.main = 1.75)
abline(h=0, lwd=2, col="grey")


#######################################################
# stochastic MPLSP - alpha service level
#######################################################
c.or <- 50
c.sh <- 1
mu.vec <- c(20,12,24,17,12)
sig.vec <- c(6,3,7,5,2)
alpha <- 0.95

data.func.alpha <- function(c.or, c.sh, mus, sds, alpha){
  
  n <- length(mus)
  seqs <- n:1
  
  # mean mat
  mean.mat <- matrix(NA, ncol=n, nrow=n, byrow=T) 
  mean.mat[lower.tri(mean.mat, diag = T)] <- unlist(sapply(seqs, function(x) cumsum(tail(mus, x) ) ))
  mean.mat <- t(mean.mat)
  # sds
  sds.mat <- matrix(NA, ncol=n, nrow=n, byrow=T) 
  sds.mat[lower.tri(mean.mat, diag = T)] <- unlist(sapply(seqs, function(x) sqrt(cumsum(tail(sds^2, x) ) ) ))
  sds.mat <- t(sds.mat)
  
  # safety factor
  z <- qnorm(alpha)
  
  # S mat
  S.mat <- mean.mat + z*sds.mat
  
  # stock holding multiplier
  mult.mat <- matrix(NA, ncol=n, nrow=n) 
  mult.mat[lower.tri(mult.mat, diag = T)] <- unlist( sapply(1:n, function(x) tail(1:n,n-x+1)-x+1 ) )
  mult.mat <- t(mult.mat)
  # cumulated means
  mmeans.mat <- matrix(NA, ncol=n, nrow=n) 
  mmeans.mat[lower.tri(mmeans.mat, diag = T)] <- unlist(sapply(seqs, function(x) sapply(1:x , function(y) sum(y:1 * head(tail(mu.vec, x),y) ) )))
  mmeans.mat <- t(mmeans.mat)
  
  # calculate costs
  cost.mat <- c.or + c.sh * (mult.mat * S.mat -  mmeans.mat)
  
  # return results
  return(list(mean.mat = mean.mat, sd.mat = sds.mat, S.mat = S.mat, cost.mat = cost.mat))
  
}

data.alpha <- data.func.alpha(c.or = c.or, c.sh = c.sh,mus = mu.vec, sds = sig.vec, alpha = alpha )

# wagner within algorithm ################################
ww.opt.func <- function(c.mat){
  n <- nrow(c.mat)
  res.mat <- matrix(NA, ncol=n, nrow=n)
  res.id <- numeric(n)
  # initializing
  res.mat[1,] <- c.mat[1,]
  for(i in 2:n){
    res.id[i-1] <- which.min(res.mat[,i-1])
    res.mat[i,] <- c.mat[i,] + res.mat[res.id[i-1],i-1]
  }
  res.id[n] <- which.min(res.mat[,n])
  polic <- res.id[n]
  repeat{
    polic <- c(polic , res.id[tail(polic,1)-1] )
    if(tail(polic,1) == 1) break
  }
  polic <- rev(polic)
  
  return(list(cum.cost = res.mat, ord.per = polic, ids = res.id))
}

# result alpha service level
res.alpha <- ww.opt.func(data.alpha$cost.mat)

#######################################################
# stochastic MPLSP - beta service level
#######################################################

beta <- 0.975

data.func.beta <- function(c.or, c.sh, mus, sds, beta){

    # beta-Servicelevel
  v.fun <- function(x, mu = 0, sigma = 1) {
    sigma * (dnorm((x-mu)/sigma) - (x-mu)/sigma * (1-pnorm((x-mu)/sigma) ))
  }
  # optimizer beta SL
  v.fun.opt.loss <- function(x, loss = loss.norm, ...) (loss - v.fun(x, ...))^2 
  
  n <- length(mus)
  seqs <- n:1
  
  # mean mat
  mean.mat <- matrix(NA, ncol=n, nrow=n, byrow=T) 
  mean.mat[lower.tri(mean.mat, diag = T)] <- unlist(sapply(seqs, function(x) cumsum(tail(mus, x) ) ))
  mean.mat <- t(mean.mat)
  # sds
  sds.mat <- matrix(NA, ncol=n, nrow=n, byrow=T) 
  sds.mat[lower.tri(mean.mat, diag = T)] <- unlist(sapply(seqs, function(x) sqrt(cumsum(tail(sds^2, x) ) ) ))
  sds.mat <- t(sds.mat)
  
  # S mat
  S.mat <- matrix(NA, ncol=n, nrow=n) 
  S.mat[upper.tri(S.mat, diag = T)] <- 
    sapply(which(upper.tri(mean.mat, diag = T)), function(x) {
      optim(v.fun.opt.loss, loss = mean.mat[x]*(1-beta), mu = mean.mat[x], sigma = sds.mat[x], par = mean.mat[x], lower = 0, upper = mean.mat[x]+5*sds.mat[x], method="L-BFGS-B")$par
    })
  
  # stock holding multiplier
  mult.mat <- matrix(NA, ncol=n, nrow=n) 
  mult.mat[lower.tri(mult.mat, diag = T)] <- unlist( sapply(1:n, function(x) tail(1:n,n-x+1)-x+1 ) )
  mult.mat <- t(mult.mat)
  # cumulated means
  mmeans.mat <- matrix(NA, ncol=n, nrow=n) 
  mmeans.mat[lower.tri(mmeans.mat, diag = T)] <- unlist(sapply(seqs, function(x) sapply(1:x , function(y) sum(y:1 * head(tail(mu.vec, x),y) ) )))
  mmeans.mat <- t(mmeans.mat)
  
  # calculate costs
  cost.mat <- c.or + c.sh * (mult.mat * S.mat -  mmeans.mat)
  
  # return results
  return(list(mean.mat = mean.mat, sd.mat = sds.mat, S.mat = S.mat, cost.mat = cost.mat))
  
}

data.beta <- data.func.beta(c.or = c.or, c.sh = c.sh,mus = mu.vec, sds = sig.vec, beta = beta )

# result beta service level
res.beta <- ww.opt.func(data.beta$cost.mat)

###############################################################
# Cost - optimization
###############################################################

c.bo <- 5

# data preparation cost
data.func.cost <- function(c.or, c.sh, mus, sds, c.bo){
  
  
  # beta-Servicelevel
  v.fun <- function(x, mu = 0, sigma = 1) {
    sigma * (dnorm((x-mu)/sigma) - (x-mu)/sigma * (1-pnorm((x-mu)/sigma) ))
  }
  # optimizer beta SL
  v.fun.opt.loss <- function(x, loss = loss.norm, ...) (loss - v.fun(x, ...))^2 
  
  # newsvendor obj. func.
  zf.news <- function(x, cu, co, mu, sigma){
    tmp1 <- function(y) co*(x-y)*dnorm(y, mean=mu, sd = sigma)
    tmp2 <- function(y) cu*(y-x)*dnorm(y, mean=mu, sd = sigma)
    integrate(tmp1, lower = 0, upper=x)$value + integrate(tmp2, lower = x, upper=Inf)$value
  } 
  
  si.opt.fun <- function(x, csh, cbo, mus, sigs){
    n <- length(mus)
    (sum(sapply(1:n , function(y) pnorm(x, mean = mus[y], sd = sigs[y]))) - n*cbo/(cbo+csh))^2
  }
  
  cost.fun <- function(s, csh, cbo, mus, sigs){
    sum(sapply(1:length(mus), function(x) zf.news(x= s, mu = mus[x], sigma  = sigs[x], co = csh, cu = cbo )))
  }
  
  
  
  n <- length(mus)
  seqs <- n:1
  
  # mean mat
  mean.mat <- matrix(NA, ncol=n, nrow=n, byrow=T) 
  mean.mat[lower.tri(mean.mat, diag = T)] <- unlist(sapply(seqs, function(x) cumsum(tail(mus, x) ) ))
  mean.mat <- t(mean.mat)
  # sds
  sds.mat <- matrix(NA, ncol=n, nrow=n, byrow=T) 
  sds.mat[lower.tri(mean.mat, diag = T)] <- unlist(sapply(seqs, function(x) sqrt(cumsum(tail(sds^2, x) ) ) ))
  sds.mat <- t(sds.mat)
  
  # S mat
  S.mat <-  matrix(NA, ncol=ncol(mean.mat), nrow=nrow(mean.mat)) 
  cost.mat <-  matrix(NA, ncol=ncol(mean.mat), nrow=nrow(mean.mat)) 
  for(i in 1:nrow(mean.mat) ){
    tmp.mu <- mean.mat[i,]
    tmp.sig <- sds.mat[i,]
    id.j <- which(!is.na(tmp.mu))
    tmp.mu <- tmp.mu[!is.na(tmp.mu)]
    tmp.sig <- tmp.sig[!is.na(tmp.sig)]
    for(j in 1:length(tmp.sig)){
      S.mat[ i, id.j[j] ] <- optim(si.opt.fun, csh = c.sh, cbo = c.bo, mus = tmp.mu[1:j], sigs = tmp.sig[1:j], par = tmp.mu[j], lower = 0, upper = 5*sum(tmp.mu[1:j]), method="L-BFGS-B")$par
      cost.mat[ i, id.j[j] ] <- c.or + cost.fun(s = S.mat[ i, id.j[j] ], csh = c.sh, cbo = c.bo, mus = tmp.mu[1:j], sigs = tmp.sig[1:j])
    }
  } 
  
  
  
  # return results
  return(list(mean.mat = mean.mat, sd.mat = sds.mat, S.mat = S.mat, cost.mat = cost.mat))
  
}

# data cost optimization
data.cost <- data.func.cost(c.or = c.or, c.sh = c.sh,mus = mu.vec, sds = sig.vec, c.bo =  c.bo )

# result cost optimization
res.cost <- ww.opt.func(data.cost$cost.mat)


###########################################################
# simulation
###########################################################

# alpha SL
c.bo <- 0
S1 <- q.mat[1,2]
S3 <- q.mat[3,5]
# beta SL
c.bo <- 0
S1 <- S.mat[1,2]
S3 <- S.mat[3,5]
# BO cost
S1 <- S.mat[1,2]
S3 <- S.mat[3,5]


reps <- 10000
res <- NULL
  
for(i in 1:reps){
  tmp.dem <- sapply(1:5, function(x) rnorm(1, mean = mu.vec[x], sd = sig.vec[x]))
  i.vec <- numeric(5)
  i.vec[1:2] <- S1 - cumsum(tmp.dem[1:2])
  i.vec[3:5] <- S3 - cumsum(tmp.dem[3:5])
  tmp.cost <- sum(i.vec[i.vec > 0]*c.sh) + 2*c.or - sum(i.vec[i.vec < 0]*c.bo)
  tmp.sl.alpha <- sum(i.vec[c(2,5)] >= 0)/2
  tmp.sl.beta <- ((sum(tmp.dem[1:2]) + ifelse(i.vec[2] < 0, i.vec[2],0 ) )/sum(tmp.dem[1:2]) + (sum(tmp.dem[3:5]) + ifelse(i.vec[5] < 0,i.vec[5],0 ) )/sum(tmp.dem[3:5])) /2
  res <- rbind(res, c(tmp.cost, tmp.sl.alpha, tmp.sl.beta, i.vec, tmp.dem))
  }

colMeans(res)

# beta service level
(1+sum(res[res[,5] < 0, 5])/ reps/sum(mu.vec[1:2]) +
1+sum(res[res[,8] < 0, 8])/ reps/sum(mu.vec[3:5]) )/2


