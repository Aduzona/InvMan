#####################################
# ELSP
##################################

# Bombergers stamp problem ##############################
c.sh <- c(0.2708,	7.396,	5.313, 4.167,	116,	11.15,	62.50,	245.8,	37.50,	1.667)*10^-5
c.or <- c(15,	20,	30,	10	,110,	50,	310,	130,	200,	5)
y.vec <- c(40,	40,	80,	160,	8,	8,	2.4,	34,	34,	40)*10 # mult by 10 for real values
p.vec <- c(3000,	800,	950,	750,	200,	600,	240,	130,	200,	1500)*10 # mult by 10 for real values
rho.vec <- y.vec/p.vec
s.vec <- c(0.0125,	0.0125	,0.025,	0.0125	,0.05,	0.025,	0.1,	0.05,	0.075,	0.0125)*10 # mult by 10 for real values

kappa <- 1-sum(y.vec/p.vec)

# objective function
elsp.obj <- function(csh, cor, y, p,t) cor/t + csh*y*(1-y/p)*t/2

# independent solution ##############################
T.vec <- sqrt(2*c.or/c.sh/y.vec/(1-y.vec/p.vec))
bat.vec <- T.vec*y.vec/p.vec + s.vec

T.vec[c(4,8,9)]
bat.vec[c(4,8,9)]

# cost of independent solution
c.cost.ind <- elsp.obj(csh = c.sh, cor = c.or, y = y.vec, p = p.vec , t = T.vec)
sum(c.cost.ind)

# common cycle solution ##############################

T.com.opt <- sqrt(2*sum(c.or)/sum(y.vec*c.sh*(1-y.vec/p.vec)))
T.com.min <- sum(s.vec)/(1-sum(y.vec/p.vec))

bat.vec.com <- T.com.opt*y.vec/p.vec + s.vec
# cost of independent solution
c.cost.com <- elsp.obj(csh = c.sh, cor = c.or, y = y.vec, p = p.vec , t = T.com.opt)
sum(c.cost.com)

# power of two heuristic  ##############################

## 1st iteration
B <- min(T.vec)
m.vec <- sqrt(2*c.or/B^2/c.sh/y.vec/(1-rho.vec))
# round to closest power of 2
m.vec <- 2^round(log(m.vec)/log(2))
# update B
B <- sqrt(2*sum(c.or/m.vec)/sum(c.sh*y.vec*m.vec*(1-rho.vec)))
c.cost.po2.1 <- elsp.obj(csh = c.sh, cor = c.or, y = y.vec, p = p.vec , t = m.vec*B)
sum(c.cost.po2.1)

## 2nd iteration
###change m vector by hand
m.vec <- c(8,2,2,1,2,4,8,1,2,2)
B <- sqrt(2*sum(c.or/m.vec)/sum(c.sh*y.vec*m.vec*(1-rho.vec)))
c.cost.po2.2 <- elsp.obj(csh = c.sh, cor = c.or, y = y.vec, p = p.vec , t = m.vec*B)
sum(c.cost.po2.2)

##############################################################
# Knapsack problem
##############################################################

library(ompr)
library(dplyr)
library(ompr.roi)
library(ROI.plugin.glpk)

n.vec[9] <- 4

M <- sum(n.vec)
K <- max(n.vec)

b.vec2 <- rep(b.vec , n.vec)

model <- MIPModel() %>%
  # create decision variable
  add_variable(a[m, k], m = 1:M, k = 1:K, type = "integer", lb = 0, ub = 1) %>%
  # maximize assignment
  set_objective(sum_expr( a[m, k], m = 1:M, k = 1:K), "max") %>%
  add_constraint(sum_expr(a[m, k], k = 1:K) == 1, m = 1:M) %>%
  add_constraint(sum_expr(a[m, k]*b.vec2[m], m = 1:M) <= B, k = 1:K)

# solve model with GLPK
result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))

result$solution

#####################################################
# Joint replenishment problme
#####################################################
n <- 10
c.or0 <- 100
set.seed(1234)
c.or <- round(runif(n, 5,100))
c.sh <- round(rnorm(n, 2, .5) , 2 )

jrp.obj.fun <- function(m , B, cor, csh, cor0) cor0/B + sum(cor/B/m) + sum(csh*B*m)

# calculate cycle times
T.vec <- sqrt(c.or/c.sh)
# order products
reo.id <- order(T.vec)
T.vec <- T.vec[reo.id]
c.or <- c.or[reo.id]
c.sh <- c.sh[reo.id]

# calculate T^2 and cumu. cost shares
res.mat <- t(cbind(c.or/c.sh,(c.or0 + cumsum(c.or))/cumsum(c.sh)))
# identify break
id.comb <- min(which(res.mat[1,] > res.mat[2,])) - 1

# calculate B
B <- min(T.vec)

# solution with m - integers #######################
m.vec.int <- round(T.vec/B,0)
m.vec.int[1:id.comb] <- 1

# re-optimize B 
B.int <- sqrt(sum(c.or/m.vec.int)/sum(m.vec.int*c.sh))
# total cost 
jrp.obj.fun( m = m.vec.int, B=B.int, cor = c.or,csh = c.sh, cor0 = c.or0)


# solution with m - power of 2 #######################
# Alternative 1: round automatically to closest power of  2
m.vec.pow2 <- T.vec/B
pows2 <- 2^(0:5)
m.vec.pow2 <- sapply(m.vec.pow2, function(x){
  pows2[which.min((pows2 - x)^2)]
})
m.vec.pow2[1:id.comb] <- 1

# re-optimize B 
B.pow2 <- sqrt(sum(c.or/m.vec.pow2)/sum(m.vec.pow2*c.sh))
# total cost 
jrp.obj.fun( m = m.vec.pow2, B=B.pow2, cor = c.or,csh = c.sh,cor0 = c.or0)

# Alternative 2: round manually
m.vec.pow2 <- rep(2, 10)
m.vec.pow2[1:id.comb] <- 1
m.vec.pow2[10] <- 4

# re-optimize B 
B.pow2 <- sqrt(sum(c.or/m.vec.pow2)/sum(m.vec.pow2*c.sh))
# total cost 
jrp.obj.fun( m = m.vec.pow2, B=B.pow2, cor = c.or,csh = c.sh,cor0 = c.or0)

