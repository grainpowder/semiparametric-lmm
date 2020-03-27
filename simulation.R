source("./misc/make_Z.R") # Builds design matrix for random effects
f1 = function(x) 2*x+sin(pi*x)
f2 = function(x) 2*x+sin(2*pi*x)
f3 = function(x) 2*x+sin(3*pi*x)

# 1. BSAR -----------------------------------------------------------------
source("./bsar.R")
library(bsamGP)

# 1-1. Nonparametric Regression -------------------------------------------
set.seed(10)
N = 300
x = runif(N)*3; ord = order(x)

# 1-1-1. y = f1(x) -----------------------------------------------------------
y = f1(x) + rnorm(N, sd=0.5)

# VB
vb_result = bsar_vb(y,x)

# MCMC
fout = bsar(y ~ fs(x), nbasis=10, shape="Free")
fit_fout = fitted(fout)

# result
result111 = list(
  beta = data.frame(VB=vb_result$mubeta.q, MCMC=mean(fit_fout$wbeta$mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,1,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - vb_result$mubeta.q
plot(x, res, main="VB: 111", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean
plot(x, res, main="MCMC: 111", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 1-1-2. y = f2(x) -----------------------------------------------------------
y = f2(x) + rnorm(N, sd=0.5)

# VB
vb_result = bsar_vb(y,x)

# MCMC
fout = bsar(y ~ fs(x), nbasis=10, shape="Free")
fit_fout = fitted(fout)

# result
result112 = list(
  beta = data.frame(VB=vb_result$mubeta.q, MCMC=mean(fit_fout$wbeta$mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,1,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - vb_result$mubeta.q
plot(x, res, main="VB: 112", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean
plot(x, res, main="MCMC: 112", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 1-1-3. y = f3(x) -----------------------------------------------------------
y = f3(x) + rnorm(N, sd=0.5)


# VB
vb_result = bsar_vb(y,x)

# MCMC
fout = bsar(y ~ fs(x), nbasis=10, shape="Free")
fit_fout = fitted(fout)

# result
result113 = list(
  beta = data.frame(VB=vb_result$mubeta.q, MCMC=mean(fit_fout$wbeta$mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,1,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - vb_result$mubeta.q
plot(x, res, main="VB: 113", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean
plot(x, res, main="MCMC: 113", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)



# 1-2. Semiparametric regression ------------------------------------------
set.seed(10)
N = 300; D = 5
x = runif(N)*3; ord = order(x)
w = matrix(rnorm(N*D), ncol=D)
beta = rnorm(D+1)

# 1-2-1. y = wb + f1(x) ---------------------------------------------------
y = cbind(1,w)%*%beta + f1(x) + rnorm(N, sd=0.5)

# VB
vb_result = bsar_vb(y,x,w)

# MCMC
fout = bsar(y ~ w + fs(x), nbasis=10, shape="Free")
fit_fout = fitted(fout)
# result
result121 = list(
  beta = data.frame(VB=vb_result$mubeta.q, MCMC=apply(apply(fit_fout$mcmc.draws$beta,2,mean))),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,1,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q
plot(x, res, main="VB: 121", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean
plot(x, res, main="MCMC: 121", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 1-2-2. y = wb + f2(x) ---------------------------------------------------
y = cbind(1,w)%*%beta + f2(x) + rnorm(N, sd=0.5)

# VB
vb_result = bsar_vb(y,x,w)

# MCMC
fout = bsar(y ~ w + fs(x), nbasis=10, shape="Free")
fit_fout = fitted(fout)

# result
result122 = list(
  beta = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,1,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q
plot(x, res, main="VB: 122", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean
plot(x, res, main="MCMC: 122", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 1-2-3. y = wb + f3(x) ---------------------------------------------------
y = cbind(1,w)%*%beta + f3(x) + rnorm(N, sd=0.5)

# VB
vb_result = bsar_vb(y,x,w)

# MCMC
fout = bsar(y ~ w + fs(x), nbasis=10, shape="Free")
fit_fout = fitted(fout)

# result
result123 = list(
  beta = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,1,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q
plot(x, res, main="VB: 123", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean
plot(x, res, main="MCMC: 123", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)


# 2. BSA-Random Intercept -------------------------------------------------
source("./RandomIntercept/normal.R")
library(bspmmGP)

# 2-1. Normal -------------------------------------------------------------
set.seed(10)
N = 80; T = 6; D = 5
x = runif(N*T)*3; ord = order(x)
u = rnorm(N)
w = matrix(rnorm(N*T*D), ncol=D)
Z = make_Z(rep(T,N))
beta = rnorm(D+1)

# 2-1-1. y = wb + zu + f1(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f1(x) + rnorm(N, sd=0.5)

# VB
vb_result = bsa(y,x,w,Z,10)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bspmr(y = ylist, x = xlist, w = wlist, nbasis = 10, random = 'gaussian', shape = 'Free', mcmc = mcmc)
saveRDS(fout, "../../simulation_rds/fout_211")
# fout = readRDS("./simulation_rds/fout_211")
# fit_fout = fitted(fout)

# result
result211 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 211", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 211", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 2-1-2. y = wb + zu + f2(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f2(x) + rnorm(N, sd=0.5)

# VB
vb_result = bsa(y,x,w,Z,10)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bspmr(y = ylist, x = xlist, w = wlist, nbasis = 10, random = 'gaussian', shape = 'Free', mcmc = mcmc)
saveRDS(fout, "../../simulation_rds/fout_212")
# fout = readRDS("./simulation_rds/fout_212")
# fit_fout = fitted(fout)

# result
result212 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 212", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 212", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 2-1-3. y = wb + zu + f3(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f3(x) + rnorm(N, sd=0.5)

# VB
vb_result = bsa(y,x,w,Z,10)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bspmr(y = ylist, x = xlist, w = wlist, nbasis = 10, random = 'gaussian', shape = 'Free', mcmc = mcmc)
saveRDS(fout, "../../simulation_rds/fout_213")
# fout = readRDS("./simulation_rds/fout_213")
# fit_fout = fitted(fout)

# result
result213 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 213", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 213", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)


# 2-2. Normal Mixture -----------------------------------------------------
set.seed(10)
N = 80; T = 6; D = 5
x = runif(N*T)*3; ord = order(x)
u = runif(N, -3, 3)
w = matrix(rnorm(N*T*D), ncol=D)
Z = make_Z(rep(T,N))
beta = rnorm(D+1)

# 2-2-1. y = wb + zu + f1(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f1(x) + rnorm(N, sd=0.5)

# VB
vb_result = mixture_bsa(y,x,w,Z,10,10)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bspmr(y = ylist, x = xlist, w = wlist, nbasis = 10, random = 'dpm', shape = 'Free', mcmc = mcmc)
saveRDS(fout, "./simulation_rds/fout_221")
# fout = readRDS("./simulation_rds/fout_221")
# fit_fout = fitted(fout)

# result
result221 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 221", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 221", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 2-2-2. y = wb + zu + f2(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f2(x) + rnorm(N, sd=0.5)

# VB
vb_result = mixture_bsa(y,x,w,Z,10,10)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bspmr(y = ylist, x = xlist, w = wlist, nbasis = 10, random = 'dpm', shape = 'Free', mcmc = mcmc)
saveRDS(fout, "./simulation_rds/fout_222")
# fout = readRDS("./simulation_rds/fout_222")
# fit_fout = fitted(fout)

# result
result222 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 222", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 222", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 2-2-3. y = wb + zu + f3(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f3(x) + rnorm(N, sd=0.5)

# VB
vb_result = mixture_bsa(y,x,w,Z,10,10)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bspmr(y = ylist, x = xlist, w = wlist, nbasis = 10, random = 'dpm', shape = 'Free', mcmc = mcmc)
saveRDS(fout, "./simulation_rds/fout_223")
# fout = readRDS("./simulation_rds/fout_223")
# fit_fout = fitted(fout)

# result
result223 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 223", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - fit_fout$wbeta$mean - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 223", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 3. BSA-Stochastic Frontier ----------------------------------------------
source("./StochasticFrontier/ex.R")
setwd("./StochasticFrontier/MCMC/")
source("bsfr_main.R")

# 3-1. Exponential --------------------------------------------------------
set.seed(10)
N = 80; T = 6; D = 5
x = runif(N*T)*3; ord = order(x)
u = rexp(N, rate=2)
w = matrix(rnorm(N*T*D), ncol=D)
Z = make_Z(rep(T,N))
beta = rnorm(D+1)


# 3-1-1. y = wb + zu + f1(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f1(x) + rnorm(N, sd=0.5)

# VB
vb_result = ex_bsa(y,x,w,Z,10,productivity=FALSE)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bsfr.main(y=ylist, w=wlist, x=xlist, nbasis=10, random='exp', mcmc=mcmc,shape='Free')
saveRDS(fout, "../../simulation_rds/fout_311")
# fout = readRDS("../../simulation_rds/fout_311")
# fit_fout = fitted.bsfm(fout)

# result
result311 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 311", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - cbind(1,w)%*%apply(fit_fout$mcmc.draws$beta,2,mean) - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 311", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 3-1-2. y = wb + zu + f2(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f2(x) + rnorm(N, sd=0.5)

# VB
vb_result = ex_bsa(y,x,w,Z,10,productivity=FALSE)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bsfr.main(y=ylist, w=wlist, x=xlist, nbasis=10, random='exp', mcmc=mcmc,shape='Free')
saveRDS(fout, "../../simulation_rds/fout_312")
# fout = readRDS("../../simulation_rds/fout_312")
# fit_fout = fitted.bsfm(fout)

# result
result312 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 312", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - cbind(1,w)%*%apply(fit_fout$mcmc.draws$beta,2,mean) - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 312", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 3-1-3. y = wb + zu + f3(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f3(x) + rnorm(N, sd=0.5)

# VB
vb_result = ex_bsa(y,x,w,Z,10,productivity=FALSE)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bsfr.main(y=ylist, w=wlist, x=xlist, nbasis=10, random='exp', mcmc=mcmc,shape='Free')
saveRDS(fout, "../../simulation_rds/fout_313")
# fout = readRDS("../../simulation_rds/fout_313")
# fit_fout = fitted.bsfm(fout)

# result
result313 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 313", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - cbind(1,w)%*%apply(fit_fout$mcmc.draws$beta,2,mean) - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 313", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 3-2. Exponenetial Mixture -----------------------------------------------
set.seed(10)
N = 80; T = 6; D = 5
x = runif(N*T)*3; ord = order(x)
u = runif(N, max=6)
w = matrix(rnorm(N*T*D), ncol=D)
Z = make_Z(rep(T,N))
beta = rnorm(D+1)

# 3-2-1. y = wb + zu + f1(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f1(x) + rnorm(N, sd=0.5)

# VB
vb_result = exmixture_bsa(y,x,w,Z,10,productivity=FALSE)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bsfr.main(y=ylist, w=wlist, x=xlist, nbasis=10, random='dpm', mcmc=mcmc,shape='Free')
saveRDS(fout, "../../simulation_rds/fout_321")
# fout = readRDS("../../simulation_rds/fout_321")
# fit_fout = fitted.bsfm(fout)

# result
result321 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 321", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - cbind(1,w)%*%apply(fit_fout$mcmc.draws$beta,2,mean) - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 321", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 3-2-2. y = wb + zu + f2(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f2(x) + rnorm(N, sd=0.5)

# VB
vb_result = exmixture_bsa(y,x,w,Z,10,productivity=FALSE)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bsfr.main(y=ylist, w=wlist, x=xlist, nbasis=10, random='dpm', mcmc=mcmc,shape='Free')
saveRDS(fout, "../../simulation_rds/fout_322")
# fout = readRDS("../../simulation_rds/fout_322")
# fit_fout = fitted.bsfm(fout)

# result
result322 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 322", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - cbind(1,w)%*%apply(fit_fout$mcmc.draws$beta,2,mean) - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 322", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)

# 3-2-3. y = wb + zu + f3(x) ----------------------------------------------
y = cbind(1,w)%*%beta + Z%*%u + f3(x) + rnorm(N, sd=0.5)

# VB
vb_result = exmixture_bsa(y,x,w,Z,10,productivity=FALSE)

# MCMC
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
  ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
  xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
  wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
fout = bsfr.main(y=ylist, w=wlist, x=xlist, nbasis=10, random='dpm', mcmc=mcmc,shape='Free')
saveRDS(fout, "../../simulation_rds/fout_323")
# fout = readRDS("../../simulation_rds/fout_323")
# fit_fout = fitted.bsfm(fout)

# result
result323 = list(
  fixed = data.frame(VB=vb_result$mubeta.q, MCMC=apply(fit_fout$mcmc.draws$beta,2,mean)),
  random = data.frame(VB=vb_result$muu.q, MCMC=apply(fit_fout$mcmc.draws$b,2,mean)),
  theta = data.frame(VB=vb_result$mutheta.q, MCMC=apply(fit_fout$mcmc.draws$theta,2,mean)[-1]))
par(mfrow=c(1,2),mar=c(2,3.9,3,1))
res = y - cbind(1,w)%*%vb_result$mubeta.q - Z%*%vb_result$muu.q
plot(x, res, main="VB: 323", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)
res = y - cbind(1,w)%*%apply(fit_fout$mcmc.draws$beta,2,mean) - Z%*%fit_fout$bhat$mean
plot(x, res, main="MCMC: 323", xlab="", ylab="residual")
lines(fit_fout$xgrid, fit_fout$fxgrid$mean, lwd=3, col=2)
lines(fit_fout$xgrid, fit_fout$fxgrid$lower, lwd=2, col=3)
lines(fit_fout$xgrid, fit_fout$fxgrid$upper, lwd=2, col=3)
