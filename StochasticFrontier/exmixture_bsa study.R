source("../StochasticFrontier/ex.R")
source("../misc/make_Z.R")
library(matrixStats)
library(truncnorm)


# Data generation ---------------------------------------------------------
# Simulation size setting
set.seed(1); R = 10 # R : Truncation level
N = 50; T = 4; D = 8
Z = make_Z(rep(T, N))

# Data generation
beta = rnorm(D+1)
w = matrix(rnorm(N*T*D), N*T, D)
f = function(x) -3*x*sin(pi*x)
x = 3*runif(N*T); ord = order(x)
# Represent inefficiency as mixture of 2 uniform distributions
assigner = runif(N)
u1 = runif(N)
u2 = 3+runif(N)
u = c(u1[assigner > 0.5], u2[assigner <= 0.5])
u = u[sample(1:N,N)] # Randomly shuffle the mixed inefficiencies
rm("u1","u2")
y = cbind(1, w)%*%beta - Z%*%u + f(x) + rnorm(nrow(Z))

# Estimation:VB -----------------------------------------------------------
# Result emitting
start = as.numeric(Sys.time())
result = exmixture_bsa(y,x,w,Z,23,R,eps=1e-4) # truncate at R=10
print(paste(round(as.numeric(Sys.time())-start, 4), "seconds elapsed."))

# Visualizing the result --------------------------------------------------
# Manually exported as 800 * 320 in png type
par(mar=c(4,3,2,1), mfrow=c(1,2))
# Fixed effect
beta_ord = order(beta[-1])
upper = qnorm(0.975,mean=result$mubeta.q[-1], sd=sqrt(diag(result$sigbeta.q)[-1]))
lower = qnorm(0.025,mean=result$mubeta.q[-1], sd=sqrt(diag(result$sigbeta.q)[-1]))
plot(1:D, beta[-1][beta_ord],
     xlab="index",ylab="",
     ylim=c(min(c(lower, beta)),max(c(upper, beta))),
     pch=19, col=2,
     main="Fixed effect(VB)")
for (idx in 1:N) lines(c(idx,idx), c(upper[beta_ord][idx],lower[beta_ord][idx]))

# Random effect
u_ord = order(u)
upper = qtruncnorm(0.975,a=0,mean=result$muu.q, sd=sqrt(result$sigu.q))
lower = qtruncnorm(0.025,a=0,mean=result$muu.q, sd=sqrt(result$sigu.q))
plot(1:N, u[u_ord],xlab="idx",ylab="",
     pch=19,col=2,
     ylim=c(min(c(0,lower,u)),max(c(upper, u))),
     main="Inefficiency(VB)")
for (idx in 1:N) lines(c(idx,idx), c(upper[u_ord][idx],lower[u_ord][idx]))

# Nonparametric
# Manually exported as 600 * 400 in png type
par(mfrow=c(1,1))
res = y-(cbind(1,w)%*%result$mubeta.q - Z%*%result$muu.q)
plot(x,res, main="Fitted mean curve(VB)", ylab="")
lines(x[ord],result$post_curve[ord],lwd=3,col=2)
lines(x[ord],result$post_upper[ord],lwd=2,lty=2)
lines(x[ord],result$post_lower[ord],lwd=2,lty=2)

# Estimation: MCMC --------------------------------------------------------
# Result emitting

csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
        ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
        xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
        wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
source("./MCMC/bsfr_main.R")
setwd("./MCMC/")
mcmc = list(nblow=1000,nskip=5,smcmc=2000,ndisp=5000)
start = as.numeric(Sys.time())
mcmcresult = bsfr.main(y=ylist,w=wlist,x=xlist,nbasis=5,random='dpm',mcmc=mcmc,shape='Free')
print(paste("MCMC:", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))
beta0 = coef.bsfm(mcmcresult)$beta[1]
mcmc_fit = fitted.bsfm(mcmcresult)
res_MC = y-beta0-mcmc_fit$wbeta$mean-Z%*%mcmc_fit$bhat$mean

# Visualizing the result --------------------------------------------------
# Manually exported as 800 * 320 in png type
par(mar=c(4,3,2,1), mfrow=c(1,2))
# Fixed effect
upper = apply(mcmc_fit$mcmc.draws$beta, 2, quantile, 0.975)[-1]
lower = apply(mcmc_fit$mcmc.draws$beta, 2, quantile, 0.025)[-1]
plot(1:D, beta[-1][beta_ord],
     xlab="index",ylab="",
     ylim=c(min(c(lower, beta)),max(c(upper, beta))),
     pch=19, col=2,
     main="Fixed effect(MCMC)")
for (idx in 1:N) lines(c(idx,idx), c(upper[beta_ord][idx],lower[beta_ord][idx]))

# Random effect
upper = rev(apply(mcmc_fit$mcmc.draws$b, 2, quantile, 0.975)[u_ord])
lower = rev(apply(mcmc_fit$mcmc.draws$b, 2, quantile, 0.025)[u_ord])
plot(1:N, u[u_ord],xlab="idx",ylab="",
     pch=19,col=2,
     ylim=c(min(c(0,lower,u)),max(c(upper, u))),
     main="Inefficiency(MCMC)")
for (idx in 1:N) lines(c(idx,idx), c(upper[idx],lower[idx]))

# Nonparametric
# Manually exported as 600 * 400 in png type
par(mfrow=c(1,1))
plot(x,res_MC, main="Fitted mean curve(MCMC)", ylab="")
lines(mcmc_fit$xgrid,mcmc_fit$fxgrid$mean,lwd=3,col=2)
lines(mcmc_fit$xgrid,mcmc_fit$fxgrid$lower,lwd=2,lty=2)
lines(mcmc_fit$xgrid,mcmc_fit$fxgrid$upper,lwd=2,lty=2)

# # Hospital ----------------------------------------------------------------
# library(coda)
# library(dplyr)
# library(ggplot2)
# hospital = read.csv("./MCMC/hospital.csv", header=TRUE)
# Z = make_Z(as.numeric(table(hospital$index)))
# y = unlist(hospital$ln_cost)
# w = hospital %>% select(-ln_cost, -ln_wage, -index) %>% as.matrix()
# x = rep(1:5, nrow(w)/5); ord = order(x)
# 
# 
# # Hospital: VB ------------------------------------------------------------
# start = as.numeric(Sys.time())
# vb_fit = exmixture_bsa(y,x,w,Z,5,productivity=FALSE)
# print(paste("VB:", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))
# res_VB = y - drop(cbind(1,w)%*%vb_fit$mubeta.q+Z%*%vb_fit$muu.q)
# 
# # Hospital: MCMC ----------------------------------------------------------
# csrows = cumsum(apply(Z,2,sum))
# csrows = c(0,csrows)
# ylist = xlist = wlist = list()
# for(i in 1:ncol(Z)){
#         ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
#         xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
#         wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
# source("./MCMC/bsfr_main.R")
# setwd("./MCMC/")
# mcmc = list(nblow=1000,nskip=5,smcmc=2000,ndisp=5000)
# start = as.numeric(Sys.time())
# mcmcresult = bsfr.main(y=ylist,w=wlist,x=xlist,nbasis=5,random='dpm',mcmc=mcmc,shape='Free')
# print(paste("MCMC:", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))
# beta0 = coef.bsfm(mcmcresult)$beta[1]
# mcmc_fit = fitted.bsfm(mcmcresult)
# res_MC = y-beta0-mcmc_fit$wbeta$mean-Z%*%mcmc_fit$bhat$mean
# 
# # MCMC is not mixed well in this dataset