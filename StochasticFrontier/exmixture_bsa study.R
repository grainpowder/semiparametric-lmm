source("../StochasticFrontier/ex.R")
source("../misc/make_Z.R")
library(matrixStats)
library(truncnorm)

# Simulation study --------------------------------------------------------

# Simulation size setting
set.seed(1); R = 10 # R : Truncation level
N = 80; T = 12; D = 8
csrows = rep(T, N)
csrows[round(N/4)] = T-2
Z = make_Z(rep(T, N))

# Data generation
beta = rnorm(D+1)
w = matrix(rnorm(N*T*D), N*T, D)*4
f = function(x) 0.5*exp(x)
x = 3*runif(N*T); ord = order(x)
# Represent inefficiency as mixture of 2 uniform distributions
assigner = runif(N)
u1 = runif(N)
u2 = 3+runif(N)
u = c(u1[assigner > 0.5], u2[assigner <= 0.5])
u = u[sample(1:N,N)] # Randomly shuffle the mixed inefficiencies
rm("u1","u2")
y = cbind(1, w)%*%beta - Z%*%u + f(x) + rnorm(nrow(Z))

# VB
start = as.numeric(Sys.time())
result = exmixture_bsa(y,x,w,Z,23,R,eps=1e-4) # truncate at R=10
print(paste(round(as.numeric(Sys.time())-start, 4), "seconds elapsed."))

# VB results
# Manually exported as 800 * 320 in png type
par(mar=c(4,3,2,1), mfrow=c(1,2))
# VB: Fixed effect
beta_ord = order(beta[-1])
upper = qnorm(0.975,mean=result$mubeta.q[-1], sd=sqrt(diag(result$sigbeta.q)[-1]))
lower = qnorm(0.025,mean=result$mubeta.q[-1], sd=sqrt(diag(result$sigbeta.q)[-1]))
plot(1:D, beta[-1][beta_ord],
     xlab="index",ylab="",
     ylim=c(min(c(lower, beta)),max(c(upper, beta))),
     pch=19, col=2,
     main="Fixed effect(VB)")
for (idx in 1:N) lines(c(idx,idx), c(upper[beta_ord][idx],lower[beta_ord][idx]))

# VB: Random effect
u_ord = order(u)
upper = qtruncnorm(0.975,a=0,mean=result$muu.q, sd=sqrt(result$sigu.q))
lower = qtruncnorm(0.025,a=0,mean=result$muu.q, sd=sqrt(result$sigu.q))
plot(1:N, u[u_ord],xlab="idx",ylab="",
     pch=19,col=2,
     ylim=c(min(c(0,lower,u)),max(c(upper, u))),
     main="Inefficiency(VB)")
for (idx in 1:N) lines(c(idx,idx), c(upper[u_ord][idx],lower[u_ord][idx]))

# VB: Nonparametric
# Manually exported as 600 * 400 in png type
par(mfrow=c(1,1))
res = y-(cbind(1,w)%*%result$mubeta.q - Z%*%result$muu.q)
plot(x,res, main="Fitted mean curve(VB)", ylab="")
lines(x[ord],result$post_curve[ord],lwd=3,col=2)
lines(x[ord],result$post_upper[ord],lwd=2,lty=2)
lines(x[ord],result$post_lower[ord],lwd=2,lty=2)

# MCMC
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

# MCMC results
# Manually exported as 800 * 320 in png type
par(mar=c(4,3,2,1), mfrow=c(1,2))
# MCMC: Fixed effect
upper = apply(mcmc_fit$mcmc.draws$beta, 2, quantile, 0.975)[-1]
lower = apply(mcmc_fit$mcmc.draws$beta, 2, quantile, 0.025)[-1]
plot(1:D, beta[-1][beta_ord],
     xlab="index",ylab="",
     ylim=c(min(c(lower, beta)),max(c(upper, beta))),
     pch=19, col=2,
     main="Fixed effect(MCMC)")
for (idx in 1:N) lines(c(idx,idx), c(upper[beta_ord][idx],lower[beta_ord][idx]))

# MCMC: Random effect
upper = rev(apply(mcmc_fit$mcmc.draws$b, 2, quantile, 0.975)[u_ord])
lower = rev(apply(mcmc_fit$mcmc.draws$b, 2, quantile, 0.025)[u_ord])
plot(1:N, u[u_ord],xlab="idx",ylab="",
     pch=19,col=2,
     ylim=c(min(c(0,lower,u)),max(c(upper, u))),
     main="Inefficiency(MCMC)")
for (idx in 1:N) lines(c(idx,idx), c(upper[idx],lower[idx]))

# MCMC: Nonparametric
# Manually exported as 600 * 400 in png type
par(mfrow=c(1,1))
plot(x,res_MC, main="Fitted mean curve(MCMC)", ylab="")
lines(mcmc_fit$xgrid,mcmc_fit$fxgrid$mean,lwd=3,col=2)
lines(mcmc_fit$xgrid,mcmc_fit$fxgrid$lower,lwd=2,lty=2)
lines(mcmc_fit$xgrid,mcmc_fit$fxgrid$upper,lwd=2,lty=2)


# University --------------------------------------------------------------

library(dplyr)
data = read.csv("../univ.csv")
data = data %>% select(univ_id,lneduc,lngrad,lnunder,lnw_ac,lnw_nac,
                       d_nohos,d_sci,d_arts,d_med,d_educ,lnkaken,p_sub)%>% 
        mutate(id=factor(data$univ_id,labels=1:80)) %>% 
        mutate(year=c(rep(1:12,35),1:10,rep(1:12,44))) # 36th obs has 10 years data, others have 12 years
head(data)

nobs  = nrow(data)
id    = as.numeric(data$id)
nid   = length(unique(data$univ_id))
y = data$lneduc
x = data$lnunder; ord = order(x)
w = data %>% select(lngrad, lnw_ac, lnw_nac, d_nohos, d_sci, d_arts, d_med, d_educ) %>% as.matrix()
head(w)
Z = data %>% 
        group_by(id) %>% 
        summarize(count=n()) %>% 
        select(count) %>% 
        unlist() %>% 
        make_Z()

# University: VB ----------------------------------------------------------

start = as.numeric(Sys.time())
result = ex_bsa(y,x,w,Z,J=2,productivity=FALSE)
print(paste(round(as.numeric(Sys.time())-start, 4), "seconds elapsed."))

# VB: Fixed effect
mean_fixed = result$mubeta.q[-1]; beta_ord=order(mean_fixed)
upper = qnorm(0.975,mean=result$mubeta.q[-1], sd=sqrt(diag(result$sigbeta.q)[-1]))[beta_ord]
lower = qnorm(0.025,mean=result$mubeta.q[-1], sd=sqrt(diag(result$sigbeta.q)[-1]))[beta_ord]
plot(1:ncol(w), mean_fixed[beta_ord],
     xlab="index",ylab="",
     ylim=c(min(c(lower, mean_fixed)),max(c(upper, mean_fixed))),
     pch=19, col=2,
     main="Fixed effect(VB)")
for (idx in 1:ncol(w)) lines(c(idx,idx), c(upper[idx],lower[idx]))

# Random effect
mean_random = result$muu.q; u_ord=order(mean_random)
upper = qtruncnorm(0.975,a=0,mean=result$muu.q, sd=sqrt(result$sigu.q))[u_ord]
lower = qtruncnorm(0.025,a=0,mean=result$muu.q, sd=sqrt(result$sigu.q))[u_ord]
plot(1:ncol(Z),mean_random[u_ord],xlab="idx",ylab="",
     pch=19,col=2,
     ylim=c(min(c(0,lower,mean_random)),max(c(upper, mean_random))),
     main="Inefficiency(VB)")
for (idx in 1:ncol(Z)) lines(c(idx,idx), c(upper[idx],lower[idx]))
        
# Nonparametric
res = y-(cbind(1,w)%*%result$mubeta.q + Z%*%result$muu.q)
plot(x,res, main="Fitted mean curve(VB)", ylab="")
lines(x[ord],result$post_curve[ord],lwd=3,col=2)

# University: MCMC --------------------------------------------------------
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
        ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
        xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
        wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
start = as.numeric(Sys.time())
mcmcresult = bsfr.main(y=ylist,w=wlist,x=xlist,nbasis=2,random='exp',mcmc=mcmc,shape='Free')
print(paste("MCMC:", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))
beta0 = coef.bsfm(mcmcresult)$beta[1]
mcmc_fit = fitted.bsfm(mcmcresult)
res_MC = y-beta0-mcmc_fit$wbeta$mean-Z%*%mcmc_fit$bhat$mean

# MCMC results
# Manually exported as 800 * 320 in png type
par(mar=c(4,3,2,1), mfrow=c(1,2))
# MCMC: Fixed effect
mean_fixed = apply(mcmc_fit$mcmc.draws$beta,2,mean)[-1]; beta_ord=order(mean_fixed)
upper = apply(mcmc_fit$mcmc.draws$beta, 2, quantile, 0.975)[-1][beta_ord]
lower = apply(mcmc_fit$mcmc.draws$beta, 2, quantile, 0.025)[-1][beta_ord]
plot(1:ncol(w), mean_fixed[beta_ord],
     xlab="index",ylab="",
     ylim=c(min(c(lower, mean_fixed)),max(c(upper, mean_fixed))),
     pch=19, col=2,
     main="Fixed effect(MCMC)")
for (idx in 1:ncol(w)) lines(c(idx,idx), c(upper[idx],lower[idx]))

# MCMC: Random effect
mean_random = apply(mcmc_fit$mcmc.draws$b,2,mean); u_ord=order(mean_random)
upper = apply(mcmc_fit$mcmc.draws$b, 2, quantile, 0.975)[u_ord]
lower = apply(mcmc_fit$mcmc.draws$b, 2, quantile, 0.025)[u_ord]
plot(1:ncol(Z),mean_random[u_ord],xlab="idx",ylab="",
     pch=19,col=2,
     ylim=c(min(c(0,lower,mean_random)),max(c(upper, mean_random))),
     main="Inefficiency(MCMC)")
for (idx in 1:ncol(Z)) lines(c(idx,idx), c(upper[idx],lower[idx]))

# MCMC: Nonparametric
# Manually exported as 600 * 400 in png type
par(mfrow=c(1,1))
upper=mcmc_fit$fxgrid$upper
lower=mcmc_fit$fxgrid$lower
plot(x,res_MC,ylim=c(min(res_MC, lower), max(res_MC,upper)),main="Fitted mean curve(MCMC)", ylab="")
lines(mcmc_fit$xgrid,mcmc_fit$fxgrid$mean,lwd=3,col=2)
lines(mcmc_fit$xgrid,lower,lwd=2,lty=2)
lines(mcmc_fit$xgrid,upper,lwd=2,lty=2)
