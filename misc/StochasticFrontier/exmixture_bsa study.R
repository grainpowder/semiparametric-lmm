source("../StochasticFrontier/ex.R")
source("../misc/make_Z.R")
library(matrixStats)
library(truncnorm)
library(ggplot2)
library(gridExtra)

# Simulation study --------------------------------------------------------

# Simulation size setting
set.seed(5); R = 10 # R : Truncation level
N = 50; T = 4; D = 8
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

# Estimation
start = as.numeric(Sys.time())
result = exmixture_bsa(y,x,w,Z,5,R) # truncate at R=10
print(paste(round(as.numeric(Sys.time())-start, 4), "seconds elapsed."))

# Result visualization
# Fixed effect(black point: true value, golden line: corresponding 95% credible interval)
beta_ord = order(beta[-1])
ggdata = data.frame(
        upper = qnorm(0.975,mean=result$mubeta.q[-1], sd=sqrt(diag(result$sigbeta.q)[-1]))[beta_ord],
        lower = qnorm(0.025,mean=result$mubeta.q[-1], sd=sqrt(diag(result$sigbeta.q)[-1]))[beta_ord],
        beta = beta[-1][beta_ord],
        index = 1:D)
fe = ggplot() +
        geom_errorbar(data=ggdata, aes(x=index, ymin=lower, ymax=upper), size=1.1, color="darkorange3", width=0.3) +
        geom_point(data=ggdata, aes(x=index, y=beta), size=1.1) +
        theme_bw() +
        xlab("Arranged index") +
        ylab("") +
        ggtitle("Fixed effects")

# Random effect(black point: true value, golden line: corresponding 95% credible interval)
u_ord = order(u)
ggdata = data.frame(
        upper = qnorm(0.975,mean=result$muu.q, sd=sqrt(result$sigu.q))[u_ord],
        lower = qnorm(0.025,mean=result$muu.q, sd=sqrt(result$sigu.q))[u_ord],
        u = u[u_ord],
        index = 1:N)
re = ggplot() +
        geom_errorbar(data=ggdata, aes(x=index, ymin=lower, ymax=upper), size=1.1, color="darkorange3", width=0.3) +
        geom_point(data=ggdata, aes(x=index, y=u), size=1.1) +
        theme_bw() +
        xlab("Arranged index") +
        ylab("") +
        ggtitle("Random effects")

# Nonparametric(black point: residual, golden line: mean curve)
vb_polygon = data.frame(
        variate = c(x[ord], rev(x[ord])),
        area = c(result$post_upper[ord], rev(result$post_lower[ord])))
residual = data.frame(
        variate = x,
        residual = y-(cbind(1,w)%*%result$mubeta.q - Z%*%result$muu.q))
mc = ggplot() +
        geom_point(data=residual,
                   mapping=aes(x=variate,y=residual),size=1.1) +
        geom_polygon(data=vb_polygon,aes(x=variate,y=area), fill="darkorange3", alpha=0.3) +
        geom_line(data=data.frame(x=x[ord], y=result$post_curve[ord]),
                  mapping=aes(x=x,y=y),color="darkorange3",size=1.1) +
        xlab("Nonparametric variate")+
        ylab("Residuals")+
        theme_bw() +
        ggtitle("Mean curve")

# Export
png("./figures/simulation.png",width=900,height=320)
grid.arrange(fe, re, mc, nrow=1)
dev.off()

# University --------------------------------------------------------------
library(dplyr)
data = read.csv("./univ.csv")
data = data %>% 
        mutate(id=factor(data$univ_id,labels=1:80)) %>% 
        filter(year >= 7) %>% 
        select(univ_id,lneduc,lngrad,lnunder,lnw_ac,lnw_nac,
                       d_nohos,d_sci,d_arts,d_med,d_educ)
nobs  = nrow(data)
id    = as.numeric(data$id)
nid   = length(unique(data$univ_id))
y = data$lneduc
x = data$lnunder; ord = order(x)
w = data %>% select(lngrad, lnw_ac, lnw_nac, d_nohos, d_sci, d_arts, d_med, d_educ) %>% as.matrix()
Z = data %>% group_by(univ_id) %>% summarize(count=n()) %>% select(count) %>% unlist() %>% make_Z()

# VB Estimation
start = as.numeric(Sys.time())
vbresult_dpm = exmixture_bsa(y,x,w,Z,J=10,productivity=FALSE)
print(paste(round(as.numeric(Sys.time())-start, 4), "seconds elapsed."))

# MCMC estimation
source("./MCMC/bsfr_main.R")
setwd("./MCMC/")
set.seed(10)
csrows = cumsum(apply(Z,2,sum))
csrows = c(0,csrows)
ylist = xlist = wlist = list()
for(i in 1:ncol(Z)){
        ylist[[i]] = y[(csrows[i]+1):csrows[i+1]]
        xlist[[i]] = x[(csrows[i]+1):csrows[i+1]]
        wlist[[i]] = w[(csrows[i]+1):csrows[i+1],]}
# MCMC parameters
mcmc = list(nblow=10000,nskip=5,smcmc=20000,ndisp=5000)
start = as.numeric(Sys.time())
foutr_dpm = bsfr.main(y=ylist,w=wlist,x=xlist,nbasis=10,random='dpm',mcmc=mcmc,shape='Free')
print(paste("MCMC:", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))
beta0 = coef.bsfm(foutr_dpm)$beta[1]
ffitr_dpm = fitted.bsfm(foutr_dpm)
res_MC = y-beta0-ffitr_dpm$wbeta$mean-Z%*%ffitr_dpm$bhat$mean

# Values to compare through visualization
vb_fixed = vbresult_dpm$mubeta.q[-1]
mc_fixed = apply(ffitr_dpm$mcmc.draws$beta,2,mean)[-1]
means = c(vb_fixed, mc_fixed)
upper = c(qnorm(0.025, vb_fixed, sqrt(diag(vbresult_dpm$sigbeta.q)[-1])), apply(ffitr_dpm$mcmc.draws$beta, 2, quantile, 0.025)[-1])
lower = c(qnorm(0.975, vb_fixed, sqrt(diag(vbresult_dpm$sigbeta.q)[-1])), apply(ffitr_dpm$mcmc.draws$beta, 2, quantile, 0.975)[-1])
fixed = data.frame(
        names = factor(dimnames(w)[[2]],levels=rev(dimnames(w)[[2]])),
        mean = means, upper = upper, lower = lower,
        method = factor(rep(c("VB","MCMC"), each=length(means)/2), levels=c("MCMC","VB")))
random = data.frame(
        MCMC = apply(ffitr_dpm$mcmc.draws$b, 2, mean), 
        VB = vbresult_dpm$muu.q)

# Fixed effects
fe = ggplot(fixed, aes(x=names, y=mean, color=method)) + 
        geom_errorbar(data=fixed, mapping=aes(ymin=lower, ymax=upper), position="dodge") +
        geom_point(position=position_dodge(width=0.9), size=1.6) +
        geom_hline(yintercept = 0, color="black", size=0.7, linetype="dotted") +
        scale_color_manual(values=c("deepskyblue4","darkorange3")) +
        xlab("") + ylab("") +
        coord_flip() +
        theme_bw() +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size=11),
              panel.grid = element_blank(),
              legend.position = c(0.17,0.83)) +
        ggtitle("Fixed effects")

# Random effects
re = ggplot(random, aes(x=MCMC, y=VB)) +
        geom_point() +
        geom_abline(slope=1, intercept=0) +
        theme_classic() +
        ggtitle("Random effects")

# Nonlinear curve
mcmc_polygon = data.frame(
        variate = c(ffitr_dpm$xgrid, rev(ffitr_dpm$xgrid)),
        area = c(ffitr_dpm$fxgrid$upper, rev(ffitr_dpm$fxgrid$lower)))
vb_polygon = data.frame(
        variate = c(x[ord], rev(x[ord])),
        area = c(vbresult_dpm$post_upper[ord], rev(vbresult_dpm$post_lower[ord])))
residual = data.frame(
        variate=rep(x,2),
        residual=c(unlist(y)-ffitr_dpm$wbeta$mean-ffitr_dpm$bhat$mean[id],
                   y-cbind(1,w)%*%vbresult_dpm$mubeta.q-Z%*%vbresult_dpm$muu.q),
        method=rep(c("MCMC","VB"),each=length(x)))
mc = ggplot() +
        geom_point(data=residual,
                   mapping=aes(x=variate,y=residual,color=method,shape=method)) +
        scale_color_manual(values=c("deepskyblue4","darkorange3"))+
        geom_polygon(data=mcmc_polygon,aes(x=variate,y=area), fill="deepskyblue4", alpha=0.3) +
        geom_polygon(data=vb_polygon,aes(x=variate,y=area), fill="darkorange3", alpha=0.3) +
        geom_line(data=data.frame(x=ffitr_dpm$xgrid, y=ffitr_dpm$fxgrid$mean),
                  mapping=aes(x=x,y=y),color="deepskyblue4",size=1.5) +
        geom_line(data=data.frame(x=x[ord], y=vbresult_dpm$post_curve[ord]),
                  mapping=aes(x=x,y=y),color="darkorange3",size=1.5) +
        xlab("lnunder")+
        ylab("Residuals")+
        theme_bw()+
        ggtitle("Mean curve") +
        theme(legend.position = c(0.85,0.15),
              panel.grid = element_blank())

# Export
png("../figures/university.png",width=900,height=320)
grid.arrange(fe, re, mc, nrow=1)
dev.off()