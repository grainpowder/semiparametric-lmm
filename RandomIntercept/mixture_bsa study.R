source("../RandomIntercept/normal.R")
source("../misc/make_Z.R")
library(ggplot2)
library(gridExtra)


# Simulation Study --------------------------------------------------------
# Simulation size setting
set.seed(50); R = 10 # R : Truncation level
N = 50; T = 5; D = 8
Z = make_Z(rep(T, N))

# Data generation
beta = rnorm(D+1)
w = matrix(rnorm(N*T*D), N*T, D)
f = function(x) 2*x*sin(pi*x)
x = 3*runif(N*T); ord = order(x)
# Represent random intercept as mixture of 3 normal distributions
assigner = runif(N)
mu1 = 3; mu2 = 0; mu3 = -3
u = rnorm(N, mu1, 0.36)
u[assigner > 0.33] = rnorm(sum(assigner > 0.33), mu2, 0.36)
u[assigner > 0.66] = rnorm(sum(assigner > 0.66), mu3, 0.36)
y = cbind(1, w)%*%beta + Z%*%u + f(x) + rnorm(nrow(Z))

# Estimation
start = as.numeric(Sys.time())
result = mixture_bsa(y,x,w,Z,23,R,tol=1e-6) # truncate at R=10
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
  residual = y-(cbind(1,w)%*%result$mubeta.q + Z%*%result$muu.q))
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


# Cadmium -----------------------------------------------------------------
library(bspmmGP)
library(dplyr)
library(Matrix)
library(latex2exp)

study = cadmium$studycode   # study code
w = cbind(ifelse(cadmium$ethnicity == 1, 1, 0),
          ifelse(cadmium$age >= 50, 1, 0),
          ifelse(cadmium$gender == 1, 1, 0),
          ifelse(cadmium$gender == 0, 1, 0))
colnames(w)=c('asian','age>=50','female','male')
y = log(cadmium$b2_GM)
x = log(cadmium$Ucd_GM)


# VB estimation
# sort the data w.r.t. studycode, to make correspondence with MCMC
cadmium_bind = cbind.data.frame(study,w,y,x) %>% arrange(study)
study = cadmium_bind$study
y = cadmium_bind$y
x = cadmium_bind$x; ord = order(x)
w = cadmium_bind %>% dplyr::select(-study,-y,-x) %>% as.matrix()
label = unique(study)
Z = matrix(0, nrow(w), length(label))
for(idx in 1:length(label)) Z[which(study==label[idx]),idx] = 1

# Emitting result
start = as.numeric(Sys.time())
vbresult_dpm = mixture_bsa(y,x,w,Z,30,10)
print(paste("VB", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))


# MCMC estimation
study = cadmium$studycode   # study code
ustudy = names(table(study))
wdata = cbind(ifelse(cadmium$ethnicity == 1, 1, 0),
              ifelse(cadmium$age >= 50, 1, 0),
              ifelse(cadmium$gender == 1, 1, 0),
              ifelse(cadmium$gender == 0, 1, 0))
colnames(wdata)=c('asian','age>=50','female','male')
ylist = xlist = wlist = list()
for(i in ustudy){
  indices = which(i==study)
  ylist[[i]] = log(cadmium$b2_GM[indices])  
  xlist[[i]] = log(cadmium$Ucd_GM[indices]) 
  wlist[[i]] = wdata[indices,]
}
# mcmc parameters
mcmc = list(nblow0 = 10000, nblow = 100000, nskip = 50,
            smcmc = 2000, maxmodmet = 10, ndisp = 2000)
start = as.numeric(Sys.time())
foutr_dpm = bspmr(y = ylist, x = xlist, w = wlist, nbasis = 50, random = 'dpm', shape = 'Free', mcmc = mcmc)
ffitr_dpm = fitted(foutr_dpm)
print(paste("MCMC:", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))

means = c(
  vbresult_dpm$mubeta.q,
  apply(ffitr_dpm$mcmc.draws$beta, 2, mean)
)
upper = c(
  qnorm(0.025, vbresult_dpm$mubeta.q, sqrt(diag(vbresult_dpm$sigbeta.q))),
  apply(ffitr_dpm$mcmc.draws$beta, 2, quantile, 0.025)
)
lower = c(
  qnorm(0.975, vbresult_dpm$mubeta.q, sqrt(diag(vbresult_dpm$sigbeta.q))),
  apply(ffitr_dpm$mcmc.draws$beta, 2, quantile, 0.975)
)
fixed = data.frame(
  names = factor(c("intercept",dimnames(w)[[2]]),levels=rev(c("intercept",dimnames(w)[[2]]))),
  mean = means,
  upper = upper,
  lower = lower,
  method = factor(rep(c("VB","MCMC"), each=length(means)/2), levels=c("MCMC","VB"))
)
random = data.frame(
  MCMC = apply(ffitr_dpm$mcmc.draws$b, 2, mean), 
  VB = vbresult_dpm$muu.q
)

# Fixed effects
fe = ggplot(fixed, aes(x=names, y=mean, color=method)) + 
  geom_errorbar(data=fixed, mapping=aes(ymin=lower, ymax=upper), position="dodge") +
  geom_point(position=position_dodge(width=0.9), size=1.6) +
  geom_hline(yintercept = 0, color="black", size=0.7, linetype="dotted") +
  scale_color_manual(values=c("deepskyblue4","darkorange3")) +
  scale_x_discrete(
    breaks=c("male","female","age>=50","asian","intercept"),
    labels=c("male","female",TeX("$age\\geq 50$"),"asian","intercept")) +
  xlab("") + ylab("") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=11),
        panel.grid = element_blank(),
        legend.position = c(0.83,0.15)) +
  ggtitle("Fixed effects")

# Random effects
re = ggplot(random, aes(x=MCMC, y=VB)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  theme_classic() +
  ggtitle("Random effects")


# Nonlinear curve
studyId = factor(study, labels = 1:length(table(study)))
mcmc_polygon = data.frame(
  variate = c(ffitr_dpm$xgrid, rev(ffitr_dpm$xgrid)),
  area = c(ffitr_dpm$fxgrid$upper, rev(ffitr_dpm$fxgrid$lower)))
vb_polygon = data.frame(
  variate = c(x[ord], rev(x[ord])),
  area = c(vbresult_dpm$post_upper[ord], rev(vbresult_dpm$post_lower[ord])))
residual = data.frame(
  variate=rep(x,2),
  residual=c(unlist(y)-ffitr_dpm$wbeta$mean-ffitr_dpm$bhat$mean[studyId],
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
  xlab("Ucd_GM")+
  ylab("Residuals")+
  theme_bw()+
  ggtitle("Mean curve") +
  theme(legend.position = c(0.15,0.85),
        panel.grid = element_blank())

# Export
png("./figures/cadmium.png",width=900,height=320)
grid.arrange(fe, re, mc, nrow=1)
dev.off()