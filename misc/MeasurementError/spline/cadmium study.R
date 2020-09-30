library(bspmmGP)
library(dplyr)


# VB ----------------------------------------------------------------------
source("randint.R")
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
vb_result = randint(y,w,x,Z)
print(paste("VB", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))

start = as.numeric(Sys.time())
dpm_result = randint_dpm(y,w,x,Z)
print(paste("VB", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))

plot(vb_result$mubeta.q, dpm_result$mubeta.q,xlab="Normal",ylab="DPM",main="Fixed Effects")
lines(-10:10,-10:10)
plot(vb_result$mub.q, dpm_result$mub.q,xlab="Normal",ylab="DPM",main="Random Intercepts")
lines(-10:10,-10:10)
plot(vb_result$xgrid,vb_result$post_curve,type="l")
lines(dpm_result$xgrid,dpm_result$post_curve,type="l",col=2)

# MCMC --------------------------------------------------------------------
library(rstan)
load("Spline_Stan_v3.RData")

betaHat =apply(extract(ps.spline,'beta')$beta,2,mean)
wbm = W%*%betaHat
b = colMeans(extract(ps.spline,'b')$b)


# Report ------------------------------------------------------------------
means = c(
  vb_result$mubeta.q,
  apply(extract(ps.spline,"beta")$beta, 2, mean)
)
upper = c(
  qnorm(0.025, vb_result$mubeta.q, sqrt(diag(vb_result$sigbeta.q))),
  apply(extract(ps.spline,"beta")$beta, 2, quantile, 0.025)
)
lower = c(
  qnorm(0.975, vb_result$mubeta.q, sqrt(diag(vb_result$sigbeta.q))),
  apply(extract(ps.spline,"beta")$beta, 2, quantile, 0.975)
)
fixed = data.frame(
  names = factor(c("intercept",dimnames(w)[[2]]),levels=rev(c("intercept",dimnames(w)[[2]]))),
  mean = means,
  upper = upper,
  lower = lower,
  method = factor(rep(c("VB","MCMC"), each=length(means)/2), levels=c("MCMC","VB"))
)
random = data.frame(
  MCMC = apply(extract(ps.spline,"b")$b, 2, mean), 
  VB = vb_result$mub.q
)

# ELBO
plot(vb_result$lb, xlab="Iteration", ylab="", type="l", main="ELBO")

# Fixed effects
ggplot(fixed, aes(x=names, y=mean, color=method)) + 
  geom_errorbar(data=fixed, mapping=aes(ymin=lower, ymax=upper), position="dodge") +
  geom_point(position=position_dodge(width=0.9), size=1.6) +
  geom_hline(yintercept = 0, color="black", size=0.7, linetype="dotted") +
  scale_color_manual(values=c("deepskyblue4","darkorange3")) +
  scale_x_discrete(
    breaks=c("male","female","age>=50","asian","intercept"),
    labels=c("male","female","age>=50","asian","intercept")) +
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
ggplot(random, aes(x=MCMC, y=VB)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  theme_classic() +
  ggtitle("Random effects")

# Mean curve
fxGrid.spline = extract(ps.spline, 'fxGrid')$fxGrid
fhat.spline.m = colMeans(fxGrid.spline)
fhat.spline.u = apply(fxGrid.spline, 2, function(x) quantile(x, 0.975))
fhat.spline.l = apply(fxGrid.spline, 2, function(x) quantile(x, 0.025))
betaHat =apply(extract(ps.spline,'beta')$beta,2,mean)
wbm = W%*%betaHat
r = logb2_GM - wbm - b[studycode]
plot(r ~ logUcd_GM, xlab=expression(log(UCd[GM])), ylab='Residual', main="Mean Curve(Measurement Error)")
# polygon(c(xGrid,rev(xGrid)),c(fhat.spline.l,rev(fhat.spline.u)),lty=2,col='gray')
lines(fhat.spline.m~xGrid,lwd=2,col=1)
# points(r~logUcd_GM,lwd=2)
points(x,y-cbind(1,w)%*%vb_result$mubeta.q-Z%*%vb_result$mub.q,col=2)
ord = order(vb_result$xgrid)
lines(vb_result$xgrid[ord],vb_result$post_curve[ord],lwd=2,col=2)
legend("topleft",c("MCMC","VB"),col=c("black","red"),lty=c(1,1))
