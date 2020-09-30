
# Data handling -----------------------------------------------------------

library(dplyr)
cadmium = read.csv("./cadmium.csv", header=TRUE)
study = cadmium$studycode   # study code
w = cbind(ifelse(cadmium$ethnicity == 1, 1, 0),
          ifelse(cadmium$age >= 50, 1, 0),
          ifelse(cadmium$gender == 1, 1, 0),
          ifelse(cadmium$gender == 0, 1, 0))
colnames(w)=c('asian','age>=50','female','male')
y = log(cadmium$b2_GM)
x = log(cadmium$Ucd_GM)

# sort the data w.r.t. studycode, to make correspondence with MCMC
cadmium_bind = cbind.data.frame(study,w,y,x) %>% arrange(study)
study = cadmium_bind$study
y = cadmium_bind$y
x = cadmium_bind$x; ord = order(x)
w = cadmium_bind %>% dplyr::select(-study,-y,-x) %>% as.matrix()
label = unique(study)
Z = matrix(0, nrow(w), length(label))
for(idx in 1:length(label)) Z[which(study==label[idx]),idx] = 1


# Random Intercept MER ----------------------------------------------------

source("./randintmer.R")

# VB
vb_result = randintmer(y,w,x,Z,n_intknot=2)

# MCMC
resolution = 200 # number of x's to estimate f(x)
order = 4 # To build cubic B-spline basis
knots = quantile(unique(v), seq(0,1,length=n_intknot+2)[-c(1,n_intknot+2)])

boundary = c(min(v)-sd(v)/2, max(v)+sd(v)/2)
xgrid = seq(boundary[1], boundary[2], length.out=resolution) 
bsgrid = bs(x=xgrid, knots=knots, intercept=TRUE, Boundary.knots=boundary)
data = list(
  n_obs = length(y),
  n_intknot = length(knots),
  n_regcoef = D,
  n_objects = N,
  int_knots = knots,
  boundary = boundary,
  order = order,
  resolution = resolution,
  y = y,
  v = v,
  bsgrid = bsgrid,
  w = w,
  Z = Z)
set.seed(100)
mcmc_result = stan("./stan_repo/me_randint/code.stan", data=data, iter=1e3)
write.csv(extract(mcmc_result, 'fxGrid')$fxGrid, "./stan_repo/me_randint/samples_cadmium_fxGrid.csv", row.names=FALSE)
write.csv(extract(mcmc_result, 'beta')$beta, "./stan_repo/me_randint/samples_cadmium_beta.csv", row.names=FALSE)
write.csv(extract(mcmc_result, 'b')$b, "./stan_repo/me_randint/samples_cadmium_b.csv", row.names=FALSE)
mcmc_curve = as.matrix(read.csv("./stan_repo/me_randint/samples_me_randint_fxGrid.csv", header=TRUE))
mcmc_beta = as.matrix(read.csv("./stan_repo/me_randint/samples_me_randint_beta.csv", header=TRUE))
mcmc_b = as.matrix(read.csv("./stan_repo/me_randint/samples_me_randint_b.csv", header=TRUE))

plot(vb_result$lb, main="ELBO plot", ylab="ELBO", xlab="Iterations", sub=paste("Converged at",round(vb_result$lb[length(vb_result$lb)],3)))

par(mfrow=c(1,2), mar=c(4,4,4,4))
plot(colMeans(mcmc_beta)[-1], vb_result$mubeta.q[-1], main="Fixed effects\n(Comparison)", xlab="MCMC", ylab="VB")
lines(-10:10,-10:10)
plot(colMeans(mcmc_b), vb_result$mub.q, main="Random effects\n(Comparison)", xlab="MCMC", ylab="VB")
lines(-10:10,-10:10)

vb_x = c(vb_result$xgrid, rev(vb_result$xgrid))
vb_y = vb_result$mubeta.q[1] + c(vb_result$post_lower, rev(vb_result$post_upper))
mcmc_x = c(xgrid, rev(xgrid))
mcmc_y = colMeans(mcmc_beta)[1] + c(apply(mcmc_curve,2,quantile,0.025), rev(apply(mcmc_curve,2,quantile,0.975)))

plot(x,y-w%*%vb_result$mubeta.q[-1]-Z%*%vb_result$mub.q, main="Variational Bayes", ylab="Residual", cex=0.5, pch=19)
polygon(vb_x, vb_y, col="darkgoldenrod1", lty = "blank")
points(x,y-w%*%vb_result$mubeta.q[-1]-Z%*%vb_result$mub.q, cex=0.5, pch=19)
lines(vb_result$xgrid, vb_result$post_curve+vb_result$mubeta.q[1], col="darkgoldenrod4", lwd=4)

plot(x,y-w%*%colMeans(mcmc_beta)[-1]-Z%*%vb_result$mub.q, main="MCMC", ylab="Residual", cex=0.5, pch=19)
polygon(mcmc_x, mcmc_y, col="deepskyblue1", lty="blank")
points(x,y-w%*%colMeans(mcmc_beta)[-1]-Z%*%vb_result$mub.q,cex=0.5, pch=19)
lines(xgrid, colMeans(mcmc_curve)+colMeans(mcmc_beta)[1], col="deepskyblue4", lwd=4)
