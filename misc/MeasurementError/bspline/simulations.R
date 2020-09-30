library(rstan)
library(splines)
library(lattice)
library(ggplot2)

# 1. Standard Nonparametric regression ------------------------------------
# Replicating example provided in Appendix of Wand(2008)
x = environmental$radiation
y = (environmental$ozone)^(1/3)
source("regression.R")
n_intknot = 2

# VB
vb_result = regression(y, x, n_intknot=n_intknot)

# MCMC
resolution = 200 # number of x's to estimate f(x)
order = 4 # To build cubic B-spline basis
knots = quantile(unique(x), seq(0,1,length=n_intknot+2)[-c(1,n_intknot+2)])
boundary = c(min(x)-sd(x)/2, max(x)+sd(x)/2)
xgrid = seq(boundary[1], boundary[2], length.out=resolution) 
BsGrid = bs(x=xgrid, knots=knots, intercept=TRUE, Boundary.knots=boundary)
data = list(
    num_data=length(y),
    n_intknot=n_intknot,
    knots=knots,
    boundary=boundary,
    order=order,
    resolution=resolution,
    Y=y, X=x,
    BsGrid=BsGrid)
# set.seed(100)
# mcmc_result = stan("./stan_repo/nonparametric/code.stan", data=data, iter=1e3)
# write.csv(extract(mcmc_result, 'fxGrid')$fxGrid, "./stan_repo/nonparametric/samples_nonparametric_fxGrid.csv", row.names=FALSE)
# write.csv(extract(mcmc_result, 'beta')$beta, "./stan_repo/nonparametric/samples_nonparametric_beta.csv", row.names=FALSE)
mcmc_curve = as.matrix(read.csv("./stan_repo/nonparametric/samples_nonparametric_fxGrid.csv", header=TRUE))
mcmc_beta = as.matrix(read.csv("./stan_repo/nonparametric/samples_nonparametric_beta.csv", header=TRUE))

# Result presentation
polygon_x = c(xgrid, rev(xgrid))
vb_y = vb_result$mubeta.q + c(vb_result$post_lower, rev(vb_result$post_upper))
mcmc_y = mean(mcmc_beta) + c(apply(mcmc_curve,2,quantile,0.025), rev(apply(mcmc_curve,2,quantile,0.975)))

par(mfrow=c(1,2), mar=c(4,4,4,4))
plot(x,y,xlab="radiation",ylab="cube root of ozone",main="Variational Bayes",pch=19, cex=0.5)
polygon(polygon_x, vb_y, col="darkgoldenrod1", lty = "blank")
points(x,y,pch=19,cex=0.5)
lines(xgrid, vb_result$mubeta.q + vb_result$post_curve, col="darkgoldenrod4", lwd=4)

plot(x,y,xlab="radiation",ylab="cube root of ozone",main="MCMC",pch=19, cex=0.5)
polygon(polygon_x, mcmc_y, col="deepskyblue1", lty = "blank")
points(x,y,pch=19,cex=0.5)
lines(xgrid, mean(mcmc_beta) + colMeans(mcmc_curve), col="deepskyblue4", lwd=4)

# # Result presentation(ggplot version)
# # Intercept
# vb_intercept = vb_result$mubeta.q
# mcmc_intercept = mean(mcmc_beta)
# # Polygon
# vb_polygon = data.frame(
#   variate = c(vb_result$xgrid, rev(vb_result$xgrid)),
#   area = c(vb_intercept + vb_result$post_upper, vb_intercept + rev(vb_result$post_lower)))
# mcmc_polygon = data.frame(
#   variate = c(xgrid, rev(xgrid)),
#   area = c(mcmc_intercept + apply(mcmc_curve,2,quantile,0.025), mcmc_intercept + rev(apply(mcmc_curve,2,quantile,0.975))))
# # Mean curves
# ozone = data.frame(x=x, y=y)
# mc = ggplot() +
#   geom_point(data=ozone, mapping=aes(x=x,y=y)) +
#   geom_polygon(data=mcmc_polygon,aes(x=variate,y=area), fill="deepskyblue4", alpha=0.3) +
#   geom_polygon(data=vb_polygon,aes(x=variate,y=area), fill="darkorange3", alpha=0.3) +
#   geom_line(data=data.frame(x=xgrid, y=mcmc_intercept+colMeans(mcmc_curve)),
#             mapping=aes(x=x,y=y),color="deepskyblue4",size=1.5) +
#   geom_line(data=data.frame(x=vb_result$xgrid, y=vb_intercept+vb_result$post_curve),
#             mapping=aes(x=x,y=y),color="darkorange3",size=1.5) +
#   xlab("radiation") + ylab("cuberoot of ozone") + theme_bw() + ggtitle("Estimated Mean curve")

# 2. Semiparametric regression with Measurement Error -----------------

source("semimer.R")
n_intknot = 5

# Data generation
set.seed(10)
N = 130
D = 6
RR = 0.9
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 1.5
beta = rnorm(D+1)
w = matrix(rnorm(N*D),N,D)
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
f = function(x) 2*x+2*sin(pi*x)
y = drop(cbind(1,w)%*%beta + f(x) + rnorm(N, sd=0.7))

# VB
vb_result = semimer(y,w,v,n_intknot=5)

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
    int_knots = knots,
    boundary = boundary,
    order = order,
    resolution = resolution,
    y = y,
    v = v,
    bsgrid = bsgrid,
    w = w)
# Commented out process that took 706.588 seconds
# set.seed(100)
# mcmc_result = stan("./stan_repo/me_semireg/code.stan", data=data, iter=1e3)
# write.csv(extract(mcmc_result, 'fxGrid')$fxGrid, "./stan_repo/me_semireg/samples_me_semireg_fxGrid.csv", row.names=FALSE)
# write.csv(extract(mcmc_result, 'beta')$beta, "./stan_repo/me_semireg/samples_me_semireg_beta.csv", row.names=FALSE)
mcmc_curve = as.matrix(read.csv("./stan_repo/me_semireg/samples_me_semireg_fxGrid.csv", header=TRUE))
mcmc_beta = as.matrix(read.csv("./stan_repo/me_semireg/samples_me_semireg_beta.csv", header=TRUE))

# Result presentation
par(mfrow=c(1,2), mar=c(4,4,4,4))
plot(vb_result$lb, main="ELBO plot", ylab="ELBO", xlab="Iterations")
plot(colMeans(mcmc_beta)[-1], vb_result$mubeta.q[-1], main="Fixed effects\n(Comparison)", xlab="MCMC", ylab="VB")
lines(-10:10,-10:10)

vb_x = c(vb_result$xgrid, rev(vb_result$xgrid))
vb_y = vb_result$mubeta.q[1] + c(vb_result$post_lower, rev(vb_result$post_upper))
mcmc_x = c(xgrid, rev(xgrid))
mcmc_y = colMeans(mcmc_beta)[1] + c(apply(mcmc_curve,2,quantile,0.025), rev(apply(mcmc_curve,2,quantile,0.975)))

plot(x,y-w%*%vb_result$mubeta.q[-1], main="Variational Bayes", ylab="Residual", cex=0.5, pch=19)
polygon(vb_x, vb_y, col="darkgoldenrod1", lty = "blank")
points(x,y-w%*%vb_result$mubeta.q[-1], cex=0.5, pch=19)
lines(vb_result$xgrid, vb_result$post_curve+vb_result$mubeta.q[1], col="darkgoldenrod4", lwd=4)

plot(x,y-w%*%colMeans(mcmc_beta)[-1], main="MCMC", ylab="Residual", cex=0.5, pch=19)
polygon(mcmc_x, mcmc_y, col="deepskyblue1", lty="blank")
points(x,y-w%*%colMeans(mcmc_beta)[-1],cex=0.5, pch=19)
lines(xgrid, colMeans(mcmc_curve)+colMeans(mcmc_beta)[1], col="deepskyblue4", lwd=4)

# 3. Random Intercept model with Measurement Error ------------------------

source("randintmer.R")
n_intknot = 5

# Data generation
source("../../misc/make_Z.R") # Builds design matrix for random effects
set.seed(10)
N = 80
T = 6
D = 5
RR = 0.9
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 1.5
beta = rnorm(D+1)
w = matrix(rnorm(N*T*D), ncol=D)
x = rnorm(N*T, mux, sqrt(xi2))
v = rnorm(N*T, x, sqrt(sig2v))
Z = make_Z(rep(T,N))
b = rnorm(N)
f = function(x) 2*x+cos(pi*x)
y = drop(cbind(1,w)%*%beta + Z%*%b + f(x) + rnorm(N, sd=0.5))

# VB (6.763 seconds)
vb_result = randintmer(y,w,v,Z,n_intknot=n_intknot)

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
# Commented out process that took 9688.34 seconds
# set.seed(100)
# mcmc_result = stan("./stan_repo/me_randint/code.stan", data=data, iter=1e3)
# write.csv(extract(mcmc_result, 'fxGrid')$fxGrid, "./stan_repo/me_randint/samples_me_randint_fxGrid.csv", row.names=FALSE)
# write.csv(extract(mcmc_result, 'beta')$beta, "./stan_repo/me_randint/samples_me_randint_beta.csv", row.names=FALSE)
# write.csv(extract(mcmc_result, 'b')$b, "./stan_repo/me_randint/samples_me_randint_b.csv", row.names=FALSE)
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
