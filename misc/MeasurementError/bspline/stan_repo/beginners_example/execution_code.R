library(rstan)


# 1. Estimation on true mean ----------------------------------------------

N = 130
y = rnorm(N, mean=1.6, sd=0.7)
data = list(N=N, y=y)
result = stan("./stan_repo/true_mean.stan", data = data, iter = 1e3)
summary(result@sim$samples[[1]]$mu)
summary(result@sim$samples[[1]]$sigma)


# 2. Linear Regression ----------------------------------------------------

N = 834
D = 10
w = matrix(rnorm(N*D), ncol=D)
beta = rnorm(D+1)
y = drop(cbind(1,w)%*%beta + rnorm(N))
data = list(N=N, D=D, y=y, w=w)
result = stan("./stan_repo/basic_regression.stan", data = data, iter = 1e3)
plot(beta, colMeans(extract(result, 'beta')$beta), xlab="True", ylab="Stan", main="Coefficients")
lines(-10:10,-10:10)
summary(extract(result, 'sigma')$sigma)
