source("regression.R")
set.seed(10)
N = 100
f = function(x) 2*x+sin(2*pi*x)
x = runif(N) * 2

# Nonparametric -----------------------------------------------------------

set.seed(10)
y = 2 + f(x) + rnorm(N, sd=0.5)
vb_result = regression(y,x)
resid = y - vb_result$mubeta.q
plot(x,resid,main="Nonparametric")
ord = order(x)
lines(x[ord],vb_result$post_curve[ord],lwd=2,col=2)
lines(x[ord],vb_result$post_upper[ord],lwd=2,lty=2)
lines(x[ord],vb_result$post_lower[ord],lwd=2,lty=2)
plot(vb_result$lb, ylab="ELBO", xlab="Iteration", main="ELBO plot")

# Semiparametric ----------------------------------------------------------

set.seed(10)
D = 5
w = matrix(rnorm(N*D), ncol=D)
beta = rnorm(D+1)
y = cbind(1,w)%*%beta + f(x) + rnorm(N, sd=0.5)
vb_result = regression(y,x,w)
resid = y - cbind(1,w)%*%vb_result$mubeta.q
plot(x,resid,main="Semiparametric")
ord = order(x)
lines(x[ord],vb_result$post_curve[ord],lwd=2,col=2)
lines(x[ord],vb_result$post_upper[ord],lwd=2,lty=2)
lines(x[ord],vb_result$post_lower[ord],lwd=2,lty=2)
plot(beta[-1], vb_result$mubeta.q[-1], xlab="True", ylab="Estimated", main="Fixed Effects")
lines(-10:10, -10:10)
plot(vb_result$lb, ylab="ELBO", xlab="Iteration", main="ELBO plot")