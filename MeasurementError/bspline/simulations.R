

# 1. Standard Non/Semiparametric ------------------------------------------

source("regression.R")
set.seed(10)
N = 100
f = function(x) 2*x+sin(2*pi*x)
x = runif(N) * 2

# Nonparametric

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

# Semiparametric

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


# 2. Semiparametric regression with Measurement Error -----------------

source("semimer.R")
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
y = cbind(1,w)%*%beta + f(x) + rnorm(N, sd=0.7)
vb_result = semimer(y,w,v,n_intknot=5)
plot(vb_result$lb, main="ELBO plot", ylab="ELBO", xlab="Iterations", sub=paste("Converged at",round(vb_result$lb[length(vb_result$lb)],3)))
plot(beta[-1], vb_result$mubeta.q[-1], main="Fixed effects", xlab="True", ylab="Estimated")
lines(-10:10,-10:10)
plot(x,y-cbind(1,w)%*%vb_result$mubeta.q, main="Estimated mean curve", ylab="Residual")
ord = order(vb_result$xgrid)
lines(vb_result$xgrid[ord],vb_result$post_curve[ord],lwd=3,col=2)
lines(vb_result$xgrid[ord],vb_result$post_lower[ord],lwd=2,lty=2)
lines(vb_result$xgrid[ord],vb_result$post_upper[ord],lwd=2,lty=2)


# 3. Random Intercept model with Measurement Error ------------------------

source("randintmer.R")
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
y = cbind(1,w)%*%beta + Z%*%b + f(x) + rnorm(N, sd=0.5)
vb_result = randintmer(y,w,v,Z,n_intknot=5)
plot(vb_result$lb, main="ELBO plot", ylab="ELBO", xlab="Iterations", sub=paste("Converged at",round(vb_result$lb[length(vb_result$lb)],3)))
plot(beta[-1],vb_result$mubeta.q[-1],xlab="True",ylab="Estimated",main="Regression coefficiencts")
lines(-10:10,-10:10)
plot(b,vb_result$mub.q,xlab="True",ylab="Estimated",main="Random Intercepts")
lines(-10:10,-10:10)
ord = order(vb_result$xgrid)
plot(x, y-cbind(1,w)%*%vb_result$mubeta.q-Z%*%vb_result$mub.q, ylab="Residual",main="Estimated mean curve")
lines(vb_result$xgrid[ord],vb_result$post_curve[ord],lwd=3,col=2)
lines(vb_result$xgrid[ord],vb_result$post_lower[ord],lwd=2,col=3)
lines(vb_result$xgrid[ord],vb_result$post_upper[ord],lwd=2,col=3)
