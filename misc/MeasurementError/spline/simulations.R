summon = dir()[!dir() == "simulations.R"]
for (filename in summon) source(paste0("./", filename))


# 1. Simple Linear Regression ---------------------------------------------

set.seed(10)
N = 130
RR = 0.8
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 3
beta = c(0.5, 1.1)
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
y = cbind(1,x)%*%beta + rnorm(N)
result = simplereg(y,v)
par(mfrow=c(1,2))
plot(result$lb,xlab="Iteration",ylab="",main="Evidence Lower Bound", type="l")
plot(x,result$ex,xlab="True",ylab="Estimate",main="Corrupted variable")
lines(-10:10,-10:10)
result$mubeta.q
par(mfrow=c(1,1))


# 2. Nonparametric Regression ---------------------------------------------

set.seed(10)
N = 130
RR = 0.7
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 1.5
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
f = function(x) 2*x+sin(pi*x)
y = f(x) + rnorm(N)
vb_result = nonparam(y,v)
ord = order(vb_result$ex)
plot(x,y,ylab="y",main="Original pattern(dot) vs Denoised pattern(line)")
lines(vb_result$ex[ord],vb_result$post_curve[ord],lwd=3,col=2)
lines(vb_result$ex[ord],vb_result$post_lower[ord],lwd=2,col=3)
lines(vb_result$ex[ord],vb_result$post_upper[ord],lwd=2,col=3)


# 3. Semiparametric Regression --------------------------------------------

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
f = function(x) 2*x+sin(pi*x)
y = cbind(1,w)%*%beta + f(x) + rnorm(N)
vb_result = semiparam(y,w,v)
plot(beta[-1],vb_result$mubeta.q[-1],main="Regression coefficiencts")
lines(-10:10,-10:10)
ord = order(vb_result$ex)
plot(x, y-cbind(1,w)%*%vb_result$mubeta.q, ylab="residual",main="Denoised pattern")
lines(vb_result$ex[ord],vb_result$post_curve[ord],lwd=3,col=2)
lines(vb_result$ex[ord],vb_result$post_lower[ord],lwd=2,col=3)
lines(vb_result$ex[ord],vb_result$post_upper[ord],lwd=2,col=3)


# 4. Semiparametric Random Intercept Model --------------------------------

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
u = rnorm(N)
f = function(x) 2*x+sin(pi*x)
y = cbind(1,w)%*%beta + Z%*%u + f(x) + rnorm(N, sd=0.5)
vb_result = randint(y,w,v,Z)
plot(beta[-1],vb_result$mubeta.q[-1],xlab="Original",ylab="Estimated",main="Regression coefficiencts")
lines(-10:10,-10:10)
plot(u,vb_result$muu.q,xlab="Original",ylab="Estimated",main="Random Intercepts")
lines(-10:10,-10:10)
ord = order(vb_result$ex)
plot(x, y-cbind(1,w)%*%vb_result$mubeta.q-Z%*%vb_result$muu.q, ylab="residual",main="Denoised pattern")
lines(vb_result$ex[ord],vb_result$post_curve[ord],lwd=3,col=2)
lines(vb_result$ex[ord],vb_result$post_lower[ord],lwd=2,col=3)
lines(vb_result$ex[ord],vb_result$post_upper[ord],lwd=2,col=3)
