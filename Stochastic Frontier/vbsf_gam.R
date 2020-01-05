vbsf_gam = function(y,w,Z,productivity=TRUE,prior=NULL,tol=1e-4,maxiter=500){
  # Pre-calculate objects to be used frequently
  sgn = (-1)^productivity
  n = ncol(Z)
  N = nrow(w)
  D = ncol(w)
  W = cbind(1,w)
  WtW = crossprod(W)
  ZtZ = crossprod(Z)
  ti = apply(Z,2,sum)
  # Initialize parameters
  if(is.null(prior)){
    rthe.0 = rsig.0 = 0.1
    slam.0 = sthe.0 = ssig.0 = 0.1
    sigbeta2 = 100}
  else{
    rthe.0 = prior$rthe.0
    rsig.0 = prior$rsig.0
    slam.0 = prior$slam.0
    sthe.0 = prior$sthe.0
    ssig.0 = prior$ssig.0
    sigbeta2 = prior$sigbeta2}
  # parametric
  low=1e-2
  sb2diag = diag(1/sigbeta2,D+1)
  rsig.q = rsig.0+N
  sig.ratio = lam.ratio = the.ratio = 1
  mubeta.q = solve(t(W)%*%W)%*%t(W)%*%y
  muu.q = sigu.q = lnu = rep(low,n)
  # Predefine f functions to be frequently used
  const = function(x) x-x+1
  x1 = function(x) x
  x2 = function(x) x^2
  lx = function(x) log(x)
  # Iteration
  for (iter in 1:maxiter){
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    # beta
    sigbeta.q = sig.ratio*WtW+sb2diag
    sigbeta.q = (sigbeta.q+t(sigbeta.q))/2
    sigbeta.q = solve(sigbeta.q)
    mubeta.q = sig.ratio*sigbeta.q%*%t(W)%*%(y-sgn*Z%*%muu.q)
    mubeta.q = drop(mubeta.q)
    # u(can be improved by using foreach functionality)
    res = drop(t(Z)%*%(y-W%*%mubeta.q))
    for (i in 1:n) {
      integrand = function(u) (the.ratio-1)*log(u)+(-lam.ratio+sgn*res[i])*u-(sig.ratio*ti[i]/2)*(u^2)
      lnC       = wandint(const,  integrand, a=low, b=7, k=10.1, init=muu.q[i]+1, log.value=TRUE)
      first     = exp(wandint(x1, integrand, a=low, b=7, k=10.1, init=muu.q[i]+1, log.value=TRUE)-lnC)
      second    = exp(wandint(x2, integrand, a=low, b=7, k=10.1, init=muu.q[i]+1, log.value=TRUE)-lnC)
      lnu[i]    = wandint(lx, integrand, a=low, b=7, k=10.1, init=muu.q[i]+1, log.value=FALSE)/exp(lnC)
      muu.q[i]  = first
      sigu.q[i] = second-first^2}
    # lambda
    rlam.q = (n+1)*the.ratio
    slam.q = slam.0+sum(muu.q)
    lam.ratio = rlam.q/slam.q
    llamb = -log(slam.q)+digamma(rlam.q)
    # theta
    integrand = function(th) -(n+1)*log(gamma(th)) + (sum(lnu)+(n+1)*llamb+log(slam.0)-sthe.0)*th + (rthe.0-1)*log(th)
    lnC       = wandint(const,  integrand, a=low, b=7, k=10.1, init=the.ratio, log.value=TRUE)
    the.ratio = exp(wandint(x1, integrand, a=low, b=7, k=10.1, init=the.ratio, log.value=TRUE)-lnC)
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^ 2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ZtZ%*%sigu.q)
    ssig.q = ssig.0+ssterm+trterm1+trterm2
    sig.ratio = rsig.q/ssig.q
    # Convergence
    mse1 = sqrt(mean((mubeta.q-mubeta.q.old)^2))
    mse2 = sqrt(mean((muu.q-muu.q.old)^2))
    if(mse1<tol && mse2<tol) break}
  return(list(
    mubeta.q=mubeta.q,sigbeta.q=sigbeta.q,
    muu.q=muu.q,sigu.q=sigu.q,
    rsig.q=rsig.q,ssig.q=ssig.q,
    Elam=lam.ratio,Etheta=the.ratio))}

source("../misc/wandint.R")
source("../misc/make_Z.R")
set.seed(1)
n = 50
ti = 4; D = 10
Z = make_Z(rep(ti, n))
sigma = 0.5
theta = 2
lambda = 2

betaT = rnorm(D+1)
w = matrix(rnorm(nrow(Z)*D), ncol=D)
u = rgamma(n,shape=theta,rate=lambda)
y = drop(cbind(1,w)%*%betaT-rep(u,each=ti)+rnorm(nrow(Z),sd=sigma))
result = vbsf_gam(y,w,Z)

plot(result$muu.q, u)
sqrt(result$ssig.q/result$rsig.q)
result$Elam
result$Etheta
