vbsf_exp = function(y,w,Z,productivity=TRUE,prior=NULL,tol=1e-5,maxiter=500){
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
    rsig.0 = rlam.0 = 1.1
    ssig.0 = slam.0 = 1.1
    sigbeta2 = 100}
  else{
    rsig.0 = prior$rsig.0
    rlam.0 = prior$rlam.0
    ssig.0 = prior$ssig.0
    slam.0 = prior$slam.0
    sigbeta2 = prior$sigbeta2}
  # parametric
  sb2diag = diag(1/sigbeta2,D+1)
  rsig.q = rsig.0+N
  rlam.q = rlam.0+2*n
  sig.ratio = rsig.0/ssig.0
  lam.ratio = rlam.0/slam.0
  mubeta.q = rep(0,D+1)
  muu.q = rep(0,n)
  # Iteration
  for(iter in 1:maxiter){
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    # beta
    sigbeta.q = sig.ratio*WtW+sb2diag
    sigbeta.q = solve(sigbeta.q)
    mubeta.q = sig.ratio*sigbeta.q%*%t(W)%*%(y-sgn*Z%*%muu.q)
    mubeta.q = drop(mubeta.q)
    # u
    sig2 = 1/(sig.ratio*ti)
    sig = sqrt(sig2)
    mu = -lam.ratio+sgn*sig.ratio*(t(Z)%*%(y-W%*%mubeta.q))
    mu = sig2*drop(mu)
    musigrat = mu/sig
    fres = exp(dnorm(musigrat,log=TRUE) - pnorm(musigrat,log.p=TRUE))
    muu.q = mu+sig*fres
    sigu.q = diag(sig2*(1-fres*(musigrat+fres)))
    # lambda
    slam.q = slam.0+2*sum(muu.q)
    lam.ratio = rlam.q/slam.q
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(diag(ZtZ%*%sigu.q))
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
    rlam.q=rlam.q,slam.q=slam.q))}
# # Simulation
# source("../misc/wandint.R")
# source("../misc/make_Z.R")
# set.seed(2)
# n = 50
# ti = 4; D = 10
# Z = make_Z(rep(ti, n))
# sigma = 0.5
# lambda = 2
# 
# betaT = rnorm(D+1)
# w = matrix(rnorm(nrow(Z)*D), ncol=D)
# u = rexp(n,lambda)
# y = drop(cbind(1,w)%*%betaT-rep(u,each=ti)+rnorm(nrow(Z),sd=sigma))
# result = vbsf_exp(y,w,Z)
# 
# plot(result$muu.q, u)
# lines(-5:5, -5:5)
# plot(result$mubeta.q, betaT)
# lines(-5:5, -5:5)
# sqrt(result$ssig.q/result$rsig.q)