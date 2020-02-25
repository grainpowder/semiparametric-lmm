vbsf_tn = function(y, w, Z, productivity=TRUE, prior=NULL, tol=1e-5, maxiter=500)
{
  sgn = (-1)^productivity
  N = ncol(Z)
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  ti = apply(Z, 2, sum)
  
  # Hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1))
    mu0 = 0
    sig02 = 100
    alam = 1.1
    blam = 0.8
    asig = 0.001
    bsig = 0.001
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    mu0 = prior$mu0
    sig02 = prior$sig02
    alam = prior$alam
    blam = prior$blam
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  mucoef = -0.5 * (1/sig02+N)
  musigrat0 = mu0 / sig02
  llamcoef = alam + (N/2) - 1
  asigtl = asig + nrow(Z)/2
  
  # Initialize parameters
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  mu = 0
  lambda2 = lambda1 = 1
  sig.ratio = asig / bsig
  
  # Functions to be used in wandint
  const = function(x) x-x+1
  mom2 = function(x) x^2
  mom1 = function(x) x
  mom0.5 = function(x) sqrt(x)
  
  # Update as
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = sigbeta.q %*% (sbimb0 + sig.ratio*t(W)%*%(y-sgn*Z%*%muu.q))
    mubeta.q = drop(mubeta.q)
    
    # inefficiency
    sigi2 = 1 / (lambda2+sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2 * (mu*lambda1+sgn*sig.ratio*t(Z)%*%(y-W%*%mubeta.q))
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE) - pnorm(musigrat, log.p = TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2 * (1-fres*(musigrat+fres))
    
    # mu
    integrand = function(m) mucoef*(m^2) + (musigrat0+lambda1*sum(muu.q))*m - N*pnorm(m, log.p=TRUE)
    lnC1 = wandint(const, integrand)
    mu = exp(wandint(mom1, integrand) - lnC1)
    # lambda
    integrand = function(l) llamcoef*log(l) - (blam+sum(muu.q^2+sigu.q)/2)*l + mu*sum(muu.q)*sqrt(l)
    lnC2 = wandint(const, integrand)
    lambda2 = exp(wandint(mom1, integrand) - lnC2)
    lambda1 = exp(wandint(mom0.5, integrand) - lnC2)
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ti*sigu.q)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = asigtl / bsigtl
    
    # Convergence
    bool1 = mean((mubeta.q.old - mubeta.q)^2) < tol
    bool2 = mean((muu.q.old - muu.q)^2) < tol
    if (bool1 & bool2) break
  }
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl
  ))
}


# Simulation
library(truncnorm)
source("../misc/make_Z.R")
source("../misc/wandint.R")
D = 10
N = 50
ti = rep(4, N)
Z = make_Z(ti)
w = matrix(rnorm(D*nrow(Z)), ncol=D)
u = rtruncnorm(N, a=0, mean=1)
beta = rnorm(D+1)
sgn = (-1)^TRUE
y = cbind(1,w)%*%beta + sgn*Z%*%u + rnorm(nrow(Z))
result = vbsf_tn(y,w,Z)
par(mfrow=c(1,2),mar=c(4,3.7,2,1))
plot(beta,result$mubeta.q)
lines(-10:10,-10:10)
plot(u,result$muu.q)
lines(-10:10,-10:10)
