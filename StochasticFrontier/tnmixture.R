tnmixture = function(y, w, Z, R=10, productivity=TRUE, prior=NULL, maxiter=500, eps=1e-4, acc=1000)
{
  library(matrixStats)
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
    alam = blam = aalp = balp = asig = bsig = 1e-3
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    mu0 = prior$mu0
    sig02 = prior$sig02
    alam = prior$alam
    blam = prior$blam
    aalp = prior$aalp
    balp = prior$balp
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  ldetsb0 = determinant(sigbeta.0, logarithm=TRUE)
  aalptl = aalp + R - 1
  asigtl = asig + 0.5*(sum(ti))
  
  # Initialize parameters
  # Level1
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  sig.ratio = asig / bsig
  # Level2
  b_lam2 = rep(1, R)
  b_lam1 = rep(1, R)
  b_mu1 = rep(0, R)
  kappa = matrix(1,N,R) / R
  # Level3 or above
  gam1s = rep(1, R-1)
  gam2s = rep(1, R-1)
  alp.ratio = aalp / balp
  
  # grids to be used in griddy gibbs routine
  real = seq(-8, 8, length.out=acc) # Interval wider than this would cause numerical overflow
  if (sum(real == 0) == 1) real = real[-which(real == 0)]
  neg = real < 0; pos = real > 0
  positive = seq(eps, 20, length.out=acc)
  if (sum(positive == 1) == 1) positive = positive[-which(positive == 1)]
  less = positive < 1; greater = positive > 1
  
  # functions to be used
  times = function(x, y) x*y
  revcumsum = function(x) {rev(cumsum(rev(x)))}
  
  # Iterate
  for (iter in 1:maxiter)
  {
    # beta test(see if algorithm converges to right posterior)
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = sigbeta.q %*% (sbimb0 + sig.ratio*t(W)%*%(y-sgn*Z%*%muu.q))
    mubeta.q = drop(mubeta.q)
    
    # inefficiency
    sigi2 = 1 / (apply(t(kappa)*b_lam2,2,sum) + sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2 * (apply(t(kappa)*b_lam1*b_mu1,2,sum) + sgn*sig.ratio*drop(t(Z)%*%(y-W%*%mubeta.q)))
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE)-pnorm(musigrat, log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2 * (1-fres*(musigrat+fres))
    
    # griddy Gibbs routine
    kappa_r = apply(kappa,2,sum)
    # mu
    # Hg construction
    Hg = outer(real, mu0/sig02+b_lam1*apply(kappa*muu.q,2,sum), times)
    Hg = Hg - outer(real^2, 0.5*(1/sig02+kappa_r), times)
    Hg = Hg - outer(pnorm(real, log.p=TRUE), kappa_r, times)
    # Moment calculation
    xi_mu = colLogSumExps(Hg)
    b_mu1 = exp(colLogSumExps(Hg[pos,]+log(real[pos])) - xi_mu) - exp(colLogSumExps(Hg[neg,]+log(-real[neg])) - xi_mu)
    # b_mu2 = exp(colLogSumExps(Hg+log(real^2)) - xi_mu)
    b_mu2 = b_mu1^2
    b_lpmu = -exp(colLogSumExps(Hg+log(-log(pnorm(real)))) - xi_mu)
    # ilambda2
    Hg = outer(log(positive), alam+0.5*kappa_r-1, times)
    Hg = Hg - outer(positive, blam+0.5*apply(kappa*(muu.q^2+sigu.q),2,sum), times)
    Hg = Hg + outer(sqrt(positive), b_mu1*apply(kappa*muu.q,2,sum), times)
    # Moment calcuation
    xi_lam2 = colLogSumExps(Hg)
    b_lam2 = exp(colLogSumExps(Hg+log(positive)) - xi_lam2)
    # b_lam1 = exp(colLogSumExps(Hg+0.5*log(positive)) - xi_lam2)
    b_lam1 = sqrt(b_lam2)
    b_llam2 = exp(colLogSumExps(Hg[greater,]+log(log(positive[greater]))) - xi_lam2) - exp(colLogSumExps(Hg[less,]+log(-log(positive[less]))) - xi_lam2)
    
    # kappa
    logstick = digamma(gam1s) - digamma(gam1s+gam2s)
    log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
    qstick = c(logstick, 0) + c(0, log1mstick)
    S = -0.5*outer(muu.q^2+sigu.q, b_lam2, times)
    S = S + outer(muu.q, b_lam1*b_mu1, times)
    S = t(t(S) + 0.5*(b_llam2-b_mu2) - b_lpmu + qstick)
    kappa = exp(S - rowLogSumExps(S))
    
    # stick length
    gam1s = 1 + apply(kappa,2,sum)[-R]
    revcs_kappa = kappa
    for (ridx in 1:nrow(kappa)) revcs_kappa[ridx,] = revcumsum(revcs_kappa[ridx,])
    gam2s = alp.ratio + apply(revcs_kappa[, -1],2,sum)
    
    # alpha
    balptl = balp - sum(digamma(gam2s)-digamma(gam1s+gam2s))
    alp.ratio = aalptl / balptl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ti*sigu.q)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = asigtl / bsigtl
    
    # # ELBO
    # term1 = -0.5*nrow(W) * (log(2*pi)-digamma(bsigtl)+log(asigtl))
    # term1 = term1 - 0.5*sig.ratio*(ssterm+trterm1+trterm2)
    # term2 = 0.5*()
    
    # Convergence
    bool1 = mean((mubeta.q.old - mubeta.q)^2) < eps
    bool2 = mean((muu.q.old - muu.q)^2) < eps
    if (bool1 & bool2) break
  }
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl,
    mur=b_mu1, lam2r=b_lam2
  ))
}

source("../misc/make_Z.R")
library(truncnorm)
set.seed(10)
N = 50
Z = make_Z(rep(4, N))
D = 10
w = matrix(rnorm(D*nrow(Z)), ncol=D)
# u = rtruncnorm(N, a=0, mean=-1)
u = runif(N) * 3
beta = rnorm(D+1)
y = cbind(1, w)%*%beta + ((-1)^TRUE)*Z%*%u + rnorm(nrow(Z))
result = tnmixture(y,w,Z)

plot(beta,result$mubeta.q)
lines(-10:10,-10:10)
plot(u,result$muu.q)
lines(-10:10,-10:10)
upper = qtruncnorm(0.975,a=0,mean=result$muu.q, sd=sqrt(result$sigu.q))
lower = qtruncnorm(0.025,a=0,mean=result$muu.q, sd=sqrt(result$sigu.q))
plot(1:N, u,xlab="idx",ylab="",pch=19,ylim=c(0,max(upper)),col=2,main="True values and corresponding 95% Credible interval\n(griddy Gibbs)")
for (idx in 1:N) lines(c(idx,idx), c(upper[idx],lower[idx]))
points(1:N, result$muu.q, pch=19)

