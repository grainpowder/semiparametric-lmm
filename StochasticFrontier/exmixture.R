exmixture = function(y, w, Z, R=10, productivity=TRUE, prior=NULL, maxiter=500, eps=1e-4)
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
    alam = blam = aalp = balp = asig = bsig = 1e-3
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    alam = prior$alam
    blam = prior$blam
    aalp = prior$aalp
    balp = prior$balp
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  aalptl = aalp + R - 1
  asigtl = asig + 0.5*(sum(ti))
  
  # Initialize parameters
  # Level1
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  sig.ratio = asig / bsig
  # Level2
  kappa = matrix(1,N,R) / R
  kappa_r = apply(kappa,2,sum)
  lam.ratio = rep(alam/blam, R)
  # Level3 or above
  gam1s = rep(0.01, R-1)
  gam2s = rep(0.01, R-1)
  alp.ratio = aalp / balp
  logstick = digamma(gam1s) - digamma(gam1s+gam2s)
  log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
  qstick = c(logstick, 0) + c(0, log1mstick)
  
  # functions to be used
  times = function(x, y) x*y
  revcumsum = function(x) {rev(cumsum(rev(x)))}
  
  # Update as
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
    # beta
    sigbeta.q = solve(sb0i+sig.ratio*WtW)
    mubeta.q = sigbeta.q%*%(sbimb0+sig.ratio*t(W)%*%(y-sgn*Z%*%muu.q))
    mubeta.q = drop(mubeta.q)
    
    # inefficiency
    sigi2 = 1/(sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2 * (-apply(t(kappa)*lam.ratio,2,sum) + sgn*sig.ratio*drop(t(Z)%*%(y-W%*%mubeta.q)))
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE)-pnorm(musigrat, log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2 * (1-fres*(musigrat+fres))
    
    # lambda_r
    alamtl = alam + kappa_r
    blamtl = blam + apply(kappa*muu.q, 2, sum)
    lam.ratio = alamtl / blamtl
    
    # kappa
    S = outer(muu.q, lam.ratio, times)
    S = t(qstick - log(blamtl) + digamma(alamtl) + t(S))
    kappa = exp(kappa - rowLogSumExps(kappa))
    kappa_r = apply(kappa, 2, sum)
    
    # stick length(Vr)
    gam1s = 1 + apply(kappa,2,sum)[-R]
    revcs_kappa = kappa
    for (ridx in 1:nrow(kappa)) revcs_kappa[ridx,] = revcumsum(revcs_kappa[ridx,])
    gam2s = alp.ratio + apply(revcs_kappa[, -1],2,sum)
    logstick = digamma(gam1s) - digamma(gam1s+gam2s)
    log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
    qstick = c(logstick, 0) + c(0, log1mstick)
    
    # alpha
    balptl = balp - sum(digamma(gam2s)-digamma(gam1s+gam2s))
    alp.ratio = aalptl / balptl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ti*sigu.q)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = asigtl / bsigtl
    
    # Convergence
    bool1 = mean((mubeta.q.old - mubeta.q)^2) < eps
    bool2 = mean((muu.q.old - muu.q)^2) < eps
    if (bool1 & bool2) break
  }
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl,
    alamtl=alamtl, blamtl=blamtl
  ))
}