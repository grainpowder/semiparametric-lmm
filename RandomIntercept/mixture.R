normaldpm_lmm = function(y,w,Z,R,prior=NULL,maxiter=500,tol=1e-5)
{
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  N = ncol(Z)
  T = nrow(Z)
  ti = apply(Z, 2, sum)
  # prior parameters
  # Note : be sure to make base distribution as flat as possible
  if (is.null(prior))
  {
    asig = bsig = 0.1
    aalp = balp = 0.1
    alam = blam = 0.01
    mu0 = 0
    sig02 = 100
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(100, D+1)
  }
  else
  {
    asig = prior$asig
    bsig = prior$bsig
    aalp = prior$aalp
    balp = prior$balp
    alam = prior$alam
    blam = prior$blam
    mu0 = prior$mu0
    sig02 = prior$sig02
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
  }
  # Define containers
  # Level1
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  sigu.q = rep(sig02, N)
  sig.ratio = asig / bsig
  # Level2
  kappa = matrix(runif(N*R), N, R)
  kappa = kappa / apply(kappa, 1, sum)
  mutls = rep(0, R)
  sigtls = rep(sig02, R)
  alamtls = rep(alam, R)
  blamtls = rep(blam, R)
  lam.ratio = alamtls / blamtls
  # Level3 or above
  gam1s = gam2s = rep(0.1, R-1)
  alp.ratio = aalp / balp
  # Define predetermined values
  sb0i = solve(sigbeta.0)
  sb0imb0 = sb0i %*% mubeta.0
  asigtl = asig + (N*T/2)
  aalptl = aalp + R - 1
  # Function to use during iteration
  numreplace = function(x, tol) {x[x<tol] = tol; x} # function to replace value(s) less than 1e-10 to 1e-10
  revcumsum = function(x) {rev(cumsum(rev(x)))}
  # Begin iteration
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = sigbeta.q %*% (sb0imb0+sig.ratio*t(W)%*%(y-Z%*%muu.q))
    mubeta.q = drop(mubeta.q)
    
    # sigma
    ssterm = sum((y - W%*%mubeta.q - Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW %*% sigbeta.q))
    trterm2 = sum(ti * sigu.q)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = asigtl / bsigtl
    
    # random intercept
    kappa_lam = t(kappa) * lam.ratio
    sigu.q = 1 / (apply(kappa_lam,2,sum)+ti*sig.ratio)
    muu.q = apply(kappa_lam*mutls,2,sum) + sig.ratio*t(Z)%*%(y-W%*%mubeta.q)
    muu.q = drop(sigu.q * muu.q)
    
    # kappa(prior...)
    logstick = digamma(gam1s) - digamma(gam1s+gam2s)
    log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
    qstick = c(logstick, 0) + c(0, log1mstick)
    
    Sir = (t(matrix(muu.q, N, R))-mutls)^2 + sigtls
    Sir = t(t(Sir)+sigu.q)*lam.ratio - digamma(alamtls) + log(blamtls)
    Sir = t(qstick - 0.5*Sir)
    Sir = exp(Sir - apply(Sir,1,max))
    for (cidx in 1:ncol(Sir)) Sir[, cidx] = numreplace(Sir[, cidx],tol=1e-10)
    kappa = Sir / apply(Sir,1,sum)
    
    # stick length
    gam1s = 1 + apply(kappa,2,sum)[-R]
    revcs_kappa = kappa
    for (ridx in 1:nrow(kappa)) revcs_kappa[ridx,] = revcumsum(revcs_kappa[ridx,])
    gam2s = alp.ratio + apply(revcs_kappa[, -1],2,sum)
    
    # alpha
    balptl = balp - sum(digamma(gam2s)-digamma(gam1s+gam2s))
    alp.ratio = aalptl / balptl
    
    # parameters: mean
    sigtls = 1 / (1/sig02+apply(t(kappa)*lam.ratio, 1, sum))
    mutls = sigtls * (mu0/sig02+apply(t(kappa*muu.q)*lam.ratio, 1, sum))
    
    # parameters: variance
    alamtls = alam + 0.5*apply(kappa,2,sum)
    ssterm_mu = (t(matrix(muu.q, N, R))-mutls)^2 + sigtls
    ssterm_mu = t(ssterm_mu) + sigu.q
    ssterm_mu = kappa*ssterm_mu
    blamtls = blam + 0.5*apply(ssterm_mu,2,sum)
    lam.ratio = alamtls / blamtls
    
    # Convergence
    convergence1 = sum((mubeta.q - mubeta.q.old)^2) < tol
    convergence2 = sum((muu.q - muu.q.old)^2) < tol
    if (convergence1 & convergence2) break
  }
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    mutls=mutls, sigtls=sigtls,
    asigtl=asigtl, bsigtl=bsigtl,
    kappa=kappa
  ))
}