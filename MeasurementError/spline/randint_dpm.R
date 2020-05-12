randint_dpm = function(y, w, v, Z, R=10, n_intknot=30, prior=NULL, maxiter=500, tol=1e-3, n_grids=1e3, resolution=200)
{
  
  N = ncol(Z)
  D = ncol(w)
  W = cbind(1,w)
  WtW = crossprod(W)
  ZtZ = crossprod(Z)
  ti = apply(Z, 2, sum)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2mu = sig2beta = 100
    aalp = balp = alam = blam = asig = bsig = anu = aphi = bphi = aups = bups = 1e-3
    bnu = 1e-4
  }
  else
  {
    aalp = prior$aalp
    balp = prior$balp
    sig2mu = prior$sig2mu
    alam = prior$alam
    blam = prior$blam
    sig2beta = prior$sig2beta
    asig = prior$asig
    bsig = prior$bsig
    anu = prior$anu
    bnu = prior$bnu
    aphi = prior$aphi
    bphi = prior$bphi
    aups = prior$aups
    bups = prior$bups
  }
  
  # Compose interior knots of the basis
  intknots = seq(min(v), max(v), length.out=n_intknot+2)[-c(1,n_intknot+2)] # equally spaced knots
  # intknots = quantile(range(v), seq(0,1,length.out=n_grids+2)[-c(1,n_grids+2)]) # quantile based knots
  
  # Compose spline basis of grids
  grids = seq(min(v)-sd(v)/2, max(v)+sd(v)/2, length.out=n_grids)
  vphig = matrix(grids, n_grids, n_intknot) - outer(rep(1,n_grids), intknots)
  vphig = vphig * (vphig > 0)
  
  # Initialize the variational parameters
  # Level 1
  mubeta.q = rep(0, D+1)
  mub.q = rep(0, N)
  muu.q = rep(0, n_intknot)
  sigu.q = diag(rep(aups/bups, n_intknot))
  asigtl = asig + (sum(ti)/2)
  sig.ratio = asig/bsig
  vphiq = matrix(0, nrow=sum(ti), ncol=n_intknot)
  # Level 2
  kappa = matrix(1,sum(ti),R) / R
  lam.ratio = rep(alam/blam, R)
  mutls = rep(mean(v), R)
  sigtls = rep(var(v), R)
  # Level 3
  revcumsum = function(x) {rev(cumsum(rev(x)))}
  aalptl = aalp + (R-1)
  alp.ratio = aalp/balp
  gam1s = gam2s = rep(1e-2, R-1)
  logstick = digamma(gam1s) - digamma(gam1s+gam2s)
  log1mstick = digamma(gam2s) - digamma(gam1s+gam2s)
  qstick = c(logstick, 0) + c(0, cumsum(log1mstick))
  anutl = anu + (sum(ti)/2)
  nu.ratio = anu/bnu
  aphitl = aphi + (N/2)
  phi.ratio = aphi/bphi
  aupstl = aups + (n_intknot/2)
  ups.ratio = aups/bups
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    mub.q.old = mub.q
    muu.q.old = muu.q
    qstick.old = qstick
    
    # beta
    sigbeta.q = solve(diag(rep(sig.ratio, D+1)) + sig.ratio*WtW)
    mubeta.q = drop(sig.ratio*sigbeta.q%*%t(W)%*%(y-Z%*%mub.q-vphiq%*%muu.q))
    
    # random intercepts
    sigb.q = solve(diag(rep(phi.ratio, N)) + sig.ratio*ZtZ)
    mub.q = drop(sig.ratio*sigb.q%*%t(Z)%*%(y-W%*%mubeta.q-vphiq%*%muu.q))
    
    # denoised values
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(muu.q,muu.q)+sigu.q)%*%t(vphig)) + nu.ratio*grids^2)
    lnpgrids = common + nu.ratio*outer(grids, v) + sig.ratio*outer(drop(vphig%*%muu.q), drop(y-W%*%mubeta.q-Z%*%mub.q))
    lnpgrids = lnpgrids + outer(grids, apply(t(kappa)*lam.ratio*mutls,2,sum)) - 0.5*outer(grids^2, apply(t(kappa)*lam.ratio,2,sum))
    pgrids = exp(t(lnpgrids))
    normalizers = apply(pgrids, 1, sum)
    ex = drop(pgrids%*%grids)/normalizers
    ex2 = drop(pgrids%*%(grids^2))/normalizers
    varx = ex2 - (ex)^2
    vphiq = pgrids%*%vphig/outer(normalizers,rep(1,n_intknot))
    vphiqtvphiq = t(vphig)%*%diag(apply(pgrids/normalizers,2,sum))%*%vphig
    
    # spline coefficients
    sigu.q = solve(diag(rep(ups.ratio,n_intknot))+sig.ratio*vphiqtvphiq)
    muu.q = drop(sig.ratio*sigu.q%*%t(vphiq)%*%(y-W%*%mubeta.q-Z%*%mub.q))
    
    # parameters(mean)
    sigtls = 1/(1/sig2mu + apply(t(kappa)*lam.ratio,1,sum))
    mutls = sigtls*apply(t(kappa*ex)*lam.ratio,1,sum)
    
    # parameters(variance)
    alamtl = alam + 0.5*apply(kappa,2,sum)
    blamtl = blam + 0.5*apply(kappa*(t((outer(rep(1,R),ex)-mutls)^2+sigtls)+varx),2,sum)
    lam.ratio = alamtl/blamtl
    
    # kappa
    S = t(lam.ratio*(t((outer(rep(1,R),ex)-mutls)^2+sigtls)+varx)) - digamma(alamtl) + log(blamtl)
    S = outer(rep(1,sum(ti)),qstick) - 0.5*t(S)
    kappa = exp(S - matrixStats::rowLogSumExps(S))
    
    # stick length
    gam1s = 1 + apply(kappa,2,sum)[-R]
    revcs_kappa = kappa
    for (ridx in 1:nrow(kappa)) revcs_kappa[ridx,] = revcumsum(revcs_kappa[ridx,])
    gam2s = alp.ratio + apply(revcs_kappa,2,sum)[-1]
    logstick = digamma(gam1s) - digamma(gam1s+gam2s)
    log1mstick = digamma(gam2s) - digamma(gam1s+gam2s)
    qstick = c(logstick, 0) + c(0, cumsum(log1mstick))
    
    # alpha
    balptl = balp - sum(log1mstick)
    alp.ratio = aalptl/balptl
    
    # sigma
    cpterm = sum((y-W%*%mubeta.q-Z%*%mub.q)^2) + sum(diag(WtW%*%sigbeta.q)) + sum(diag(ZtZ%*%sigb.q))
    cpterm = cpterm - 2*sum((y-W%*%mubeta.q-Z%*%mub.q)*(vphiq%*%muu.q)) + sum(diag(vphiqtvphiq%*%(sigu.q+outer(muu.q,muu.q))))
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # nu
    bnutl = bnu + 0.5*sum((v-ex)^2+varx)
    nu.ratio = anutl/bnutl
    
    # phi
    bphitl = bphi + 0.5*sum(mub.q^2+diag(sigb.q))
    phi.ratio = aphitl/bphitl
    
    # upsilon
    bupstl = bups + 0.5*sum(muu.q^2+diag(sigu.q))
    ups.ratio = aupstl/bupstl
    
    # ELBO
    
  }
  browser()
}