randint_dpm = function(y, w, v, Z, R=10, n_intknot=30, prior=NULL, maxiter=500, tol=1e-3, n_grids=1e3, resolution=200)
{
  N = ncol(Z)
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  ZtZ = crossprod(Z)
  ti = apply(Z, 2, sum)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2mu = sig2beta = 100
    aalp = balp = alam = blam = anu = aphi = bphi = aups = bups = asig = bsig = 1e-3
    bnu = 1e-4
  }
  else
  {
    aalp = prior$aalp
    balp = prior$balp
    sig2mu = prior$sig2mu
    sig2beta = prior$sig2beta
    alam = prior$alam
    blam = prior$blam
    anu = prior$anu
    bnu = prior$bnu
    aphi = prior$aphi
    bphi = prior$bphi
    aups = prior$aups
    bups = prior$bups
    asig = prior$asig
    bsig = prior$bsig
  }
  
  # Compose interior knots of the basis
  intknots = seq(min(v), max(v), length.out=n_intknot+2)[-c(1,n_intknot+2)] # equally spaced knots
  # intknots = quantile(range(v), seq(0,1,length.out=n_grids+2)[-c(1,n_grids+2)]) # quantile based knots
  
  # Compose spline basis of grids
  grids = seq(min(v)-sd(v)/2, max(v)+sd(v)/2, length.out=n_grids)
  vphig = matrix(grids, n_grids, n_intknot) - outer(rep(1,n_grids), intknots)
  vphig = vphig * (vphig > 0)
  
  # Initialize variational parameters
  aalptl = aalp + R-1
  anutl = anu + sum(ti)/2
  aphitl = aphi + N/2
  aupstl = aups + n_intknot/2
  asigtl = asig + sum(ti)/2
  alp.ratio = aalp/balp
  nu.ratio = anu/bnu
  phi.ratio = aphi/bphi
  ups.ratio = aups/bups
  sig.ratio = asig/bsig
  muu.q = rep(0, n_intknot)
  sigu.q = diag(rep(ups.ratio, n_intknot))
  vphiq = matrix(0, sum(ti), n_intknot)
  mub.q = rep(0, N)
  gam1s = gam2s = rep(1, R-1)
  logstick = digamma(gam1s) - digamma(gam1s+gam2s)
  log1mstick = digamma(gam2s) - digamma(gam1s+gam2s)
  alamtls = blamtls = rep(0.01, R)
  lam.ratio = alamtls/blamtls
  kappa = matrix(1, nrow=sum(ti), ncol=R)/R
  mutls = rep(mean(v), R)
  sigtls = rep(var(v), R)
  revcumsum = function(x) {rev(cumsum(rev(x)))}
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # beta
    sigbeta.q = solve(diag(rep(1/sig2beta,D+1))+sig.ratio*WtW)
    mubeta.q = drop(sig.ratio*sigbeta.q%*%t(W)%*%(y-Z%*%mub.q-vphiq%*%muu.q))
    
    # random intercepts
    sigb.q = solve(diag(rep(phi.ratio,N))+sig.ratio*ZtZ)
    mub.q = drop(sig.ratio*sigb.q%*%t(Z)%*%(y-W%*%mubeta.q-vphiq%*%muu.q))
    
    # denoised value
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(muu.q,muu.q)+sigu.q)%*%t(vphig)) + nu.ratio*(grids^2))
    lnpgrids = common + nu.ratio*outer(grids, v) + sig.ratio*outer(drop(vphig%*%muu.q), drop(y-W%*%mubeta.q-Z%*%mub.q))
    lnpgrids = lnpgrids + outer(grids, apply(t(kappa)*lam.ratio*mutls,2,sum)) - 0.5*outer(grids^2, apply(t(kappa)*lam.ratio,2,sum))
    pgrids = exp(t(lnpgrids))
    normalizers = apply(pgrids, 1, sum)
    ex = drop(pgrids%*%grids)/normalizers
    ex2 = drop(pgrids%*%(grids^2))/normalizers
    varx = ex2 - (ex)^2
    vphiq = pgrids%*%vphig/outer(normalizers,rep(1,n_intknot))
    vphiqtvphiq = t(vphig)%*%diag(apply(pgrids/normalizers,2,sum))%*%vphig
    
    # kappa
    lnkappa = (outer(rep(1,R),ex)-mutls)^2 + sigtls # R x sum(ti)
    lnkappa = -0.5*lam.ratio*t(t(lnkappa) + varx) - 0.5*(log(blamtls) - digamma(alamtls)) + c(logstick,0) + c(0,cumsum(log1mstick)) # R x sum(ti)
    kappa = exp(t(lnkappa) - matrixStats::colLogSumExps(lnkappa)) # sum(ti) x R
    
    # stick length
    gam1s = 1 + apply(kappa,2,sum)[-R]
    revcs_kappa = kappa
    for (ridx in 1:nrow(kappa)) revcs_kappa[ridx,] = revcumsum(revcs_kappa[ridx,])
    gam2s = alp.ratio + apply(revcs_kappa,2,sum)[-1]
    logstick = digamma(gam1s) - digamma(gam1s+gam2s)
    log1mstick = digamma(gam2s) - digamma(gam1s+gam2s)
    
    # alpha
    balptl = balp - sum(log1mstick)
    alp.ratio = aalptl/balptl
    
    # parameters(mean)
    sigtls = 1/(1/sig2mu + lam.ratio*apply(kappa,2,sum))
    mutls = sigtls*lam.ratio*apply(kappa*ex,2,sum)
    
    # parameters(variance)
    alamtls = alam + 0.5*apply(kappa,2,sum)
    blamtls = blam + 0.5*apply(kappa*(t((outer(rep(1,R),ex)-mutls)^2+sigtls)+varx),2,sum)
    lam.ratio = alamtls/blamtls
    
    # nu
    bnutl = bnu + 0.5*(sum((v-ex)^2)+sum(varx))
    nu.ratio = anutl/bnutl
    
    # phi
    bphitl = bphi + 0.5*sum((mub.q)^2+diag(sigb.q))
    phi.ratio = aphitl/bphitl
    
    # spline coefficients
    sigu.q = solve(diag(rep(ups.ratio,n_intknot))+sig.ratio*vphiqtvphiq)
    muu.q = drop(sig.ratio*sigu.q%*%t(vphiq)%*%(y-W%*%mubeta.q-Z%*%mub.q))
    
    # upsilon
    bupstl = bups + 0.5*sum((muu.q)^2+diag(sigu.q))
    ups.ratio = aupstl/bupstl
    
    # sigma
    cpterm = sum((y-W%*%mubeta.q-Z%*%mub.q)^2) + sum(diag(WtW%*%sigbeta.q)) + sum(diag(ZtZ%*%sigb.q))
    cpterm = cpterm - 2*sum((y-W%*%mubeta.q-Z%*%mub.q)*(vphiq%*%muu.q)) + sum(diag(vphiqtvphiq%*%(sigu.q+outer(muu.q,muu.q))))
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # ELBO
    lbnew = -0.5*sum(ti)*log(2*pi) - 0.5*sum(ti)*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*sum(ti)*log(2*pi) - 0.5*sum(ti)*(log(bnutl)-digamma(anutl)) - 0.5*nu.ratio*(sum((v-ex)^2)+sum(varx))
    lbnew = lbnew + sum(kappa*(-0.5*outer(rep(1,sum(ti)),log(blamtls)-digamma(alamtls)) - 0.5*t(((outer(rep(1,R),ex)-mutls)^2)/lam.ratio)))
    lbnew = lbnew - sum(kappa*(-0.5*outer(varx,rep(1,R))-0.5))
    lbnew = lbnew - 0.5*(D+1)*log(sig2beta) - 0.5*sum(mubeta.q^2+diag(sigbeta.q))/sig2beta
    lbnew = lbnew + 0.5*determinant(sigbeta.q,logarithm=TRUE)$modulus[1] + 0.5*(D+1)
    lbnew = lbnew - 0.5*N*(log(bphitl)-digamma(aphitl)) - 0.5*phi.ratio*sum(mub.q^2+diag(sigb.q))
    lbnew = lbnew + 0.5*determinant(sigb.q,logarithm=TRUE)$modulus[1] + 0.5*N
    lbnew = lbnew - 0.5*n_intknot*(log(bupstl)-digamma(aupstl)) - 0.5*ups.ratio*sum(muu.q^2+diag(sigu.q))
    lbnew = lbnew + 0.5*determinant(sigu.q,logarithm=TRUE)$modulus[1] + 0.5*n_intknot
    lbnew = lbnew + sum(-0.5*log(sig2mu)-0.5*(mutls+sigtls)/sig2mu)
    lbnew = lbnew - sum(-0.5*log(sigtls)-0.5)
    lbnew = lbnew + sum(-lgamma(alam) + alam*log(blam) - (alam+1)*(log(blamtls)-digamma(alamtls)) - lam.ratio*blam)
    lbnew = lbnew - sum(-lgamma(alamtls) + alamtls*log(blamtls) - (alamtls+1)*(log(blamtls)-digamma(alamtls)) - lam.ratio*blamtls)
    lbnew = lbnew + sum(-log(balptl)+digamma(aalptl)+(alp.ratio-1)*log1mstick)
    lbnew = lbnew - sum(lgamma(gam1s+gam2s)-lgamma(gam1s)-lgamma(gam2s)+(gam1s-1)*logstick+(gam2s-1)*log1mstick)
    lbnew = lbnew + sum(t(t(kappa)*(c(logstick,0) + c(0,cumsum(log1mstick)))) - kappa*log(kappa))
    lbnew = lbnew - lgamma(aalp) + aalp*log(balp) + (aalp-1)*(-log(balptl)+digamma(aalptl)) - alp.ratio*balp
    lbnew = lbnew + lgamma(aalptl) - aalptl*log(balptl) - (aalptl-1)*(-log(balptl)+digamma(aalptl)) + alp.ratio*balptl
    lbnew = lbnew - lgamma(anu) + anu*log(bnu) - (anu+1)*(log(bnutl)-digamma(anutl)) - nu.ratio*bnu
    lbnew = lbnew + lgamma(anutl) - anutl*log(bnutl) + (anutl+1)*(log(bnutl)-digamma(anutl)) + nu.ratio*bnutl
    lbnew = lbnew - lgamma(aphi) + aphi*log(bphi) - (aphi+1)*(log(bphitl)-digamma(aphitl)) - phi.ratio*bphi
    lbnew = lbnew + lgamma(aphitl) - aphitl*log(bphitl) + (aphitl+1)*(log(bphitl)-digamma(aphitl)) + phi.ratio*bphitl
    lbnew = lbnew - lgamma(asig) + asig*log(bsig) - (asig+1)*(log(bsigtl)-digamma(asigtl)) - sig.ratio*bsig
    lbnew = lbnew + lgamma(asigtl) - asigtl*log(bsigtl) + (asigtl+1)*(log(bsigtl)-digamma(asigtl)) + sig.ratio*bsigtl
    lbnew = lbnew - lgamma(aups) + aups*log(bups) - (aups+1)*(log(bupstl)-digamma(aupstl)) - ups.ratio*bups
    lbnew = lbnew + lgamma(aupstl) - aupstl*log(bupstl) + (aupstl+1)*(log(bupstl)-digamma(aupstl)) + ups.ratio*bupstl
    lb[iter] = lbnew
    if (abs(lbnew-lbold)<tol) break
    lbold = lbnew
  }
  lb = lb[1:iter]
  xgrid = seq(min(v)-sd(v)/2, max(v)+sd(v)/2, length.out=resolution) 
  vphi = matrix(xgrid, resolution, n_intknot) - outer(rep(1,resolution),intknots)
  vphi = vphi * (vphi > 0)
  post_curve=drop(vphi%*%muu.q)
  
  return(list(
    lb=lb, ex=ex, varx=varx, mubeta.q=mubeta.q, sigbeta.q=sigbeta.q, mub.q=mub.q, sigb.q=sigb.q, muu.q=muu.q, sigu.q=sigu.q,
    sig.ratio=sig.ratio, nu.ratio=nu.ratio, phi.ratio=phi.ratio, ups.ratio=ups.ratio, mutls=mutls, sigtls=sigtls,
    xgrid=xgrid,
    post_curve=post_curve,
    post_lower=vphi%*%qnorm(0.025,muu.q,sqrt(diag(sigu.q))),
    post_upper=vphi%*%qnorm(0.975,muu.q,sqrt(diag(sigu.q)))
  ))
}