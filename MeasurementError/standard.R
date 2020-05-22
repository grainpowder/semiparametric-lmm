randomintercept = function(y, w, x, Z, n_intknot=30, prior=NULL, maxiter=500, tol=1e-4)
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
    sig2beta = 100
    aphi = bphi = aups = bups = asig = bsig = 1e-3
  }
  else
  {
    sig2beta = prior$sig2beta
    aphi = prior$aphi
    bphi = prior$bphi
    aups = prior$aups
    bups = prior$bups
    asig = prior$asig
    bsig = prior$bsig
  }
  
  # Compose interior knots of the basis
  intknots = seq(min(x), max(x), length.out=n_intknot+2)[-c(1,n_intknot+2)] # equally spaced knots
  # intknots = quantile(range(x), seq(0,1,length.out=n_grids+2)[-c(1,n_grids+2)]) # quantile based knots
  vphi = matrix(x, length(x), n_intknot) - outer(rep(1,length(x)), intknots)
  vphi = vphi * (vphi > 0)
  vphitvphi = crossprod(vphi)
  
  # Initialize variational parameters
  asigtl = asig + (sum(ti)/2)
  aphitl = aphi + (N/2)
  aupstl = aups + (n_intknot/2)
  sig.ratio = asig/bsig
  phi.ratio = aphi/bphi
  ups.ratio = aups/bups
  muu.q = rep(0, n_intknot)
  mub.q = rep(0, N)
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # beta
    sigbeta.q = solve(diag(1/sig2beta,D+1) + sig.ratio*WtW)
    mubeta.q = drop(sig.ratio*sigbeta.q%*%t(W)%*%(y-Z%*%mub.q-vphi%*%muu.q))
    
    # random intercept
    sigb.q = solve(diag(phi.ratio,N) + sig.ratio*ZtZ)
    mub.q = drop(sig.ratio*sigb.q%*%t(Z)%*%(y-W%*%mubeta.q-vphi%*%muu.q))
    
    # phi
    bphitl = bphi + 0.5*sum(mub.q^2+diag(sigb.q))
    phi.ratio = aphitl/bphitl
    
    # spline coefficients
    sigu.q = solve(diag(ups.ratio,n_intknot)+sig.ratio*vphitvphi)
    muu.q = drop(sig.ratio*sigu.q%*%t(vphi)%*%(y-W%*%mubeta.q-Z%*%mub.q))
    
    # upsilon
    bupstl = bups + 0.5*sum(muu.q^2+diag(sigu.q))
    ups.ratio = aupstl/bupstl
    
    # sigma
    cpterm = sum((y-W%*%mubeta.q-Z%*%mub.q-vphi%*%muu.q)^2)
    cpterm = cpterm + sum(diag(WtW%*%sigbeta.q)) + sum(diag(ZtZ%*%sigb.q)) + sum(diag(vphitvphi*sigu.q))
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # ELBO
    lbnew = -0.5*sum(ti)*log(2*pi) - 0.5*sum(ti)*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*(D+1)*log(sig2beta) - 0.5*sum(mubeta.q^2+diag(sigbeta.q))/sig2beta
    lbnew = lbnew + 0.5*determinant(sigbeta.q,logarithm=TRUE)$modulus[1] + 0.5*(D+1)
    lbnew = lbnew - 0.5*N*(log(bphitl)-digamma(aphitl)) - 0.5*phi.ratio*sum(mub.q^2+diag(sigb.q))
    lbnew = lbnew + 0.5*determinant(sigb.q,logarithm=TRUE)$modulus[1] + 0.5*N
    lbnew = lbnew - 0.5*n_intknot*(log(bupstl)-digamma(aupstl)) - 0.5*ups.ratio*sum(muu.q^2+diag(sigu.q))
    lbnew = lbnew + 0.5*determinant(sigu.q,logarithm=TRUE)$modulus[1] + 0.5*n_intknot
    lbnew = lbnew - lgamma(aphi) + aphi*log(bphi) - (aphi+1)*(log(bphitl)-digamma(aphitl)) - phi.ratio*bphi
    lbnew = lbnew + lgamma(aphitl) - aphitl*log(bphitl) + (aphitl+1)*(log(bphitl)-digamma(aphitl)) + phi.ratio*bphitl
    lbnew = lbnew - lgamma(aups) + aups*log(bups) - (aups+1)*(log(bupstl)-digamma(aupstl)) - ups.ratio*bups
    lbnew = lbnew + lgamma(aupstl) - aupstl*log(bupstl) + (aupstl+1)*(log(bupstl)-digamma(aupstl)) + ups.ratio*bupstl
    lbnew = lbnew - lgamma(asig) + asig*log(bsig) - (asig+1)*(log(bsigtl)-digamma(asigtl)) - sig.ratio*bsig
    lbnew = lbnew + lgamma(asigtl) - asigtl*log(bsigtl) + (asigtl+1)*(log(bsigtl)-digamma(asigtl)) + sig.ratio*bsigtl
    lb[iter] = lbnew
    if (abs(lbnew-lbold)<tol) break
    lbold = lbnew
  }
  lb = lb[1:iter]
  post_curve=drop(vphi%*%muu.q)
  
  return(list(
    lb=lb, mubeta.q=mubeta.q, sigbeta.q=sigbeta.q, mub.q=mub.q, sigb.q=sigb.q, muu.q=muu.q, sigu.q=sigu.q,
    sig.ratio=sig.ratio, phi.ratio=phi.ratio, ups.ratio=ups.ratio,
    post_curve=post_curve,
    post_lower=vphi%*%qnorm(0.025,muu.q,sqrt(diag(sigu.q))),
    post_upper=vphi%*%qnorm(0.975,muu.q,sqrt(diag(sigu.q)))
  ))
}