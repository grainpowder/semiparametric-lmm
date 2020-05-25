regression = function(y, x, w=NULL, n_intknot=10, prior=NULL, maxiter=1000, tol=1e-4, resolution=200)
{
  library(splines)
  N = length(y)
  D = ncol(w); if (is.null(D)) D = 0
  W = cbind(matrix(1,nrow=N), w)
  WtW = crossprod(W)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2beta = 100
    aups = bups = asig = bsig = 1e-3
  }
  else
  {
    sig2beta = prior$sig2beta
    aups = prior$aups
    bups = prior$bups
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = diag(rep(1/sig2beta,D+1))
  if (D == 0) sb0i = 1/sig2beta
  
  # Build spline basis
  intKnots = quantile(unique(x), seq(0,1,length=n_intknot+2)[-c(1,n_intknot+2)])
  boundary = c(min(x)-sd(x)/2, max(x)+sd(x)/2)
  vphi = bs(x, knots=intKnots, Boundary.knots=boundary, intercept=TRUE)
  vphitvphi = crossprod(vphi)
  n_knot = ncol(vphi)
  
  # Initialize variational parameters
  aupstl = aups + (n_knot/2)
  asigtl = asig + (N/2)
  mubeta.q = solve(WtW)%*%t(W)%*%y
  ups.ratio = aups/bups
  sig.ratio = asig/bsig
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # Spline coefficients
    sigu.q = solve(diag(rep(ups.ratio,n_knot))+sig.ratio*vphitvphi)
    muu.q = drop(sig.ratio*sigu.q%*%t(vphi)%*%(y-W%*%mubeta.q))
    
    # upsilon
    bupstl = bups + 0.5*sum(muu.q^2+diag(sigu.q))
    ups.ratio = aupstl/bupstl
    
    # beta
    sigbeta.q = solve(sb0i+sig.ratio*WtW)
    mubeta.q = drop(sig.ratio*sigbeta.q%*%t(W)%*%(y-vphi%*%muu.q))
    
    # sigma
    cpterm = sum((y-W%*%mubeta.q-vphi%*%muu.q)^2)
    cpterm = cpterm + sum(diag(WtW%*%sigbeta.q)) + sum(diag(vphitvphi%*%sigu.q))
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # ELBO
    lbnew = -0.5*N*log(2*pi) - 0.5*N*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*(D+1)*log(sig2beta) - 0.5*sum(mubeta.q^2+diag(sigbeta.q))/sig2beta
    lbnew = lbnew + 0.5*determinant(sigbeta.q,logarithm=TRUE)$modulus[1] + 0.5*(D+1)
    lbnew = lbnew - 0.5*n_knot*(log(bupstl)-digamma(aupstl)) - 0.5*ups.ratio*sum(muu.q^2+diag(sigu.q))
    lbnew = lbnew + 0.5*determinant(sigu.q,logarithm=TRUE)$modulus[1] + 0.5*n_knot
    lbnew = lbnew - lgamma(asig) + asig*log(bsig) - (asig+1)*(log(bsigtl)-digamma(asigtl)) - sig.ratio*bsig
    lbnew = lbnew + lgamma(asigtl) - asigtl*log(bsigtl) + (asigtl+1)*(log(bsigtl)-digamma(asigtl)) + sig.ratio*bsigtl
    lbnew = lbnew - lgamma(aups) + aups*log(bups) - (aups+1)*(log(bupstl)-digamma(aupstl)) - ups.ratio*bups
    lbnew = lbnew + lgamma(aupstl) - aupstl*log(bupstl) + (aupstl+1)*(log(bupstl)-digamma(aupstl)) + ups.ratio*bupstl
    lb[iter] = lbnew
    if (abs(lbnew-lbold)<tol) break
    lbold = lbnew
  }
  lb = lb[1:iter]
  post_curve=drop(vphi%*%muu.q)
  return(list(
    lb=lb, mubeta.q=mubeta.q, sigbeta.q=sigbeta.q, muu.q=muu.q, sigu.q=sigu.q,
    sig.ratio=sig.ratio, ups.ratio=ups.ratio,
    post_curve=post_curve,
    post_lower=vphi%*%qnorm(0.025,muu.q,sqrt(diag(sigu.q))),
    post_upper=vphi%*%qnorm(0.975,muu.q,sqrt(diag(sigu.q)))
  ))
}