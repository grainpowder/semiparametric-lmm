nonparam = function(y, v, n_intknot=30, prior=NULL, maxiter=500, tol=1e-4, n_grids=1e3, resolution=200)
{
  # Nonparametric Regression model with error using first order polynomial basis
  # Model : yi = f(xi) + ei
  # Input
  #   y       : response variable
  #   v       : contaminated explanatory variable
  #   prior   : predefined values of hyperparameters
  #   maxiter : stopping criteria(maximum number of iterations)
  #   tol     : stopping critieria(tolerance level for change of ELBO)
  #   n_grids : level of accuracy of involved numerical integration
  
  N = length(y)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2beta = sig2mu = 100
    asig = bsig = anu = axi = bxi = aups = bups = 1e-3
    bnu = 1e-4
  }
  else
  {
    sig2beta = prior$sig2beta
    sig2mu = prior$sig2mu
    asig = prior$asig
    bsig = prior$bsig
    anu = prior$anu
    axi = prior$axi
    bxi = prior$bxi
    aups = prior$aups
    bups = prior$bups
    bnu = prior$bnu
  }
  
  # Compose interior knots of the basis
  intknots = seq(min(v), max(v), length.out=n_intknot+2)[-c(1,n_intknot+2)] # equally spaced knots
  # intknots = quantile(range(v), seq(0,1,length.out=n_grids+2)[-c(1,n_grids+2)]) # quantile based knots
  
  # Compose spline basis of grids
  grids = seq(min(v)-sd(v)/2, max(v)+sd(v)/2, length.out=n_grids)
  vphig = matrix(grids, n_grids, n_intknot) - outer(rep(1,n_grids), intknots)
  vphig = vphig * (vphig > 0)
  
  # Initialize variational parameters
  asigtl = asig + (N/2)
  anutl = anu + (N/2)
  axitl = axi + (N/2)
  aupstl = aups + (n_intknot/2)
  sig.ratio = asig/bsig
  nu.ratio = anu/bnu
  xi.ratio = axi/bxi
  ups.ratio = aups/bups
  muu.q = rep(0, n_intknot)
  sigu.q = diag(rep(ups.ratio, n_intknot))
  mutl = mean(v)
  mu0tl = 0
  vphiq = matrix(0, N, n_intknot)
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # intercept
    sig02tl = 1/(1/sig2beta+N*sig.ratio)
    mu0tl = sig02tl*sig.ratio*sum(y-vphiq%*%muu.q)
    
    # denoised values
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(muu.q,muu.q)+sigu.q)%*%t(vphig)) + (nu.ratio+xi.ratio)*(grids^2) - 2*xi.ratio*mutl*grids)
    lnpgrids = common + nu.ratio*outer(grids, v) + sig.ratio*outer(drop(vphig%*%muu.q), drop(y-mu0tl)); pgrids = exp(t(lnpgrids))
    normalizers = apply(pgrids, 1, sum)
    ex = drop(pgrids%*%grids)/normalizers
    ex2 = drop(pgrids%*%(grids^2))/normalizers
    varx = ex2 - (ex)^2
    vphiq = pgrids%*%vphig/outer(normalizers,rep(1,n_intknot))
    vphiqtvphiq = t(vphig)%*%diag(apply(pgrids/normalizers,2,sum))%*%vphig
    
    # spline coefficients
    sigu.q = solve(diag(rep(ups.ratio,n_intknot))+sig.ratio*vphiqtvphiq)
    muu.q = drop(sig.ratio*sigu.q%*%t(vphiq)%*%(y-mu0tl))
    
    # sigma
    cpterm = sum((y-mu0tl)^2) + N*sig02tl
    cpterm = cpterm - 2*sum((y-mu0tl)*(vphiq%*%muu.q)) + sum(diag(vphiqtvphiq%*%(sigu.q+outer(muu.q,muu.q))))
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # nu
    bnutl = bnu + 0.5*sum((v-ex)^2+varx)
    nu.ratio = anutl/bnutl
    
    # mu
    sig2mutl = 1/(1/sig2mu+N*xi.ratio)
    mutl = xi.ratio*sig2mutl*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    xi.ratio = axitl/bxitl
    
    # upsilon
    bupstl = bups + 0.5*sum(muu.q^2+diag(sigu.q))
    ups.ratio = aupstl/bupstl
    
    # ELBO
    lbnew = -0.5*N*log(2*pi) - 0.5*N*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*N*log(2*pi) - 0.5*N*(log(bnutl)-digamma(anutl)) - 0.5*nu.ratio*(sum((v-ex)^2)+sum(varx))
    lbnew = lbnew - 0.5*N*(log(bxitl)-digamma(axitl)) - 0.5*xi.ratio*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    lbnew = lbnew + 0.5*sum(log(varx)) + 0.5*N
    lbnew = lbnew - 0.5*log(sig2beta) - 0.5*(mu0tl^2+sig02tl)/sig2beta
    lbnew = lbnew + 0.5*log(sig02tl) + 0.5
    lbnew = lbnew - 0.5*n_intknot*(log(bupstl)-digamma(aupstl)) - 0.5*ups.ratio*sum(muu.q^2+diag(sigu.q))
    lbnew = lbnew + 0.5*determinant(sigu.q,logarithm=TRUE)$modulus[1] + 0.5*n_intknot
    lbnew = lbnew - 0.5*log(sig2mu) - 0.5*(mutl^2+sig2mutl)/sig2mu
    lbnew = lbnew + 0.5*log(sig2mutl) + 0.5
    lbnew = lbnew - lgamma(axi) + axi*log(bxi) - (axi+1)*(log(bxitl)-digamma(axitl)) - xi.ratio*bxi
    lbnew = lbnew + lgamma(axitl) - axitl*log(bxitl) + (axitl+1)*(log(bxitl)-digamma(axitl)) + xi.ratio*bxitl
    lbnew = lbnew - lgamma(anu) + anu*log(bnu) - (anu+1)*(log(bnutl)-digamma(anutl)) - nu.ratio*bnu
    lbnew = lbnew + lgamma(anutl) - anutl*log(bnutl) + (anutl+1)*(log(bnutl)-digamma(anutl)) + nu.ratio*bnutl
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
    lb=lb, ex=ex, varx=varx, mu0tl=mu0tl, sig02tl=sig02tl, muu.q=muu.q, sigu.q=sigu.q,
    sig.ratio=sig.ratio, nu.ratio=nu.ratio, xi.ratio=xi.ratio, ups.ratio=ups.ratio,
    xgrid=xgrid,
    post_curve=post_curve
  ))
}