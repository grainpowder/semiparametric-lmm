mespline = function(y, v, prior=NULL, maxiter=1000, tol=1e-4, n_grids=1e3)
{
  # Estimate parameters of yi=f(xi)+ei where (vi,yi) is observed instead of (xi,yi)
  library(splines)
  N = length(y)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2v = 0.1
    sig2beta = sig2mu = 100
    axi = bxi = atau = btau = asig = bsig = 1e-3
  }
  else
  {
    sig2v = prior$sig2v
    sig2beta = prior$sig2beta
    sig2mu = prior$sig2mu
    axi = prior$axi
    bxi = prior$bxi
    asig = prior$asig
    bsig = prior$bsig
    atau = prior$atau
    btau = prior$btau
  }
  
  # Generate grids and corresponding spline coefficients
  nknots = min(round(log(length(unique(v)))), 35)
  knots = quantile(v, seq(0,1,length.out=nknots+2)[-c(1,(nknots+2))])
  grids = seq(min(v)-sd(v), max(v)+sd(v), length.out=n_grids)
  vphig = bs(grids, knots=knots, intercept=FALSE); K = ncol(vphig); uidx = 3:(K+2)
  vphig = cbind(cbind(1,grids),vphig)
  
  # Initialize variational parameters
  mutl = mean(v)
  axitl = axi + (N/2)
  asigtl = asig + (N/2)
  atautl = atau + (K/2)
  xi.ratio = axi/bxi
  sig.ratio = asig/bsig
  tau.ratio = atau/btau
  munu.q = munu.0 = rep(0, K+2)
  signu.q = signu.0 = diag(rep(c(sig2beta,1/tau.ratio),c(2,K)))
  ldetsn0 = determinant(signu.0, logarithm=TRUE)$modulus[1]
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # denoised value
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(munu.q,munu.q)+signu.q)%*%t(vphig)) + (1/sig2v+xi.ratio)*(grids^2) - 2*mutl*xi.ratio*grids)
    lnpgrids = common + outer(grids,v)/sig2v + sig.ratio*outer(drop(vphig%*%munu.q),y); pgrids = exp(t(lnpgrids))
    normalizers = apply(pgrids, 1, sum)
    ex = pgrids%*%grids/normalizers
    ex2 = pgrids%*%(grids^2)/normalizers
    varx = ex2 - (ex)^2
    vphiq = pgrids%*%vphig/outer(normalizers,rep(1,K+2))
    vphiqtvphiq = t(vphig)%*%diag(apply(pgrids/normalizers,2,sum))%*%vphig
    
    # nu
    signu.q = solve(diag(rep(c(1/sig2beta,tau.ratio),c(2,K)))+sig.ratio*vphiqtvphiq)
    munu.q = drop(sig.ratio*signu.q%*%t(vphiq)%*%y)
    
    # mu
    sig2mutl = 1/(1/sig2mu+N*xi.ratio)
    mutl = xi.ratio*sig2mutl*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    xi.ratio = axitl/bxitl
    
    # tau
    btautl = btau + 0.5*sum((munu.q^2+diag(signu.q))[uidx])
    tau.ratio = atautl/btautl
    
    # sigma
    cpterm = sum(y^2) - 2*sum(y*drop(vphiq%*%munu.q)) + sum(diag(vphiqtvphiq%*%(outer(munu.q,munu.q)+signu.q)))
    bsigtl = bsig + 0.5*cpterm
    
    # ELBO
    lbnew = -0.5*N*log(2*pi) - 0.5*N*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*N*log(2*pi) - 0.5*N*log(sig2v) - 0.5*(sum((v-ex)^2)+sum(varx))/sig2v
    lbnew = lbnew - 0.5*N*(log(bxitl)-digamma(axitl)) - 0.5*xi.ratio*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    lbnew = lbnew + 0.5*sum(log(varx)) + 0.5*N
    lbnew = lbnew - 0.5*ldetsn0 - 0.5*(sum((munu.q-munu.0)*solve(signu.0,munu.q-munu.0))+sum(diag(solve(signu.0,signu.q))))
    lbnew = lbnew + 0.5*determinant(signu.q,logarithm=TRUE)$modulus[1] + 0.5*(K+2)
    lbnew = lbnew - 0.5*log(sig2mu) - 0.5*(mutl^2+sig2mutl)/sig2mu
    lbnew = lbnew + 0.5*log(sig2mutl) + 0.5
    lbnew = lbnew - lgamma(axi) + axi*log(bxi) - (axi+1)*(log(bxitl)-digamma(axitl)) - xi.ratio*bxi
    lbnew = lbnew + lgamma(axitl) - axitl*log(bxitl) + (axitl+1)*(log(bxitl)-digamma(axitl)) + xi.ratio*bxitl
    lbnew = lbnew - lgamma(atau) + atau*log(btau) - (atau+1)*(log(btautl)-digamma(atautl)) - tau.ratio*btau
    lbnew = lbnew + lgamma(atautl) - atautl*log(btautl) + (atautl+1)*(log(btautl)-digamma(atautl)) + tau.ratio*btautl
    lbnew = lbnew - lgamma(asig) + asig*log(bsig) - (asig+1)*(log(bsigtl)-digamma(asigtl)) - sig.ratio*bsig
    lbnew = lbnew + lgamma(asigtl) - asigtl*log(bsigtl) + (asigtl+1)*(log(bsigtl)-digamma(asigtl)) + sig.ratio*bsigtl
    lb[iter] = lbnew
    if (lbnew>lbold & lbnew-lbold<tol) break
    lbold = lbnew
  }
  lb = lb[1:iter]
  knots = quantile(ex, seq(0,1,length.out=nknots+2)[-c(1,(nknots+2))])
  vphi = bs(ex, knots=knots, intercept=FALSE)
  vphi = cbind(cbind(1,ex),vphi)
  post_curve=drop(vphi%*%munu.q)
  curve_var=drop(vphi^2%*%diag(signu.q))
  return(list(
    lb=lb, ex=ex, munu.q=munu.q, signu.q=signu.q, sig.ratio=sig.ratio, tau.ratio=tau.ratio, xi.ratio=xi.ratio,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve
  ))
}

set.seed(10)
N = 130
RR = 0.8
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 1.5
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
f1 = function(x) 2*x+sin(pi*x)
f2 = function(x) 3*exp(-2*x)
f3 = function(x) exp(0.8*x)
y1 = f1(x) + rnorm(N)
vb_result1 = mespline(y1,v)
ord = order(vb_result1$ex)
plot(x,y1,ylab="y",main="Original pattern vs Denoised pattern")
lines(vb_result1$ex[ord],vb_result1$post_curve[ord],lwd=3,col=2)
lines(vb_result1$ex[ord],vb_result1$post_lower[ord],lwd=2,col=3)
lines(vb_result1$ex[ord],vb_result1$post_upper[ord],lwd=2,col=3)
y2 = f2(x) + rnorm(N)
vb_result2 = mespline(y2,v)
plot(x,y2,ylab="y",main="Original pattern vs Denoised pattern")
lines(vb_result2$ex[ord],vb_result2$post_curve[ord],lwd=3,col=2)
lines(vb_result2$ex[ord],vb_result2$post_lower[ord],lwd=2,col=3)
lines(vb_result2$ex[ord],vb_result2$post_upper[ord],lwd=2,col=3)
y3 = f3(x) + rnorm(N)
vb_result3 = mespline(y3,v)
plot(x,y3,ylab="y",main="Original pattern vs Denoised pattern")
lines(vb_result3$ex[ord],vb_result3$post_curve[ord],lwd=3,col=2)
lines(vb_result3$ex[ord],vb_result3$post_lower[ord],lwd=2,col=3)
lines(vb_result3$ex[ord],vb_result3$post_upper[ord],lwd=2,col=3)