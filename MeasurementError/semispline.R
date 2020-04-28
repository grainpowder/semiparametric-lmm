semispline = function(y, w, v, prior=NULL, maxiter=1000, tol=1e-4, n_grids=1e3)
{
  library(splines)
  N = length(y)
  D = ncol(w)
  W = cbind(1,w)
  WtW = crossprod(W)
  
  # Hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1))
    sig2v = 0.1
    sig2mu = sig2alpha = 100
    axi = bxi = atau = btau = asig = bsig = 1e-3
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    sig2v = prior$sig2v
    sig2mu = prior$sig2mu
    sig2alpha = prior$sig2alpha
    axi = prior$axi
    bxi = prior$bxi
    atau = prior$atau
    btau = prior$btau
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  ldetsb0 = determinant(sigbeta.0, logarithm=TRUE)$modulus[1]
  
  # Generate grids and corresponding spline coefficients
  nknots = min(round(log(length(unique(v)))), 35)
  knots = quantile(v, seq(0,1,length.out=nknots+2)[-c(1,(nknots+2))])
  grids = seq(min(v)-sd(v), max(v)+sd(v), length.out=n_grids)
  vphig = bs(grids, knots=knots, intercept=FALSE); K = ncol(vphig); uidx = 3:(K+2)
  vphig = cbind(cbind(1,grids),vphig)
  
  # Initialize variational parameters
  axitl = axi + (N/2)
  atautl = atau + (K/2)
  asigtl = asig + (N/2)
  xi.ratio = axi/bxi
  tau.ratio = atau/btau
  sig.ratio = asig/bsig
  munu.q = munu.0 = rep(0, K+2)
  signu.q = signu.0 = diag(rep(c(sig2alpha,1/tau.ratio),c(2,K)))
  ldetsn0 = determinant(signu.0, logarithm=TRUE)$modulus[1]
  mutl = mean(v)
  mubeta.q = rep(0,D+1)
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # denoised value
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(munu.q,munu.q)+signu.q)%*%t(vphig)) + (1/sig2v+xi.ratio)*(grids^2) - 2*mutl*xi.ratio*grids)
    lnpgrids = common + outer(grids, v)/sig2v + sig.ratio*outer(drop(vphig%*%munu.q), drop(y-W%*%mubeta.q)); pgrids = exp(t(lnpgrids))
    normalizers = apply(pgrids, 1, sum)
    ex = pgrids%*%grids/normalizers
    ex2 = pgrids%*%(grids^2)/normalizers
    varx = ex2 - (ex)^2
    vphiq = pgrids%*%vphig/outer(normalizers,rep(1,K+2))
    vphiqtvphiq = t(vphig)%*%diag(apply(pgrids/normalizers,2,sum))%*%vphig
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = drop(sigbeta.q%*%(sbimb0 + sig.ratio*t(W)%*%(y-vphiq%*%munu.q)))
    
    # nu
    signu.q = solve(diag(rep(c(1/sig2alpha,tau.ratio),c(2,K))) + sig.ratio*vphiqtvphiq)
    munu.q = drop(sig.ratio*signu.q%*%t(vphiq)%*%(y-W%*%mubeta.q))
    
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
    cpterm_beta = drop(W%*%mubeta.q)
    cpterm_nu = drop(vphiq%*%munu.q)
    cpterm = sum(y^2) + sum(diag(WtW%*%(outer(mubeta.q,mubeta.q)+sigbeta.q))) + sum(diag(vphiqtvphiq%*%(outer(munu.q,munu.q)+signu.q)))
    cpterm = cpterm - 2*sum(y*cpterm_beta) - 2*sum(y*cpterm_nu) + 2*sum(cpterm_beta*cpterm_nu)
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # ELBO
    lbnew = -0.5*N*log(2*pi) - 0.5*N*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*N*log(2*pi) - 0.5*N*log(sig2v) - 0.5*(sum((v-ex)^2)+sum(varx))/sig2v
    lbnew = lbnew - 0.5*N*(log(bxitl)-digamma(axitl)) - 0.5*xi.ratio*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    lbnew = lbnew + 0.5*sum(log(varx)) + 0.5*N
    lbnew = lbnew - 0.5*ldetsn0 - 0.5*(sum((munu.q-munu.0)*solve(signu.0,munu.q-munu.0))+sum(diag(solve(signu.0,signu.q))))
    lbnew = lbnew + 0.5*determinant(signu.q,logarithm=TRUE)$modulus[1] + 0.5*(K+2)
    lbnew = lbnew - 0.5*ldetsb0 - 0.5*(sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0))+sum(diag(solve(sigbeta.0,sigbeta.q))))
    lbnew = lbnew + 0.5*determinant(sigbeta.q,logarithm=TRUE)$modulus[1] + 0.5*(D+1)
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
    lb=lb, ex=ex, mubeta.q=mubeta.q, sigbeta.q=sigbeta.q, munu.q=munu.q, signu.q=signu.q, 
    sig.ratio=sig.ratio, tau.ratio=tau.ratio, xi.ratio=xi.ratio,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve
  ))
}

set.seed(10)
N = 130
D = 6
RR = 0.7
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 1.5
beta = rnorm(D+1)
w = matrix(rnorm(N*D),N,D)
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
f = function(x) 2*x+sin(pi*x)
y = cbind(1,w)%*%beta + f(x) + rnorm(N)
vb_result = semispline(y,w,v)
plot(beta,vb_result$mubeta.q,main="Regression coefficiencts")
lines(-10:10,-10:10)
ord = order(vb_result1$ex)
plot(x,y,ylab="y",main="Original pattern(dot) vs Denoised pattern(line)")
lines(vb_result1$ex[ord],vb_result1$post_curve[ord],lwd=3,col=2)
lines(vb_result1$ex[ord],vb_result1$post_lower[ord],lwd=2,col=3)
lines(vb_result1$ex[ord],vb_result1$post_upper[ord],lwd=2,col=3)