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
    sig2v = 0.1
    sig2mu = 100
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1))
    axi = bxi = atau = btau = asig = bsig = 1e-3
  }
  else
  {
    sig2v = prior$sig2v
    sig2mu = prior$sig2mu
    axi = prior$axi
    bxi = prior$bxi
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
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
  vphig = bs(grids, knots=knots, intercept=FALSE)
  K = ncol(vphig)
  
  # Initialize variatonal parameters
  asigtl = asig + (N/2)
  atautl = atau + (K/2)
  axitl = axi + (N/2)
  sig.ratio = asig/bsig
  tau.ratio = atau/btau
  xi.ratio = axi/bxi
  mutheta.q = rep(0, K)
  sigtheta.q = diag(rep(1/tau.ratio, K))
  mutl = mean(v)
  mubeta.q = rep(0,D+1)
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # denoised value
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(mutheta.q,mutheta.q)+sigtheta.q)%*%t(vphig)) + (1/sig2v+xi.ratio)*grids - 2*xi.ratio*mutl*grids)
    lnpgrids = common + outer(grids, v)/sig2v + sig.ratio*outer(drop(vphig%*%mutheta.q), drop(y-W%*%mubeta.q)); pgrids = exp(t(lnpgrids))
    normalizers = apply(pgrids, 1, sum)
    ex = drop(pgrids%*%grids)/normalizers
    ex2 = drop(pgrids%*%(grids^2))/normalizers
    varx = ex2 - (ex)^2
    vphiq = pgrids%*%vphig/outer(normalizers,rep(1,K))
    vphiqtvphiq = t(vphig)%*%diag(apply(pgrids/normalizers,2,sum))%*%vphig
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = drop(sigbeta.q%*%(sbimb0 + sig.ratio*t(W)%*%(y-vphiq%*%mutheta.q)))
    
    # theta
    sigtheta.q = solve(diag(rep(tau.ratio,K)) + sig.ratio*vphiqtvphiq)
    mutheta.q = drop(sig.ratio*sigtheta.q%*%t(vphiq)%*%(y-W%*%mubeta.q))
    
    # mu
    sig2mutl = 1/(1/sig2mu + xi.ratio*N)
    mutl = sig2mutl*xi.ratio*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    xi.ratio = axitl/bxitl
    
    # tau
    btautl = btau + 0.5*sum(mutheta.q^2+diag(sigtheta.q))
    tau.ratio = atautl/btautl
    
    # sigma
    cp_beta = drop(W%*%mubeta.q)
    cp_theta = drop(vphiq%*%mutheta.q)
    cpterm = sum(y^2) + sum(diag(WtW%*%(outer(mubeta.q,mubeta.q)+sigbeta.q))) + sum(diag(vphiqtvphiq%*%(outer(mutheta.q,mutheta.q)+sigtheta.q)))
    cpterm = cpterm - 2*sum(y*cp_beta) - 2*sum(y*cp_theta) + 2*sum(cp_beta*cp_theta)
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # ELBO
    lbnew = -0.5*N*log(2*pi) - 0.5*N*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*N*log(2*pi) - 0.5*N*log(sig2v) - 0.5*(sum((v-ex)^2)+sum(varx))/sig2v
    lbnew = lbnew - 0.5*N*(log(bxitl)-digamma(axitl)) - 0.5*xi.ratio*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    lbnew = lbnew + 0.5*sum(log(varx)) + 0.5*N
    lbnew = lbnew - 0.5*ldetsb0 - 0.5*(sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0))+sum(diag(solve(sigbeta.0,sigbeta.q))))
    lbnew = lbnew + 0.5*determinant(sigbeta.q,logarithm=TRUE)$modulus[1] + 0.5*(D+1)
    lbnew = lbnew - 0.5*K*(log(btautl)-digamma(atautl)) - 0.5*tau.ratio*sum(mutheta.q^2+diag(sigtheta.q))
    lbnew = lbnew + 0.5*determinant(sigtheta.q,logarithm=TRUE)$modulus[1] + 0.5*K
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
  post_curve=drop(vphi%*%mutheta.q)
  curve_var=drop(vphi^2%*%diag(sigtheta.q))
  return(list(
    lb=lb, ex=ex, mutheta.q=mutheta.q, sigtheta.q=sigtheta.q, mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    sig.ratio=sig.ratio, tau.ratio=tau.ratio, xi.ratio=xi.ratio,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve
  ))
}

set.seed(10)
N = 130
D = 6
RR = 0.9
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
# plot(beta,vb_result$mubeta.q,main="Regression coefficiencts")
# lines(-10:10,-10:10)
# ord = order(vb_result$ex)
# plot(x,y,ylab="y",main="Original pattern(dot) vs Denoised pattern(line)")
# lines(vb_result$ex[ord],vb_result$post_curve[ord],lwd=3,col=2)
# lines(vb_result$ex[ord],vb_result$post_lower[ord],lwd=2,col=3)
# lines(vb_result$ex[ord],vb_result$post_upper[ord],lwd=2,col=3)
