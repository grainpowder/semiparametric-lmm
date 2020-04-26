mebsar = function(y, w, v, J, prior=NULL, maxiter=1000, tol=1e-4, n_grid=1e3)
{
  library(matrixStats)
  N = length(y)
  D = ncol(w)
  W = cbind(1,w)
  WtW = crossprod(W)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2v = 0.1
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1))
    sig2mu = 100
    axi = bxi = atau = btau = asig = bsig = 1e-3
    w0 = 1
  }
  else
  {
    sig2v = prior$sig2v
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    sig2mu = prior$sig2mu
    axi = prior$axi
    bxi = prior$bxi
    atau = prior$atau
    btau = prior$btau
    asig = prior$asig
    bsig = prior$bsig
    w0 = prior$w0
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  ldetsb0 = as.numeric(determinant(sigbeta.0, logarithm=TRUE)$modulus)
  
  # Grid
  grids = seq(min(v)-sd(v), max(v)+sd(v), length.out=n_grid)
  vphig = sqrt(2)*cos(outer(seq(0,1,length.out=n_grid),pi*(1:J)))
  
  # Initialize variational parameters
  # parametric
  xi.ratio = axi/bxi
  sig.ratio = asig/bsig
  tau.ratio = atau/btau
  axitl = axi + 0.5*N
  asigtl = asig + 0.5*(N+J)
  atautl = atau + 0.5*J
  mutl = mean(v)
  mubeta.q = rep(0, D+1)
  ex = ex2 = rep(0,N)
  # nonparametric
  mutheta.q = rep(0,J)
  sigtheta.q = diag(rep(1,J))
  mupsi.q = 0
  sig2psi.q = 0.01
  DQvec = DQvec_gen(J, mupsi.q, sig2psi.q)
  vphiq = matrix(0, N, J)
  vphiqtvphiq = t(vphiq)%*%vphiq
  
  # ELBO calculation setting
  lb_const = -w0+(J*(J+1))/4
  lb = rep(0, maxiter)
  lbold = -Inf
  
  # Update as
  for (iter in 1:maxiter)
  {
    # denoised value
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(mutheta.q,mutheta.q)+sigtheta.q)%*%t(vphig)) + (1/sig2v+xi.ratio)*grids^2 - 2*xi.ratio*mutl*grids)
    resid = drop(y-W%*%mubeta.q)
    lnpgrids = common + outer(grids, v)/sig2v + sig.ratio*outer(drop(vphig%*%mutheta.q), resid)
    pgrids = exp(lnpgrids)
    normalizers = apply(pgrids, 2, sum)
    ex = drop(t(pgrids)%*%grids / normalizers)
    ex2 = drop(t(pgrids)%*%(grids^2) / normalizers)
    varx = ex2 - ex^2
    vphiq = t(pgrids)%*%vphig / normalizers
    vphiqtvphiq = t(vphig)%*%diag(apply(t(pgrids)/normalizers,2,sum))%*%vphig
    
    # lnpgrids = matrix(0, N, n_grid)
    # for (n in 1:N)
    # {
    #   if (n == 45) browser()
    #   lnpgrid = common + (v[n]/sig2v)*grids + sig.ratio*resid[n]*drop(vphig%*%mutheta.q)
    #   lnpgrids[n,] = lnpgrid
    #   lognorm = logSumExp(lnpgrid)
    #   
    #   gstar = max(log(abs(grids)) + lnpgrid)
    #   ex[n] = exp(log(sum(sign(grids)*exp(log(abs(grids)) + lnpgrid - gstar))) + gstar - lognorm)
    #   ex2[n] = exp(logSumExp(2*log(abs(grids)) + lnpgrid) - lognorm)
    # 
    #   gstar = apply(log(abs(vphig)) + lnpgrid, 2, max)
    #   vphiq[n,] = exp(log(apply(sign(vphig)*exp(log(abs(vphig))+lnpgrid-outer(rep(1,n_grid),gstar)), 2, sum)) + gstar - lognorm)
    # }
    # varx = ex2 - (ex)^2
    # A = colLogSumExps(exp(lnpgrids-rowLogSumExps(lnpgrids)))
    # vphiqtvphiq = t(vphig)%*%diag(A)%*%vphig
    
    # beta
    sigbeta.q = solve(sb0i+sig.ratio*WtW)
    mubeta.q = drop(sigbeta.q%*%(sbimb0+sig.ratio*t(W)%*%(y-vphiq%*%mutheta.q)))
    
    # theta
    sigtheta.q = solve(sig.ratio*(tau.ratio*DQvec+vphiqtvphiq))
    mutheta.q = drop(sig.ratio*sigtheta.q%*%t(vphiq)%*%(y-W%*%mubeta.q))
    
    # mu
    sig2mutl = 1/((1/sig2mu)+N*xi.ratio)
    mutl = sig2mutl*xi.ratio*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    xi.ratio = axitl/bxitl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-vphiq%*%mutheta.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(diag(vphiqtvphiq%*%sigtheta.q))
    trterm3 = sum(diag((outer(mutheta.q,mutheta.q)+sigtheta.q)%*%DQvec))
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2+tau.ratio*trterm3)
    sig.ratio = asigtl/bsigtl
    
    # tau
    btautl = btau + 0.5*sig.ratio*trterm3
    tau.ratio = atautl/btautl
    
    # psi
    psi_update = ncvmp(w0, J, asigtl, bsigtl, atautl, btautl, mutheta.q, sigtheta.q, mupsi.q, sig2psi.q)
    mupsi.q = psi_update$mupsi.q
    sig2psi.q = psi_update$sig2psi.q
    
    # ELBO
    lbnew = -0.5*N*log(2*pi) - 0.5*N*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*(ssterm+trterm1+trterm2)
    lbnew = lbnew - 0.5*N*log(2*pi) - 0.5*N*log(sig2v) - 0.5*(sum((v-ex)^2)+sum(varx))/sig2v
    lbnew = lbnew - 0.5*N*(log(bxitl)-digamma(axitl)) - 0.5*xi.ratio*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    lbnew = lbnew + 0.5*sum(log(varx)) + 0.5*N
    lbnew = lbnew - 0.5*ldetsb0 - 0.5*sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0)) - 0.5*sum(diag(solve(sigbeta.0,sigbeta.q)))
    lbnew = lbnew + 0.5*determinant(sigbeta.q,logarithm=TRUE)$modulus[1] + 0.5*(D+1)
    lbnew = lbnew - 0.5*J*(log(bsigtl)-digamma(asigtl)) - 0.5*J*(log(btautl)-digamma(atautl)) - 0.5*sig.ratio*tau.ratio*trterm3
    lbnew = lbnew + 0.5*determinant(sigtheta.q,logarithm=TRUE)$modulus + 0.5*J
    lbnew = lbnew - 0.5*log(sig2mu) - 0.5*(mutl^2+sig2mutl)/sig2mu
    lbnew = lbnew + 0.5*log(sig2mutl) + 0.5
    lbnew = lbnew - lgamma(axi) + axi*log(bxi) - (axi+1)*(log(bxitl)-digamma(axitl)) - xi.ratio*bxi
    lbnew = lbnew + lgamma(axitl) - axitl*log(bxitl) + (axitl+1)*(log(bxitl)-digamma(axitl)) + xi.ratio*bxitl
    lbnew = lbnew - lgamma(asig) + asig*log(bsig) - (asig+1)*(log(bsigtl)-digamma(asigtl)) - sig.ratio*bsig
    lbnew = lbnew + lgamma(asigtl) - asigtl*log(bsigtl) + (asigtl+1)*(log(bsigtl)-digamma(asigtl)) + sig.ratio*bsigtl
    lbnew = lbnew - lgamma(atau) + atau*log(btau) - (atau+1)*(log(btautl)-digamma(atautl)) - tau.ratio*btau
    lbnew = lbnew + lgamma(atautl) - atautl*log(btautl) + (atautl+1)*(log(btautl)-digamma(atautl)) + tau.ratio*btautl
    lbnew = lbnew + lb_const*Eqabspsi_gen(mupsi.q, sig2psi.q) - 0.5*sig.ratio*tau.ratio*sum(diag((outer(mutheta.q,mutheta.q)+sigtheta.q)%*%DQvec_gen(J, mupsi.q, sig2psi.q)))
    lbnew = lbnew + 0.5*(log(2*pi)+log(sig2psi.q)+1)
    lb[iter] = lbnew
    if (abs(lbnew-lbold)<tol) break
    lbold = lbnew
  }
  lb = lb[1:iter]
  vphi = sqrt(2)*cos(outer(ex,pi*(1:J)))
  post_curve=drop(vphi%*%mutheta.q)
  curve_var=drop(vphi^2%*%diag(sigtheta.q))
  return(list(
    lb=lb, mutheta.q=mutheta.q, sigtheta.q=sigtheta.q,
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    asigtl=asigtl, bsigtl=bsigtl,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve
  ))
}

set.seed(10)
N = 136
D = 5
beta = rnorm(D+1)
RR = 0.8
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 2
w = matrix(rnorm(N*D),N,D)
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
f = function(x) x*sin(2*pi*x)
y = cbind(1,w)%*%beta + f(x) + rnorm(N)
vb_result = mebsar(y,w,v,10)
res = y - cbind(1,w)%*%vb_result$mubeta.q
ord = order(x)
plot(x, res, main="Result", xlab="", ylab="residual")
lines(x[ord], vb_result$post_curve[ord], lwd=3, col=2)
lines(x[ord], vb_result$post_upper[ord], lwd=2, col=3)
lines(x[ord], vb_result$post_lower[ord], lwd=2, col=3)

