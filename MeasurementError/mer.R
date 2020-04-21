mer = function(y, w, v, prior=NULL, maxiter=500, tol=1e-4, n_grid=1e3)
{
  library(matrixStats)
  library(splines)
  N = length(y)
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  
  # Define hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1))
    sigv2 = sigmu2 = 1e2
    axi = bxi = atau = btau = asig = bsig = 1e-3
  } else {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    sigv2 = prior$sigv2
    sigmu2 = prior$sigmu2
    axi = prior$axi
    bxi = prior$bxi
    atau = prior$atau
    btau = prior$btau
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  
  # Define grids
  grids = seq(min(v), max(v), length.out=n_grid)
  knots = quantile(grids, seq(0,1,length.out=20)[-c(1,20)])
  vphi_g = bs(grids, knots=knots, intercept=FALSE) # M*J
  J = ncol(vphi_g)
  
  # Initialize hyperparameters
  knots = quantile(v, seq(0,1,length.out=20)[-c(1,20)])
  vphiq = bs(v, knots=knots, intercept=FALSE)
  vphitvphiq = crossprod(vphiq)
  mubeta.q = rep(0, D+1)
  mutheta.q = rep(0, J)
  sigtheta.q = diag(rep(1, J))
  sig.ratio = asig/bsig
  tau.ratio = atau/btau
  xi.ratio = axi/bxi
  axitl = axi+(N/2)
  atautl = atau+(J/2)
  asigtl = asig+(N/2)
  mutl = 0
  
  # Update as
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    mutheta.q.old = mutheta.q
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = sigbeta.q%*%(sbimb0 + sig.ratio*t(W)%*%(y-vphiq%*%mutheta.q))
    mubeta.q = drop(mubeta.q)
    
    # denoised value(x)
    pmgrid = sig.ratio*vphi_g%*%(mutheta.q^2+diag(sigtheta.q)) + (xi.ratio+(1/sigv2))*(grids^2) - 2*xi.ratio*mutl*grids
    pmgrid = -0.5*drop(pmgrid) + outer(grids,v)/sigv2 + sig.ratio*outer(drop(vphi_g%*%mutheta.q), drop(y-W%*%mubeta.q))
    pmgrid = t(pmgrid) # N*M
    
    qx = exp(colLogSumExps(log(grids)+t(pmgrid)) - rowLogSumExps(pmgrid))
    qx2 = exp(colLogSumExps(2*log(grids)+t(pmgrid)) - rowLogSumExps(pmgrid))
    qvarx = qx2 - (qx)^2
    if (sum(qvarx<0) != 0) browser()
    
    vphiq = exp(pmgrid)%*%vphi_g / exp(outer(rowLogSumExps(pmgrid), rep(1,J)))
    vphitvphiq = t(vphiq) %*% vphiq
    
    # theta
    sigtheta.q = tau.ratio*diag(J) + sig.ratio*vphitvphiq
    mutheta.q = sig.ratio*sigtheta.q%*%t(vphiq)%*%(y-W%*%mubeta.q)
    mutheta.q = drop(mutheta.q)
    
    # mu
    sigmu2tl = 1/(1/sigmu2+N*xi.ratio)
    mutl = xi.ratio*sigmu2tl*sum(qx)
    
    # xi
    bxitl = bxi + 0.5*(sum((qx-mutl)^2) + sum(qvarx) + N*sigmu2tl)
    xi.ratio = axitl / bxitl
    
    # tau
    btautl = btau + 0.5*sum(mutheta.q^2+diag(sigtheta.q))
    tau.ratio = atautl / btautl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-vphiq%*%mutheta.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(diag(vphitvphiq%*%sigtheta.q))
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = asigtl / bsigtl
    
    # Convergence
    bool1 = if (mean((mubeta.q.old-mubeta.q)^2) < tol)
    bool2 = if (mean((mutheta.q.old-mutheta.q)^2) < tol)
    if (bool1&bool2) break
  }
  print(iter)
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q, qx=qx
  ))
}

# y = wb + f(x)
set.seed(10)
N = 130; D = 5
beta = rnorm(D+1)
RR = 0.8
psi2 = 0.9
sig2v = psi2/RR-psi2
mux = 3
x = rnorm(N, mux, sqrt(psi2))
v = rnorm(N, x, sqrt(sig2v))
w = matrix(rnorm(N*D), N,D)
f = function(x) x*sin(pi*x)
y = cbind(1,w)%*%beta + f(x) + rnorm(N)
plot(v,y,main="Noise",pch=19)
plot(v,y,main="Noise and Overshadowed Pattern",pch=19)
points(x,y, col=2,pch=19,cex=0.7)

result = mer(y,w,v)