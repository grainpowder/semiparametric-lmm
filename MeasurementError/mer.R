mer = function(y, w, v, prior=NULL, maxiter=500, tol=1e-4, n_grid=1e3, nknots=10)
{
  library(splines)
  library(matrixStats)
  N = length(y)
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  
  # Define hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1))
    sigv2 = 1
    sigmu2 = 100
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
  err = diff(range(v))
  grids = seq(min(v)-err/10, max(v)+err/10, length.out=n_grid)
  knots = quantile(grids, seq(0,1,length.out=nknots)[-c(1,nknots)])
  vphi_g = bs(grids, knots=knots, intercept=FALSE) # M*J
  J = ncol(vphi_g)
  
  # Initialize hyperparameters
  knots = quantile(v, seq(0,1,length.out=nknots)[-c(1,nknots)])
  vphiq = bs(v, knots=knots, intercept=FALSE)
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
  qx = qx2 = qvarx = rep(0,N)
  
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
    common = sig.ratio*diag(vphi_g%*%(outer(mutheta.q,mutheta.q)+sigtheta.q)%*%t(vphi_g)) + (xi.ratio+1/sigv2)*(grids^2) - 2*xi.ratio*mutl*grids
    resid = drop(y-W%*%mubeta.q)
    for (n in 1:N)
    {
      pmgrid = -0.5*common + v[n]*grids/sigv2 + sig.ratio*resid[n]*drop(vphi_g%*%mutheta.q)
      lognorm = logSumExp(pmgrid)
      gstar = max(log(abs(grids))+pmgrid)
      colstar = apply(abs(vphi_g)+pmgrid,2,max)
      
      qx[n] = exp(log(sum(sign(grids) * exp(log(abs(grids))+pmgrid-gstar))) + gstar - lognorm)
      qx2[n] = exp(logSumExp(2*log(abs(grids))+pmgrid) - lognorm)
      vphiq[n,] = exp(log(apply(sign(vphi_g) * exp(log(abs(vphi_g))+pmgrid-outer(rep(1,n_grid), colstar)),2,sum)) + colstar - lognorm)
    }
    qvarx = qx2-(qx)^2
    vphitvphiq = t(vphiq)%*%vphiq
    
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
    bool1 = mean((mubeta.q.old-mubeta.q)^2) < tol
    bool2 = mean((mutheta.q.old-mutheta.q)^2) < tol
    if (bool1&bool2) {
      break
    } else {
      print("")
      print(mean((mubeta.q.old-mubeta.q)^2))
      print(mean((mutheta.q.old-mutheta.q)^2))
      print("")
    }
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

result = mer(y,w,x)