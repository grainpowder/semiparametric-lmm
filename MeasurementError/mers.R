mers = function(y, w, v, prior=NULL, maxiter=500, tol=1e-4, n_grid=1e3)
{
  library(matrixStats)
  library(splines)
  N = length(y)
  D = ncol(w)
  W = cbind(1,w)
  WtW = crossprod(W)
  nknots = min(round(log(length(unique(v)))), 35)
  
  # Hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(1e3, D+1))
    sig2v = 100
    sig2mu = 1e8
    axi = bxi = atau = btau = asig = bsig = 1e-3
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    sig2v = prior$sig2v
    sig2mu = prior$sig2mu
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
  deviance = diff(range(v))
  grid = seq(min(v)-deviance/10, max(v)+deviance/10, length.out=n_grid)
  
  # Construct spline bases on these grids
  knots = quantile(grid, seq(0,1,length=(nknots+2))[-c(1,nknots+2)])
  vphi_g = bs(grid, knots=knots, degree=3, intercept=TRUE) # M*J
  J = ncol(vphi_g)
  
  # Initialize variational parameters
  sig.ratio = tau.ratio = xi.ratio = 1
  mubeta.q = rep(0, D+1)
  mutheta.q = rep(0, J)
  sigtheta.q = diag(rep(10, J))
  asigtl = asig + (N/2)
  atautl = atau + (J/2)
  axitl = axi + (N/2)
  mutl = 0
  vphiq = matrix(0, N, J)
  ex = ex2 = rep(0, N)
  
  # Update as
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    mutheta.q.old = mutheta.q
    
    # denoised variate(x)
    common = -0.5*(sig.ratio*diag(vphi_g%*%(outer(mutheta.q, mutheta.q)+sigtheta.q)%*%t(vphi_g)) + (xi.ratio+(1/sig2v))*grid^2 - 2*mutl*xi.ratio*grid) # M*1
    resid = drop(y-W%*%mubeta.q)
    for (n in 1:N)
    {
      # Assume a>b>c. Then log(exp(a)-exp(b)+exp(c)) = log(1-exp(b-a)+exp(b-c))+a
      pmgrid = common + v[n]*grid/sig2v + sig.ratio*resid[n]*vphi_g%*%mutheta.q # M*1
      pmgrid = drop(pmgrid)
      lognormalizer = logSumExp(pmgrid)
      # normalizer = sum(exp(pmgrid))
      max_exp = max(log(abs(grid))+pmgrid)
      ex[n] = exp(log(sum(sign(grid)*exp(log(abs(grid))+pmgrid-max_exp)))+max_exp - lognormalizer)
      # ex[n] = sum(grid*exp(pmgrid))/normalizer
      ex2[n] = exp(logSumExp(2*log(abs(grid))+pmgrid) - lognormalizer)
      # ex2[n] = sum(grid^2*exp(pmgrid))/normalizer
      colmax = apply(log(abs(vphi_g))+pmgrid,2,max)
      exponent =  sign(vphi_g)*exp(log(abs(vphi_g))+pmgrid-outer(rep(1,n_grid),colmax))
      vphiq[n,] = exp(log(apply(exponent,2,sum))+colmax-lognormalizer)
      # vphiq[n,] = drop(t(vphi_g)%*%exp(pmgrid))/exp(lognormalizer)
    }
    knots = quantile(ex, seq(0,1,length=(nknots+2))[-c(1,nknots+2)])
    vphiq = bs(ex, knots=knots, degree=3, intercept=TRUE) # M*J
    vphitvphiq = t(vphiq)%*%vphiq
    varx = ex2 - (ex)^2
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = sigbeta.q %*% (sbimb0 + sig.ratio*t(W)%*%(y-vphiq%*%mutheta.q))
    mubeta.q = drop(mubeta.q)
    
    # theta
    sigtheta.q = solve(tau.ratio*diag(J) + sig.ratio*vphitvphiq)
    mutheta.q = sig.ratio * sigtheta.q%*%t(vphiq)%*%(y-W%*%mubeta.q)
    mutheta.q = drop(mutheta.q)
    
    # mu
    sig2mutl = 1/((1/sig2mu)+xi.ratio*N)
    mutl = sig2mutl*xi.ratio*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    xi.ratio = axitl/bxitl
    
    # tau
    btautl = btau + 0.5*(sum(mutheta.q^2)+sum(diag(sigtheta.q)))
    tau.ratio = atautl/btautl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-vphiq%*%mutheta.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(diag(vphitvphiq%*%sigtheta.q))
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = asigtl/bsigtl
    
    # Convergence
    bool1 = sqrt(mean((mubeta.q-mubeta.q.old)^2)) < tol
    bool2 = sqrt(mean((mutheta.q-mutheta.q.old)^2)) < tol
    if (bool1&bool2) 
    {
      break
    }
    else
    {
      cat(iter,"\tbeta :", sqrt(mean((mubeta.q-mubeta.q.old)^2)), "\ttheta :", sqrt(mean((mutheta.q-mutheta.q.old)^2)), "\n")
    }
  }
  knots = quantile(ex, seq(0,1,length=(nknots+2))[-c(1,nknots+2)])
  curve = bs(ex, knots=knots, degree=3, intercept=TRUE)%*%mutheta.q
  browser()
  plot(ex,curve)
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
# plot(v,y,main="Noise",pch=19)
# plot(v,y,main="Noise and Overshadowed Pattern",pch=19)
# points(x,y, col=2,pch=19,cex=0.7)

result = mers(y,w,x)
