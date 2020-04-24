mer = function(y, w, v, prior=NULL, maxiter=500, tol=1e-4, n_grid=1e3)
{
  # Fit yi = Wi*beta + gamma*xi where vi~N(xi, sig2v)
  # Main point : order of the optimization
  
  N = length(y)
  D = ncol(w)
  W = cbind(1,w)
  WtW = crossprod(W)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2v = var(v)
    sig2beta = sig2gam = sig2mu = 100
    axi = bxi = asig = bsig = 1e-3
  }
  else
  {
    sig2v = prior$sig2v
    sig2beta = prior$sig2beta
    sig2gam = prior$sig2gam
    sig2mu = prior$sig2mu
    axi = prior$axi
    bxi = prior$bxi
    asig = prior$asig
    bsig = prior$bsig
  }
  
  # Initialize variational parameters
  mubeta.q = solve(crossprod(W))%*%t(W)%*%y
  ex = v
  varx = var(v)
  axitl = axi + (N/2)
  asigtl = asig + (N/2)
  xi.ratio = axi/bxi
  sig.ratio = asig/bsig
  mutl = mean(v)
  gamtl = 0
  sig2gamtl = 100
  
  # Update as
  for (iter in 1:maxiter)
  {
    
    mubeta.q.old = mubeta.q
    ex.old = ex
    
    # beta
    sigbeta.q = solve(diag(rep(1/sig2beta,D+1)) + sig.ratio*WtW)
    mubeta.q = sig.ratio*sigbeta.q%*%t(W)%*%(y-gamtl*ex)
    mubeta.q = drop(mubeta.q)
    
    # gamma
    sig2gamtl = 1/(1/sig2gam+sig.ratio*sum(ex^2+varx))
    gamtl = sig2gamtl*sig.ratio*sum(ex*(y-W%*%mubeta.q))
    
    # mu
    sig2mutl = 1/(N*xi.ratio+1/sig2mu)
    mutl = sig2mutl*xi.ratio*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2)+N*varx+N*sig2mutl)
    xi.ratio = axitl/bxitl
    
    # sigma
    bsigtl = bsig + 0.5*(sum((y-W%*%mubeta.q-gamtl*ex)^2)+sum(diag(WtW%*%sigbeta.q))+sum(ex^2+varx)*sig2gamtl)
    sig.ratio = asigtl/bsigtl
    
    # denoised value
    varx = 1/(xi.ratio+(1/sig2v)+sig.ratio*(gamtl^2+sig2gamtl))
    ex = varx*(xi.ratio*mutl + v/sig2v + sig.ratio*gamtl*(y-W%*%mubeta.q))
    ex = drop(ex)
    
    # Convergence
    mse1 = sqrt(mean((mubeta.q-mubeta.q.old)^2))
    mse2 = sqrt(mean((ex-ex.old)^2))
    if (mse1<tol & mse2<tol) {
      break
    } else {
      cat("\nbeta :", mse1, "\tx :", mse2)
    }
  }
  browser()
}

# y = wb + f(x)
set.seed(10)
N = 130; D = 5
beta = rnorm(D+1)
gam = 1.1
RR = 0.7
xi2 = 0.9
sig2v = xi2/RR-xi2
mux = 3
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
w = matrix(rnorm(N*D), N,D)
y = cbind(1,w)%*%beta + gam*x + rnorm(N)
for (idx in 1:D) plot(w[,idx],y,main=paste("idx =",idx),xlab="")
plot(v,y,main="Noise and Signal",pch=19)
points(x,y,cex=0.7,col=2,pch=19)
result = mer(y,w,x)
