mer = function(y, w, v, prior=NULL, maxiter=500, tol=1e-4, n_grid=1e3)
{
  # Fit yi = Wi*beta + gamma*xi where vi~N(xi, sig2v)
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
  ex = v
  mutl = mean(v)
  varx = var(v)
  asigtl = asig + (N/2)
  axitl = axi + (N/2)
  gamtl = 0
  sig.ratio = asig/bsig
  xi.ratio = axi/bxi
  lb = rep(0, maxiter)
  lbold = -Inf
  
  mubeta.q = rep(0,D+1)
  sig2gamtl = 10
  # Update as
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    ex.old = ex
    
    # denoised value(x)
    varx = 1/(xi.ratio+(1/sig2v)+sig.ratio*(gamtl^2+sig2gamtl))
    ex = varx*(xi.ratio*mutl+v/sig2v+sig.ratio*gamtl*drop(y-W%*%mubeta.q))
    
    # mu
    sig2mutl = 1/(1/sig2mu+xi.ratio*N)
    mutl = sig2mutl*xi.ratio*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2)+N*(sig2mutl+varx))
    xi.ratio = axitl/bxitl
    
    # beta
    sigbeta.q = solve(diag(rep(1/sig2beta,D+1)) + sig.ratio*WtW)
    mubeta.q = sig.ratio*sigbeta.q%*%t(W)%*%(y-gamtl*ex)
    mubeta.q = drop(mubeta.q)
    
    # gamma
    sig2gamtl = 1/(1/sig2gam+sig.ratio*sum(ex^2+varx))
    gamtl = sig2gamtl*sig.ratio*sum(ex*drop(y-W%*%mubeta.q))
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-gamtl*ex)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = (gamtl^2+sig2gamtl)*N*varx
    trterm3 = sig2gamtl*sum(ex^2)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2+trterm3)
    sig.ratio = asigtl/bsigtl
    
    # ELBO
    lbnew = -0.5*log(2*pi) - (N/2)*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*(ssterm+trterm1+trterm2+trterm3)
    lbnew = lbnew - 0.5*log(2*pi) - (N/2)*log(sig2v) - (sum((v-ex)^2)+N*varx)/(2*sig2v)
    lbnew = lbnew - 0.5*N*(log(bxitl)-digamma(axitl)) - 0.5*xi.ratio*(sum((ex-mutl)^2) + N*(varx+sig2mutl))
    lbnew = lbnew + 0.5*N*(varx+1)
    lbnew = lbnew - 0.5*(D+1)*log(sig2beta) - sum(mubeta.q^2+diag(sigbeta.q))/(2*sig2beta)
    lbnew = lbnew + 0.5*(determinant(sigbeta.q)$modulus[1]+D+1)
    lbnew = lbnew - 0.5*log(sig2gam) - 0.5*(gamtl^2+sig2gamtl)/sig2gam + 0.5*(log(sig2gamtl)+1)
    lbnew = lbnew - 0.5*log(sig2mu) - 0.5*(mutl^2+sig2mutl)/sig2mu + 0.5*(log(sig2mutl)+1)
    lbnew = lbnew + lgamma(axitl)-lgamma(axi)+axitl*log(bxitl)-axi*log(bxi)+(axitl+1)*(log(bxitl)-digamma(axitl))-(axi+1)*(log(bxi)-digamma(axi))+axitl-bxi*xi.ratio
    lbnew = lbnew + lgamma(asigtl)-lgamma(asig)+asigtl*log(bsigtl)-asig*log(bsig)+(asigtl+1)*(log(bsigtl)-digamma(asigtl))-(asig+1)*(log(bsig)-digamma(asig))+asigtl-bsig*sig.ratio
    
    # Convergence
    if (abs(lbnew-lbold) < tol) break
    lbold = lb[iter] = lbnew
    mse1 = sqrt(mean((mubeta.q-mubeta.q.old)^2))
    mse2 = sqrt(mean((ex-ex.old)^2))
    cat("\nbeta :",mse1, "\tx :", mse2)
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
# for (idx in 1:D) plot(w[,idx],y,main=paste("idx =",idx),xlab="")
# plot(v,y,main="Noise and Signal",pch=19)
# points(x,y,cex=0.7,col=2,pch=19)
result = mer(y,w,v)
