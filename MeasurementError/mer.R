mer = function(y, v, prior=NULL, maxiter=1000, tol=1e-4, n_grid=1e3)
{
  # Estimate parameters of yi=b0+b1*xi+ei where (vi,yi) is observed instead of (xi,yi)
  N = length(y)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2v = 0.1
    sig2beta = sig2mu = 100
    axi = bxi = asig = bsig = 1e-3
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
  }
  
  # Initialize variational parameters
  mutl = mean(v)
  axitl = axi + (N/2)
  asigtl = asig + (N/2)
  xi.ratio = axi/bxi
  sig.ratio = asig/bsig
  mubeta.q = rep(0, 2)
  sigbeta.q = diag(rep(10, 2))
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # denoised value
    varx = 1/(xi.ratio+(1/sig2v)+sig.ratio*(mubeta.q^2+diag(sigbeta.q))[2])
    ex = varx*(xi.ratio*mutl + v/sig2v + sig.ratio*mubeta.q[2]*(y-mubeta.q[1]))
    EX = cbind(1, ex)
    EXX = matrix(c(N, sum(ex), sum(ex), sum(ex^2+varx)),2,2)
    
    # beta
    sigbeta.q = solve(diag(rep(1/sig2beta,2))+sig.ratio*EXX)
    mubeta.q = drop(sig.ratio*sigbeta.q%*%t(EX)%*%y)
    
    # mu
    sig2mutl = 1/((1/sig2mu)+N*xi.ratio)
    mutl = sig2mutl*xi.ratio*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2)+N*(varx+sig2mutl))
    xi.ratio = axitl/bxitl
    
    # sigma
    cpterm = sum(y^2) - 2*sum(y*(EX%*%mubeta.q)) + sum(diag(EXX%*%(outer(mubeta.q,mubeta.q)+sigbeta.q)))
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # ELBO
    lbnew = -0.5*N*log(2*pi) - 0.5*N*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*N*log(2*pi) - 0.5*N*log(sig2v) - 0.5*(sum((v-ex)^2)+N*varx)/sig2v
    lbnew = lbnew - 0.5*N*(log(bxitl)-digamma(axitl)) - 0.5*xi.ratio*(sum((ex-mutl)^2)+N*(varx+sig2mutl))
    lbnew = lbnew + 0.5*N*log(varx) + 0.5*N
    lbnew = lbnew - log(sig2beta) - 0.5*sum(mubeta.q^2+diag(sigbeta.q))/sig2beta
    lbnew = lbnew + 0.5*determinant(sigbeta.q,logarithm=TRUE)$modulus[1] + 1
    lbnew = lbnew - 0.5*log(sig2mu) - 0.5*(mutl^2+sig2mutl)/sig2mu
    lbnew = lbnew + 0.5*log(sig2mutl) + 0.5
    lbnew = lbnew - lgamma(axi) + axi*log(bxi) - (axi+1)*(log(bxitl)-digamma(axitl)) - xi.ratio*bxi
    lbnew = lbnew + lgamma(axitl) - axitl*log(bxitl) + (axitl+1)*(log(bxitl)-digamma(axitl)) + xi.ratio*bxitl
    lbnew = lbnew - lgamma(asig) + asig*log(bsig) - (asig+1)*(log(bsigtl)-digamma(asigtl)) - sig.ratio*bsig
    lbnew = lbnew + lgamma(asigtl) - asigtl*log(bsigtl) + (asigtl+1)*(log(bsigtl)-digamma(asigtl)) + sig.ratio*bsigtl
    lb[iter] = lbnew
    if (lbnew>lbold & lbnew-lbold<tol) break
    lbold = lbnew
  }
  lb = lb[1:iter]
  return(list(lb=lb, ex=ex, varx=varx, mubeta.q=mubeta.q, sigbeta.q=sigbeta.q))
}

# y = b0 + b1x
set.seed(10)
N = 130
RR = 0.8
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 3
beta = c(0.5, 1.1)
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
y = cbind(1,x)%*%beta + rnorm(N)
result = mer(y,v)
par(mfrow=c(1,2))
plot(result$lb,xlab="Iteration",ylab="",main="Evidence Lower Bound", type="l")
plot(x,result$ex,xlab="True",ylab="Estimate",main="Corrupted variable")
lines(-10:10,-10:10)
result$mubeta.q