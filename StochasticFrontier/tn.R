tn = function(y, w, Z, productivity=TRUE, prior=NULL, tol=1e-6, eps=1e-8, maxiter=500, acc=1000)
{
  if (!("matrixStats" %in% row.names(installed.packages()))) stop("Install 'matrixStats' package first")
  library(matrixStats)
  # Fitting Stochastic Frontier for truncated normally distributed inefficiency using VB, GG
  sgn = (-1)^productivity
  N = ncol(Z)
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  ti = apply(Z, 2, sum)
  
  # Hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1)) 
    mu0 = 0
    sig02 = 100 
    alam = blam = asig = bsig = 0.001
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    mu0 = prior$mu0
    sig02 = prior$sig02
    alam = prior$alam
    blam = prior$blam
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  mu0sig02rat = mu0 / sig02
  
  # Initialize parameters
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  bar_mu = mu0
  bar_l2 = alam / blam
  bar_l1 = alam / blam
  sig.ratio = asig / bsig
  asigtl = asig + (sum(ti)/2)
  
  # Update as
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = sigbeta.q %*% (sbimb0 + sig.ratio*t(W)%*%(y-sgn*Z%*%muu.q))
    mubeta.q = drop(mubeta.q)
    
    # inefficiency
    sigi2 = 1 / (bar_l2+sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2*(bar_l1*bar_mu+sgn*sig.ratio*t(Z)%*%(y-W%*%mubeta.q))
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE)-pnorm(musigrat, log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2*(1-fres*(musigrat+fres))
    
    # mu
    # mu is defined on the real line, so both lower(L), upper(U) bound should be determined.
    h_mu = function(x) -0.5*(1/sig02+N)*(x^2) + (mu0sig02rat+bar_l1*sum(muu.q))*x - N*pnorm(x,log.p=TRUE)
    L = -1; while (exp(h_mu(L)) > eps) L = L - 1
    U = 1;  while (exp(h_mu(U)) > eps) U = U + 1
    # Generating adaptive grid as described in C.Ritter(1992)
    grid_candidates = seq(L, U, length.out=acc)
    cdf_mu = cumsum(exp(h_mu(grid_candidates)))
    cdf_mu = cdf_mu / cdf_mu[acc]
    probs = seq(eps, 1-eps, length.out=acc)
    grids = rep(0, length=acc)
    for (idx in 1:acc) grids[idx] = grid_candidates[sum(cdf_mu <= probs[idx])]
    grids = unique(grids)
    # Approximate the expectation using griddy Gibbs
    grids_pos = grids[grids > 0]
    grids_neg = grids[grids < 0]
    griddy_pos = exp(logSumExp(log(grids_pos)+h_mu(grids_pos)) - logSumExp(h_mu(grids_pos)))
    if (length(grids_pos) == 0) griddy_pos = 0
    griddy_neg = exp(logSumExp(log(-grids_neg)+h_mu(grids_neg)) - logSumExp(h_mu(grids_neg)))
    if (length(grids_neg) == 0) griddy_neg = 0
    bar_mu = griddy_pos - griddy_neg
    
    # lambda
    # Determine upper bound of the integral
    h_lambda = function(x) (alam+(N/2)-1)*log(x) - (blam+0.5*sum(muu.q^2+sigu.q))*x + bar_mu*sum(muu.q)*sqrt(x)
    U = 1; while (exp(h_lambda(U)) > eps) U = U + 1
    # Generating adaptive grid as described in C.Ritter(1992)
    grid_candidates = seq(eps, U, length.out=acc)
    cdf_lambda = cumsum(exp(h_lambda(grid_candidates)))
    cdf_lambda = cdf_lambda / cdf_lambda[acc]
    probs = seq(eps, 1-eps, length.out=acc)
    grids = rep(0, length=acc)
    for (idx in 1:acc) grids[idx] = grid_candidates[sum(cdf_lambda <= probs[idx])]
    grids = unique(grids)
    # Approximate the expectation using griddy Gibbs
    bar_l2 = exp(logSumExp(log(grids)+h_lambda(grids)) - logSumExp(h_lambda(grids)))
    bar_l1 = exp(logSumExp(0.5*log(grids)+h_lambda(grids)) - logSumExp(h_lambda(grids)))
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ti*sigu.q)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    
    # Convergence
    bool1 = mean((mubeta.q.old-mubeta.q)^2) < tol
    bool2 = mean((muu.q.old-muu.q)^2) < tol
    if (bool1 & bool2) break
  }
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl
  ))
}

source("../misc/make_Z.R")
library(truncnorm)
set.seed(10)
N = 50
Z = make_Z(rep(4, N))
D = 10
w = matrix(rnorm(D*nrow(Z)), ncol=D)
u = rtruncnorm(N, a=0, mean=-1)
beta = rnorm(D+1)
y = cbind(1, w)%*%beta + ((-1)^TRUE)*Z%*%u + rnorm(nrow(Z))
result = tn(y,w,Z)

plot(beta,result$mubeta.q)
lines(-10:10,-10:10)
plot(u,result$muu.q)
lines(-10:10,-10:10)
