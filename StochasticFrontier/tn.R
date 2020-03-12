tn = function(y, w, Z, productivity=TRUE, prior=NULL, tol=1e-2, eps=1e-6, maxiter=500, acc=1000)
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
  real_line = seq(-acc/10, acc/10, len=acc)
  positive_line = seq(eps, acc/10, len=acc)
  
  # Initialize parameters
  asigtl = asig + 0.5*sum(ti)
  mu2coef = -0.5 * ((1/sig02)+N)
  ll2coef = -(alam + 0.5*N + 1)
  sig.ratio = asig / bsig
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  bar_l2 = 1
  bar_l1 = 1
  bar_mu = 0
  
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
    sigi2 = 1/(bar_l2 + sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2*(bar_mu*bar_l1 + sgn*sig.ratio*t(Z)%*%(y-W%*%mubeta.q))
    mui = drop(mui)
    
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE) - pnorm(musigrat, log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2*(1-fres*(musigrat+fres))
    
    
    # mu
    h_mu = function(x) mu2coef*(x^2) + (mu0sig02rat+bar_l1*sum(muu.q))*x - N*pnorm(x,log.p=TRUE)
    # Determining lower(L), upper(U) bound of grids
    lnxi_mu = logSumExp(h_mu(real_line))
    L = -1; while (exp(h_mu(L) - lnxi_mu) > eps) L = L-1
    U = 1 ; while (exp(h_mu(U) - lnxi_mu) > eps) U = U+1
    # Generating adaptive grid as described in C.Ritter(1992)
    grid_candidates = seq(L, U, length.out=acc)
    cdf_mu = cumsum(exp(h_mu(grid_candidates)-lnxi_mu))
    cdf_mu = cdf_mu / cdf_mu[acc]
    probs = seq(eps, 1-eps, length.out=acc)
    grids = rep(0, length=acc)
    for (idx in 1:acc) grids[idx] = grid_candidates[sum(cdf_mu < probs[idx])]
    grids = unique(grids)
    # # Diagnosis
    # plot(grid_candidates, exp(h_mu(grid_candidates) - lnxi_mu), type="l")
    # points(grids, rep(0, length(grids)), col=2, pch=18)
    # Approximate required moment w.r.t mu using griddy Gibbs
    grids_pos = grids[grids > 0]
    grids_neg = grids[grids < 0]
    lnxi_mu = logSumExp(h_mu(grids))
    griddy_pos = exp(logSumExp(log(grids_pos)+h_mu(grids_pos)) - lnxi_mu)
    if (length(grids_pos) == 0) griddy_pos = 0
    griddy_neg = exp(logSumExp(log(-grids_neg)+h_mu(grids_neg)) - lnxi_mu)
    if (length(grids_neg) == 0) griddy_neg = 0
    bar_mu = griddy_pos - griddy_neg
    
    # lambda2
    h_lam2 = function(x) ll2coef*log(x) - (blam+0.5*(sum(muu.q+sigu.q)))*(1/x) + bar_mu*sum(muu.q)*sqrt(1/x)
    # Determining upper(U) bound of grids
    lnxi_lam2 = logSumExp(h_lam2(positive_line))
    U = 1; while (exp(h_lam2(U) - lnxi_lam2) > eps) U = U+1
    # Generating adaptive grid as described in C.Ritter(1992)
    grid_candidates = seq(eps, U, length.out=acc)
    cdf_lam2 = cumsum(exp(h_lam2(grid_candidates)-lnxi_lam2))
    cdf_lam2 = cdf_lam2 / cdf_lam2[acc]
    probs = seq(eps, 1-eps, length.out=acc)
    grids = rep(0, length=acc)
    for (idx in 1:acc) grids[idx] = grid_candidates[sum(cdf_lam2 < probs[idx])]
    grids = unique(grids)
    # # Diagnosis
    # plot(grid_candidates, exp(h_lam2(grid_candidates) - lnxi_lam2), type="l")
    # points(grids, rep(0, length(grids)), col=2, pch=18)
    # Approximate required moment w.r.t lambda2 using griddy Gibbs
    lnxi_l2 = logSumExp(h_lam2(grids))
    bar_l2 = exp(logSumExp(-log(grids)+h_lam2(grids)) - lnxi_l2)
    bar_l1 = exp(logSumExp(-0.5*log(grids)+h_lam2(grids)) - lnxi_l2)
    
    # sigma2
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ti*sigu.q)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = asigtl / bsigtl
    
    # Convergence
    bool1 = sqrt(mean((mubeta.q.old-mubeta.q)^2)) < tol
    bool2 = sqrt(mean((muu.q.old-muu.q)^2)) < tol
    if (bool1 & bool2) break
  }
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl
  ))
}