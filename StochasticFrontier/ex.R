ex = function(y, w, Z, productivity=TRUE, prior=NULL, tol=1e-4, maxiter=500)
{
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
    alam = blam = asig = bsig = 0.001
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    alam = prior$alam
    blam = prior$blam
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  
  # Initialize parameters
  asigtl = asig + 0.5*sum(ti)
  alamtl = alam + N
  sig.ratio = asig / bsig
  lam.ratio = alam / blam
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  
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
    sigi2 = 1/(sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2*(-lam.ratio + sgn*sig.ratio*t(Z)%*%(y-W%*%mubeta.q))
    mui = drop(mui)
    
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE) - pnorm(musigrat, log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2*(1-fres*(musigrat+fres))
    
    # lambda
    blamtl = blam + sum(muu.q)
    lam.ratio = alamtl / blamtl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ti*sigu.q)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = alamtl / blamtl
    
    # Convergence
    bool1 = mean((mubeta.q.old - mubeta.q)^2) < tol
    bool2 = mean((muu.q.old - muu.q)^2) < tol
    if (bool1 & bool2) break
  }
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl,
    alamtl=alamtl, blamtl=blamtl
  ))
}