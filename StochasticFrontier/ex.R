times = function(x, y) x*y
revcumsum = function(x) {rev(cumsum(rev(x)))}


# Exponential inefficiency ------------------------------------------------

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


# BSA ---------------------------------------------------------------------

ex_bsa = function(y, x, w, Z, J, productivity=TRUE, prior=NULL, maxiter=500, tol=1e-4)
{
  # Pre-calculate objects to be used frequently
  sgn = (-1)^productivity
  N = ncol(Z)
  D = ncol(w)
  W = cbind(1,w)
  WtW = crossprod(W)
  ti = apply(Z,2,sum)
  # zero-one standardization of x
  minx = min(x)
  maxx = max(x)
  if(minx < 0 | maxx > 1) x = (x-minx)/(maxx-minx)
  
  # Hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1))
    alam = blam = asig = bsig = atau = btau = 1e-2 
    w0 = 2
    mupsi.q.start = 1
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    alam = prior$alam
    blam = prior$blam
    asig = prior$asig
    bsig = prior$bsig
    atau = prior$atau
    btau = prior$btau
    w0 = prior$w0
    mupsi.q.start = prior$mupsi.q.start
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  ldetsb0 = as.numeric(determinant(sigbeta.0, logarithm=TRUE)$modulus)
  
  # Initialize variational parameters
  # parametric
  mubeta.q = mubeta.0
  muu.q = rep(0, N)
  sig.ratio = asig/bsig
  tau.ratio = atau/btau
  lam.ratio = alam/blam
  alamtl = alam + N
  # nonparametric
  mutheta.q = rep(0,J)
  mupsi.q = mupsi.q.start
  mu2psi.q = mupsi.q^2
  sig2psi.q = mu2psi.q/100
  sigpsi.q = sqrt(sig2psi.q)
  Jfull = J
  vphifull = sqrt(2)*cos(outer(x,pi*(1:J)))
  vphitvphifull = t(vphifull)%*%vphifull
  vphi = vphifull[,1:J]
  vphitvphi = vphitvphifull[1:J,1:J]
  
  # ELBO calculation setting
  lbold = -Inf
  dif = tol+1
  lb = c()
  
  # Iteration
  for (iter in 1:maxiter)
  {
    mubeta.q.old=mubeta.q
    muu.q.old=muu.q
    
    # theta
    J.old = J
    J = min(floor(-15/(-0.5*mupsi.q)),Jfull)
    mutheta.q.old = mutheta.q[1:J]
    asigtl = asig+(sum(ti)+J)*0.5
    atautl = atau+(J/2)
    bindices = (1:J)
    bindices2= bindices^2
    mfactor = -J*(J+1)/(4*w0)
    
    # Set up spectral design matrix
    vphi = vphifull[,1:J]
    vphitvphi = vphitvphifull[1:J,1:J]
    
    # Update variational distribution parameters for 
    # theta
    term1 = exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2 = exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    sigpsi.q = sqrt(sig2psi.q)
    Qvec = term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec = Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices))
    DQvec = diag(Qvec)
    
    sigtheta.q = sig.ratio*(tau.ratio*DQvec+vphitvphi)
    sigtheta.q = solve(sigtheta.q)
    mutheta.q = sig.ratio*sigtheta.q%*%t(vphi)%*%(y-W%*%mubeta.q-sgn*Z%*%muu.q)
    mutheta.q = drop(mutheta.q)
    
    # beta
    sigbeta.q = solve(sb0i+sig.ratio*WtW)
    mubeta.q = sigbeta.q%*%(sbimb0+sig.ratio*t(W)%*%(y-sgn*Z%*%muu.q-vphi%*%mutheta.q))
    mubeta.q = drop(mubeta.q)
    
    # inefficiency
    sigi2 = 1/(sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2*(-lam.ratio+sgn*sig.ratio*t(Z)%*%(y-W%*%mubeta.q-vphi%*%mutheta.q))
    mui = drop(mui)
    musigrat = mui/sigi
    fres = exp(dnorm(musigrat,log=TRUE)-pnorm(musigrat,log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2*(1-fres*(musigrat+fres))
    
    # lambda
    blamtl = blam+sum(muu.q)
    lam.ratio = alamtl/blamtl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q-vphi%*%mutheta.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ti*sigu.q)
    trterm3 = sum(diag(vphitvphi%*%sigtheta.q))
    trterm4 = sum(diag((outer(mutheta.q,mutheta.q)+sigtheta.q)%*%DQvec))
    bsigtl = bsig+0.5*(ssterm+trterm1+trterm2+trterm3+tau.ratio*trterm4)
    sig.ratio = asigtl/bsigtl
    
    # tau
    btautl = btau+0.5*sig.ratio*trterm4
    tau.ratio = atautl/btautl
    lntau = log(btautl) - digamma(atautl)
    lnsig = log(bsigtl) - digamma(asigtl)
    lnlam = -log(blamtl) + digamma(alamtl)
    
    # ELBO calculation(before NCVMP)
    lbnew = 0
    sigpsi.q = sqrt(sig2psi.q)
    mu2psi.q = mupsi.q^2
    const1 = sqrt(2/pi)
    # likelihood
    lbnew = lbnew - (sum(ti)/2*log(2*pi)) - (sum(ti)/2)*lnsig - sig.ratio*(ssterm+trterm1+trterm2+trterm3)/2
    # beta
    lbnew = lbnew - (ldetsb0 - determinant(sigbeta.q,logarithm=TRUE)$modulus)/2 - (sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0))+sum(diag(solve(sigbeta.0,sigbeta.q)))-(D+1))/2
    # inefficiency
    lbnew = lbnew + sum(lnlam - lam.ratio*muu.q + (log(2*pi)+log(sigu.q)+1)/2 + log(1-pnorm(0,muu.q,sqrt(sigu.q))))
    # theta
    lbnew = lbnew - J*lnsig/2 - J*lntau/2 - sig.ratio*tau.ratio*trterm4/2 + determinant(sigtheta.q,logarithm=TRUE)$modulus/2 + J/2
    # lambda
    lbnew = lbnew + (-lgamma(alam)+lgamma(alamtl)) + alam*log(blam)-alamtl*log(blamtl) + lnlam*(alam-alamtl) - lam.ratio*(blam-blamtl)
    # sigma
    lbnew = lbnew + (-lgamma(asig)+lgamma(asigtl)) + asig*log(bsig)-asigtl*log(bsigtl) - lnsig*(asig-asigtl) - sig.ratio*(bsig-bsigtl)
    # tau
    lbnew = lbnew + (-lgamma(atau)+lgamma(atautl)) + atau*log(btau)-atautl*log(btautl) - lntau*(atau-atautl) - tau.ratio*(btau-btautl)
    # psi(tentative)
    lbnew.wopsibits = lbnew
    lbnew = lbnew + 0.5 * (log(2*pi)+log(sig2psi.q)+1)
    S1 =  -w0 * (sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q)) + mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    lbnew = lbnew + T*(T+1)/4 * (sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    Qvec = term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec = Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)) 
    lbnew = lbnew - 0.5*tau.ratio*sig.ratio*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(Qvec)))
    lbnew = lbnew + log(w0/2) + S1
    lbfull = lbnew
    
    # NCVMP
    sig3psi.q = sigpsi.q^3
    sig4psi.q = sigpsi.q^4
    mu2psi.q = mupsi.q^2
    const1 = sqrt(2/pi)
    const2 = exp(-mu2psi.q/(2*sig2psi.q))
    pd1sig = -w0*((1/(2*sig2psi.q)+sigpsi.q*mu2psi.q/(2*sig4psi.q))*const1*const2-mu2psi.q/sig3psi.q*dnorm(-mupsi.q/sigpsi.q))
    pd1mu = -w0*(-mupsi.q/sigpsi.q*const1*const2+(1-2*pnorm(-mupsi.q/sigpsi.q))+2*mupsi.q/sigpsi.q*dnorm(-mupsi.q/sigpsi.q))
    
    term1 = exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2 = exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    term3 = mupsi.q/(2*sig3psi.q)-bindices/(2*sig2psi.q)
    term4 = mupsi.q/(2*sig3psi.q)+bindices/(2*sig2psi.q)
    term5 = dnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term6 = dnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term7 = 1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term8 = 1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)
    DQmu = term1*term5*1/sigpsi.q+bindices*term1*term7-term2*term6*1/sigpsi.q-bindices*term2*term8    
    DQsig = -term1*term5*term3+bindices2/2*term1*term7+term2*term6*term4+bindices2/2*term2*term8    
    pd2sig = mfactor*pd1sig-0.5*sig.ratio*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*DQsig)
    pd2mu = mfactor*pd1mu-0.5*sig.ratio*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*DQmu)
    sig2psi.q.old = sig2psi.q
    sig2psi.q = -0.5/(pd1sig+pd2sig)
    mupsi.q.old = mupsi.q
    
    # ELBO calculation for psi
    sigpsi.q = sqrt(sig2psi.q)
    mu2psi.q = mupsi.q^2
    term1 = exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2 = exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    lbnew = lbnew.wopsibits
    lbnew = lbnew-(-0.5*log(2*pi)-0.5*log(sig2psi.q)-0.5)
    S1 =  -w0*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    lbnew = lbnew+T*(T+1)/4*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    Qvec = term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec = Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)) 
    lbnew = lbnew-0.5*tau.ratio*sig.ratio*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(Qvec)))
    lbnew = lbnew+log(w0/2)+S1
    dif = lbnew-lbfull
    
    # Modify step size until ELBO increases
    if(dif < 0)
    {
      step = 1
      dif.try = dif
      while(dif.try < 0) 
      {
        step = step*0.5
        sig2psi.q.try = 1/(1/sig2psi.q.old+step*(1/sig2psi.q-1/sig2psi.q.old))
        sigpsi.q.try = sqrt(sig2psi.q.try)
        mupsi.q.try = sig2psi.q.try*(mupsi.q.old/sig2psi.q.old+step*(mupsi.q/sig2psi.q-mupsi.q.old/sig2psi.q.old))
        mu2psi.q.try = mupsi.q.try^2
        term1 = exp(sig2psi.q.try*bindices2/2+mupsi.q.try*bindices)
        term2 = exp(sig2psi.q.try*bindices2/2-mupsi.q.try*bindices)
        lbnew = lbnew.wopsibits
        lbnew = lbnew-(-0.5*log(2*pi)-0.5*log(sig2psi.q.try)-0.5)
        S1 =  -w0*(sigpsi.q.try*const1*exp(-mu2psi.q.try/(2*sig2psi.q.try))+mupsi.q.try*(1-2*pnorm(-mupsi.q.try/sigpsi.q.try)))
        lbnew = lbnew+T*(T+1)/4*(sigpsi.q.try*const1*exp(-mu2psi.q.try/(2*sig2psi.q.try))+mupsi.q.try*(1-2*pnorm(-mupsi.q.try/sigpsi.q.try)))
        Qvec = term1*(1-pnorm(-mupsi.q.try/sigpsi.q.try-sigpsi.q.try*bindices))
        Qvec = Qvec+term2*(1-pnorm(mupsi.q.try/sigpsi.q.try-sigpsi.q.try*bindices)) 
        lbnew = lbnew-0.5*tau.ratio*sig.ratio*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(Qvec)))
        lbnew = lbnew+log(w0/2)+S1
        dif.try = lbnew-lbfull
      }
      sigpsi.q = sigpsi.q.try
      sig2psi.q = sig2psi.q.try
      mupsi.q = mupsi.q.try
      mu2psi.q = mu2psi.q.try
    }    
    
    dif = lbnew-lbold
    dif = dif/abs(lbnew)
    lbold = lbnew
    lb = c(lb,lbnew)
    if (dif < tol) break
  }
  post_curve=drop(vphi%*%mutheta.q)
  curve_var=drop(vphi^2%*%diag(sigtheta.q))
  return(list(
    lb=lb,
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl,
    alamtl=alamtl, blamtl=blamtl,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve
  ))
}

# Exponential mixture inefficiency ----------------------------------------

exmixture = function(y, w, Z, R=10, productivity=TRUE, prior=NULL, maxiter=500, eps=1e-4)
{
  library(matrixStats)
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
    alam = blam = aalp = balp = asig = bsig = 1e-3
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    alam = prior$alam
    blam = prior$blam
    aalp = prior$aalp
    balp = prior$balp
    asig = prior$asig
    bsig = prior$bsig
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  aalptl = aalp + R - 1
  asigtl = asig + 0.5*(sum(ti))
  
  # Initialize parameters
  # Level1
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  sig.ratio = asig / bsig
  # Level2
  kappa = matrix(1,N,R) / R
  kappa_r = apply(kappa,2,sum)
  lam.ratio = rep(alam/blam, R)
  # Level3 or above
  gam1s = rep(0.01, R-1)
  gam2s = rep(0.01, R-1)
  alp.ratio = aalp / balp
  logstick = digamma(gam1s) - digamma(gam1s+gam2s)
  log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
  qstick = c(logstick, 0) + c(0, log1mstick)
  
  # Update as
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
    # beta
    sigbeta.q = solve(sb0i+sig.ratio*WtW)
    mubeta.q = sigbeta.q%*%(sbimb0+sig.ratio*t(W)%*%(y-sgn*Z%*%muu.q))
    mubeta.q = drop(mubeta.q)
    
    # inefficiency
    sigi2 = 1/(sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2 * (-apply(t(kappa)*lam.ratio,2,sum) + sgn*sig.ratio*drop(t(Z)%*%(y-W%*%mubeta.q)))
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE)-pnorm(musigrat, log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2 * (1-fres*(musigrat+fres))
    
    # lambda_r
    alamtl = alam + kappa_r
    blamtl = blam + apply(kappa*muu.q, 2, sum)
    lam.ratio = alamtl / blamtl
    
    # kappa
    S = outer(muu.q, lam.ratio, times)
    S = t(qstick - log(blamtl) + digamma(alamtl) + t(S))
    kappa = exp(kappa - rowLogSumExps(kappa))
    kappa_r = apply(kappa, 2, sum)
    
    # stick length(Vr)
    gam1s = 1 + apply(kappa,2,sum)[-R]
    revcs_kappa = kappa
    for (ridx in 1:nrow(kappa)) revcs_kappa[ridx,] = revcumsum(revcs_kappa[ridx,])
    gam2s = alp.ratio + apply(revcs_kappa[, -1],2,sum)
    logstick = digamma(gam1s) - digamma(gam1s+gam2s)
    log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
    qstick = c(logstick, 0) + c(0, log1mstick)
    
    # alpha
    balptl = balp - sum(digamma(gam2s)-digamma(gam1s+gam2s))
    alp.ratio = aalptl / balptl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(ti*sigu.q)
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2)
    sig.ratio = asigtl / bsigtl
    
    # Convergence
    bool1 = mean((mubeta.q.old - mubeta.q)^2) < eps
    bool2 = mean((muu.q.old - muu.q)^2) < eps
    if (bool1 & bool2) break
  }
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl,
    alamtl=alamtl, blamtl=blamtl
  ))
}


# BSA: Exponential mixture inefficiency -----------------------------------

exmixture_bsa = function(y, x, w, Z, J, R=10, productivity=TRUE, prior=NULL, maxiter=500, eps=1e-4)
{
  library(matrixStats)
  sgn = (-1)^productivity
  N = ncol(Z)
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  ti = apply(Z, 2, sum)
  # zero-one standardization of x
  minx = min(x)
  maxx = max(x)
  if(minx < 0 | maxx > 1) x = (x-minx)/(maxx-minx)
  
  # Hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(rep(100, D+1))
    w0 = 2
    mupsi.q.start = 1
    alam = blam = aalp = balp = atau = btau = asig = bsig = 1e-3
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    alam = prior$alam
    blam = prior$blam
    aalp = prior$aalp
    balp = prior$balp
    asig = prior$asig
    bsig = prior$bsig
    atau = prior$atau
    btau = prior$btau
    w0 = prior$w0
    mupsi.q.start = prior$mupsi.q.start
  }
  sb0i = solve(sigbeta.0)
  sbimb0 = sb0i %*% mubeta.0
  aalptl = aalp + R - 1
  asigtl = asig + 0.5*(sum(ti))
  
  # Initialize parameters
  # Level1
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  sig.ratio = asig / bsig
  # Level2
  kappa = matrix(1,N,R) / R
  kappa_r = apply(kappa,2,sum)
  lam.ratio = rep(alam/blam, R)
  tau.ratio = atau / btau
  # Level3 or above
  gam1s = rep(0.01, R-1)
  gam2s = rep(0.01, R-1)
  alp.ratio = aalp / balp
  logstick = digamma(gam1s) - digamma(gam1s+gam2s)
  log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
  qstick = c(logstick, 0) + c(0, log1mstick)
  # Nonparametric
  mutheta.q = rep(0,J)
  mupsi.q = mupsi.q.start
  mu2psi.q = mupsi.q^2
  sig2psi.q = mu2psi.q/100
  sigpsi.q = sqrt(sig2psi.q)
  Jfull = J
  vphifull = sqrt(2)*cos(outer(x,pi*(1:J)))
  vphitvphifull = t(vphifull)%*%vphifull
  vphi = vphifull[,1:J]
  vphitvphi = vphitvphifull[1:J,1:J]
  
  # Update as
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
    # Determine J and corresponding values
    J = min(floor(-15/(-0.5*mupsi.q)),Jfull)    
    mutheta.q.old = mutheta.q[1:J]
    asigtl = asig + 0.5*(nrow(w)+J)
    atautl = atau + 0.5*J
    bindices = (1:J)
    bindices2= bindices^2
    mfactor = -J*(J+1)/(4*w0)
    
    # Set up spectral design matrix
    vphi = vphifull[,1:J]
    vphitvphi = vphitvphifull[1:J,1:J]
    
    # Update variational distribution parameters for 
    # theta
    term1 = exp(sig2psi.q*bindices2/2 + mupsi.q*bindices)
    term2 = exp(sig2psi.q*bindices2/2 - mupsi.q*bindices)
    sigpsi.q = sqrt(sig2psi.q)
    Qvec = term1 * (1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec = Qvec + term2 * (1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices))
    DQvec = diag(Qvec)
    
    sigtheta.q = sig.ratio * (tau.ratio*DQvec+vphitvphi)
    sigtheta.q = solve(sigtheta.q)
    mutheta.q = sig.ratio * sigtheta.q %*% t(vphi) %*% (y-W%*%mubeta.q-sgn*Z%*%muu.q)
    mutheta.q = drop(mutheta.q)
    
    # beta
    sigbeta.q = solve(sb0i+sig.ratio*WtW)
    mubeta.q = sigbeta.q %*% (sbimb0 + sig.ratio*t(W)%*%(y-sgn*Z%*%muu.q-vphi%*%mutheta.q))
    mubeta.q = drop(mubeta.q)
    
    # inefficiency
    sigi2 = 1/(sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2 * (-apply(t(kappa)*lam.ratio,2,sum) + sgn*sig.ratio*drop(t(Z)%*%(y-W%*%mubeta.q-vphi%*%mutheta.q)))
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE)-pnorm(musigrat, log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2 * (1-fres*(musigrat+fres))
    
    # lambda_r
    alamtl = alam + kappa_r
    blamtl = blam + apply(kappa*muu.q, 2, sum)
    lam.ratio = alamtl / blamtl
    
    # kappa
    S = outer(muu.q, lam.ratio, times)
    S = t(qstick - log(blamtl) + digamma(alamtl) + t(S))
    kappa = exp(kappa - rowLogSumExps(kappa))
    kappa_r = apply(kappa, 2, sum)
    
    # stick length(Vr)
    gam1s = 1 + apply(kappa,2,sum)[-R]
    revcs_kappa = kappa
    for (ridx in 1:nrow(kappa)) revcs_kappa[ridx,] = revcumsum(revcs_kappa[ridx,])
    gam2s = alp.ratio + apply(revcs_kappa[, -1],2,sum)
    logstick = digamma(gam1s) - digamma(gam1s+gam2s)
    log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
    qstick = c(logstick, 0) + c(0, log1mstick)
    
    # alpha
    balptl = balp - sum(digamma(gam2s)-digamma(gam1s+gam2s))
    alp.ratio = aalptl / balptl
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q-vphi%*%mutheta.q)^2)
    trterm1 = sum(diag(WtW %*% sigbeta.q))
    trterm2 = sum(ti * sigu.q)
    trterm3 = sum(diag(vphitvphi %*% sigtheta.q))
    trterm4 = sum(diag((outer(mutheta.q, mutheta.q)+sigtheta.q)%*%DQvec))
    bsigtl = bsig + 0.5*(ssterm+trterm1+trterm2+trterm3+tau.ratio*trterm4)
    sig.ratio = asigtl / bsigtl
    
    # tau
    btautl = btau + 0.5*(sig.ratio*trterm4)
    tau.ratio = atautl / btautl
    
    # psi(NCVMP)
    sig3psi.q = sigpsi.q^3
    sig4psi.q = sigpsi.q^4
    mu2psi.q = mupsi.q^2
    const1 = sqrt(2/pi)
    const2 = exp(-mu2psi.q/(2*sig2psi.q))
    pd1sig = -w0*((1/(2*sig2psi.q)+sigpsi.q*mu2psi.q/(2*sig4psi.q))*const1*const2-mu2psi.q/sig3psi.q*dnorm(-mupsi.q/sigpsi.q))
    pd1mu = -w0*(-mupsi.q/sigpsi.q*const1*const2+(1-2*pnorm(-mupsi.q/sigpsi.q))+2*mupsi.q/sigpsi.q*dnorm(-mupsi.q/sigpsi.q))
    
    term1 = exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2 = exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    term3 = mupsi.q/(2*sig3psi.q)-bindices/(2*sig2psi.q)
    term4 = mupsi.q/(2*sig3psi.q)+bindices/(2*sig2psi.q)
    term5 = dnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term6 = dnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term7 = 1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term8 = 1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)
    DQmu = term1*term5*1/sigpsi.q+bindices*term1*term7-term2*term6*1/sigpsi.q-bindices*term2*term8    
    DQsig =  -term1*term5*term3+bindices2/2*term1*term7+term2*term6*term4+bindices2/2*term2*term8    
    pd2sig =  mfactor*pd1sig-0.5*sig.ratio*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*DQsig)
    pd2mu =  mfactor*pd1mu-0.5*sig.ratio*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*DQmu)
    sig2psi.q.old = sig2psi.q
    sig2psi.q =  -0.5/(pd1sig+pd2sig)
    mupsi.q.old = mupsi.q
    
    # Convergence
    bool1 = mean((mubeta.q.old - mubeta.q)^2) < eps
    bool2 = mean((muu.q.old - muu.q)^2) < eps
    if (bool1 & bool2) break
  }
  post_curve=drop(vphi%*%mutheta.q)
  curve_var=drop(vphi^2%*%diag(sigtheta.q))
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl,
    alamtl=alamtl, blamtl=blamtl,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve
  ))
}
