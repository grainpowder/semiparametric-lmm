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
  
  # functions to be used
  times = function(x, y) x*y
  revcumsum = function(x) {rev(cumsum(rev(x)))}
  
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