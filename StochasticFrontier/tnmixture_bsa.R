tnmixture_bsa = function(y, x, w, Z, J, R=10, productivity=TRUE, prior=NULL, maxiter=500, eps=1e-4, acc=1000)
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
    mu0 = 0
    sig02 = 100
    w0 = 2
    mupsi.q.start = 1
    alam = blam = aalp = balp = asig = bsig = atau = btau = 1e-3
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    mu0 = prior$mu0
    sig02 = prior$sig02
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
  ldetsb0 = determinant(sigbeta.0, logarithm=TRUE)$modulus
  aalptl = aalp + R - 1
  asigtl = asig + 0.5*(sum(ti))
  
  # Initialize parameters
  # Level1
  mubeta.q = rep(0, D+1)
  muu.q = rep(0, N)
  sig.ratio = asig / bsig
  # Level2
  b_lam2 = rep(1, R)
  b_lam1 = rep(1, R)
  b_mu1 = rep(0, R)
  kappa = matrix(1,N,R) / R
  kappa_r = apply(kappa,2,sum)
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
  
  # grids to be used in griddy gibbs routine
  real = seq(-8, 8, length.out=acc) # Interval wider than this would cause numerical overflow
  if (sum(real == 0) == 1) real = real[-which(real == 0)]
  neg = real < 0; pos = real > 0
  positive = seq(eps, 50, length.out=acc)
  if (sum(positive == 1) == 1) positive = positive[-which(positive == 1)]
  less = positive < 1; greater = positive > 1
  
  # functions to be used
  times = function(x, y) x*y
  revcumsum = function(x) {rev(cumsum(rev(x)))}
  
  # ELBO container
  elbo = rep(0, maxiter+1)
  elbo[1] = -Inf
  
  # Iterate
  for (iter in 1:maxiter)
  {
    # beta test(see if algorithm converges to right posterior)
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
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = sigbeta.q %*% (sbimb0 + sig.ratio*t(W)%*%(y-sgn*Z%*%muu.q-vphi%*%mutheta.q))
    mubeta.q = drop(mubeta.q)
    
    # inefficiency
    sigi2 = 1 / (apply(t(kappa)*b_lam2,2,sum) + sig.ratio*ti)
    sigi = sqrt(sigi2)
    mui = sigi2 * (apply(t(kappa)*b_lam1*b_mu1,2,sum) + sgn*sig.ratio*drop(t(Z)%*%(y-W%*%mubeta.q-vphi%*%mutheta.q)))
    musigrat = mui / sigi
    fres = exp(dnorm(musigrat, log=TRUE)-pnorm(musigrat, log.p=TRUE))
    muu.q = mui + sigi*fres
    sigu.q = sigi2 * (1-fres*(musigrat+fres))
    
    # griddy Gibbs routine
    # mu
    # Hg construction
    Hg = outer(real, mu0/sig02+b_lam1*apply(kappa*muu.q,2,sum), times)
    Hg = Hg - outer(real^2, 0.5*(1/sig02+kappa_r), times)
    Hg = Hg - outer(pnorm(real, log.p=TRUE), kappa_r, times)
    # Moment calculation
    xi_mu = colLogSumExps(Hg)
    b_mu1 = exp(colLogSumExps(Hg[pos,]+log(real[pos])) - xi_mu) - exp(colLogSumExps(Hg[neg,]+log(-real[neg])) - xi_mu)
    b_mu2 = b_mu1^2
    b_lpmu = -exp(colLogSumExps(Hg+log(-log(pnorm(real)))) - xi_mu)
    # ilambda2
    Hg = outer(log(positive), alam+0.5*kappa_r-1, times)
    Hg = Hg - outer(positive, blam+0.5*apply(kappa*(muu.q^2+sigu.q),2,sum), times)
    Hg = Hg + outer(sqrt(positive), b_mu1*apply(kappa*muu.q,2,sum), times)
    # Moment calcuation
    xi_lam2 = colLogSumExps(Hg)
    b_lam2 = exp(colLogSumExps(Hg+log(positive)) - xi_lam2)
    b_lam1 = sqrt(b_lam2)
    b_llam2 = exp(colLogSumExps(Hg[greater,]+log(log(positive[greater]))) - xi_lam2) - exp(colLogSumExps(Hg[less,]+log(-log(positive[less]))) - xi_lam2)
    
    # kappa
    S = -0.5*outer(muu.q^2+sigu.q, b_lam2, times)
    S = S + outer(muu.q, b_lam1*b_mu1, times)
    S = t(t(S) + 0.5*(b_llam2-b_mu2) - b_lpmu + qstick)
    kappa = exp(S - rowLogSumExps(S))
    kappa_r = apply(kappa,2,sum)
    
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
    
    # # ELBO
    # term1 = -0.5*nrow(W)*(log(2*pi)-digamma(bsigtl)+log(asigtl))-0.5*sig.ratio*(ssterm+trterm1+trterm2)
    # term2 = -0.5*(ldetsb0+sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0))+sum(diag(sb0i%*%sigbeta.q)))
    # term2 = term2 + 0.5*(determinant(sigbeta.q,logarithm=TRUE)$modulus+sum(diag(sigbeta.q%*%sb0i)))
    # term3 = -0.5*outer(b_lam2,muu.q^2+sigu.q,times)+outer(b_lam1*b_mu1,muu.q,times)+0.5*(b_lam2-b_mu2)-b_lpmu
    # term3 = sum(diag(kappa%*%term3)) + 0.5*(sum(sigu.q)+N*(log(2*pi)+1)+sum(exp(dnorm(muu.q/sqrt(sigu.q),log=TRUE)-pnorm(muu.q/sqrt(sigu.q),log.p=TRUE))))
    # term4 = -0.5*R*(log(2*pi)+log(sig02))-0.5*sum((b_mu1-mu0)^2)/sig02
    # term4 = term4 + sum(xi_mu + 0.5*(1/sig02+kappa_r)*b_mu2 - (mu0/sig02+b_lam1*apply(kappa*muu.q,2,sum))*b_mu1 + kappa_r*b_lpmu)
    # term5 = -lgamma(alam)+alam*log(blam)+(alam-1)*sum(b_llam2)-blam*sum(b_lam2)
    # term5 = term5 + sum(xi_lam2-(alam+0.5*kappa_r-1)*b_llam2+(blam+0.5*apply(kappa*(muu.q^2+sigu.q),2,sum))*b_lam2-b_mu1*apply(kappa*muu.q,2,sum)*b_lam1)
    # term6 = sum(t(kappa)*(qstick-log(t(kappa))))
    # term8 = (digamma(aalptl)-log(balptl))*(aalp-aalptl)-aalptl/balptl*(balp-balptl)-(lgamma(aalp)-lgamma(aalptl))+aalp*log(balp)-aalptl*log(balptl)
    # term9 = (digamma(asigtl)-log(bsigtl))*(asig-asigtl)-asigtl/bsigtl*(bsig-bsigtl)-(lgamma(asig)-lgamma(asigtl))+asig*log(bsig)-asigtl*log(bsigtl)
    # elbo[iter+1] = term1+term2+term3+term4+term5+term6+term8+term9
    # if (abs(elbo[iter+1]-elbo[iter]) < eps) break
    
    # Convergence
    bool1 = mean((mubeta.q.old - mubeta.q)^2) < eps
    bool2 = mean((muu.q.old - muu.q)^2) < eps
    if (bool1 & bool2) break
  }
  # elbo=elbo[-1]
  # elbo=elbo[-which(elbo==0)]
  # browser()
  post_curve=drop(vphi%*%mutheta.q)
  curve_var=drop(vphi^2%*%diag(sigtheta.q))
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    asigtl=asigtl, bsigtl=bsigtl,
    mur=b_mu1, lam2r=b_lam2,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve
  ))
}