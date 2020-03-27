bsar_vb = function(y, x, w=NULL, J=10, prior=NULL, maxiter=500, tol=1e-6)
{
  # Code to estimate parameters of simple semiparametric regression using Variational Inference
  # Uses Gaussian process prior to estimate nonlinear relationship between x and y
  # 
  # y : response vector
  # x : explanatory variable(s) whose nonlinear relationship has to be estimated
  # w : explanatory variable(s) whose linear relationship has to be estimated
  N = length(y)
  if (is.null(w)) {D = 0} else {D = ncol(w)}
  W = cbind(matrix(1, nrow=N), w)
  WtW = crossprod(W)
  
  # Hyperparameters
  if (is.null(prior))
  {
    mubeta.0 = rep(0, D+1)
    if (D == 0) {sigbeta.0 = matrix(100,1,1)} else {sigbeta.0 = diag(rep(100, D+1))}
    asig = bsig = atau = btau = 1e-3
    w0 = 2
    mupsi.q.start = 1
  }
  else
  {
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
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
  
  # zero-one standardization of x
  minx = min(x)
  maxx = max(x)
  if(minx < 0 | maxx > 1) x = (x-minx)/(maxx-minx)
  
  # Initialize variational parameters
  # parametric
  sig.ratio = asig/bsig
  tau.ratio = atau/btau
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
    J.old = J
    J = min(floor(-15/(-0.5*mupsi.q)),Jfull)
    mutheta.q.old = mutheta.q[1:J]
    asigtl = asig+(N+J)*0.5
    atautl = atau+(J/2)
    bindices = (1:J)
    bindices2= bindices^2
    mfactor = -J*(J+1)/(4*w0)
    
    # Set up spectral design matrix
    vphi = vphifull[,1:J]
    vphitvphi = vphitvphifull[1:J,1:J]
    
    # Update variational distribution parameters for 
    # beta
    sigbeta.q = solve(sb0i+sig.ratio*WtW)
    mubeta.q = sigbeta.q%*%(sbimb0+sig.ratio*t(W)%*%(y-vphi%*%mutheta.q))
    mubeta.q = drop(mubeta.q)
    
    # theta
    term1 = exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2 = exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    sigpsi.q = sqrt(sig2psi.q)
    Qvec = term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec = Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices))
    DQvec = diag(Qvec)
    
    sigtheta.q = sig.ratio*(tau.ratio*DQvec+vphitvphi)
    sigtheta.q = solve(sigtheta.q)
    mutheta.q = sig.ratio*sigtheta.q%*%t(vphi)%*%(y-W%*%mubeta.q)
    mutheta.q = drop(mutheta.q)
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-vphi%*%mutheta.q)^2)
    trterm1 = sum(diag(WtW%*%sigbeta.q))
    trterm2 = sum(diag(vphitvphi%*%sigtheta.q))
    trterm3 = sum(diag((outer(mutheta.q,mutheta.q)+sigtheta.q)%*%DQvec))
    bsigtl = bsig+0.5*(ssterm+trterm1+trterm2+tau.ratio*trterm3)
    sig.ratio = asigtl/bsigtl
    
    # tau
    btautl = btau+0.5*sig.ratio*trterm3
    tau.ratio = atautl/btautl
    lntau = log(btautl) - digamma(atautl)
    lnsig = log(bsigtl) - digamma(asigtl)
    
    # ELBO calculation(before NCVMP)
    lbnew = 0
    sigpsi.q = sqrt(sig2psi.q)
    mu2psi.q = mupsi.q^2
    const1 = sqrt(2/pi)
    # likelihood
    lbnew = lbnew - (N/2*log(2*pi)) - (N/2)*lnsig - sig.ratio*(ssterm+trterm1+trterm2)/2
    # beta
    lbnew = lbnew - (ldetsb0 - determinant(sigbeta.q,logarithm=TRUE)$modulus)/2 - (sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0))+sum(diag(solve(sigbeta.0,sigbeta.q)))-(D+1))/2
    # theta
    lbnew = lbnew - J*lnsig/2 - J*lntau/2 - sig.ratio*tau.ratio*trterm3/2 + determinant(sigtheta.q,logarithm=TRUE)$modulus/2 + J/2
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
    lb=lb, mutheta.q=mutheta.q, sigtheta.q=sigtheta.q,
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    asigtl=asigtl, bsigtl=bsigtl,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve
  ))
}
