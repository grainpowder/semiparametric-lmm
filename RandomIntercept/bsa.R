bsa = function(y,x,w,Z,J,prior=NULL,maxiter=500,tol=1e-5){
  # Pre-calculate objects to be used frequently
  N = nrow(w)
  D = ncol(w)
  n = ncol(Z)
  W = cbind(1,w)
  C = cbind(W,Z)
  bidx = 1:(D+1)
  uidx = (D+2):(D+n+1)
  CtC = crossprod(C)
  # zero-one standardization of x
  minx = min(x)
  maxx = max(x)
  if(minx < 0 | maxx > 1) x = (x-minx)/(maxx-minx)
  
  # Initialize parameters
  if(is.null(prior)){
    rsig.0 = rtau.0 = rlam.0 = 0.1
    ssig.0 = stau.0 = slam.0 = 1
    sigbeta2 = 1e3
    w0 = 2
    mupsi.q.start = 1}
  else{
    rsig.0 = prior$rsig.0
    rtau.0 = prior$rtau.0
    rlam.0 = prior$rlam.0
    ssig.0 = prior$ssig.0
    stau.0 = prior$stau.0
    slam.0 = prior$slam.0
    sigbeta2 = prior$sigbeta2
    w0 = prior$w0
    mupsi.q.start = prior$mupsi.q.start}
  # parametric
  sig.ratio = rsig.0/ssig.0
  rlam.q = rsig.0+n
  lam.ratio = rlam.0/slam.0
  tau.ratio = rtau.0/stau.0
  muphi.q = rep(0,D+n+1)
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
  
  # Iteration
  for(iter in 1:maxiter){
    muphi.q.old = muphi.q
    
    # Determine J and corresponding values
    J = min(floor(-15/(-0.5*mupsi.q)),Jfull)    
    mutheta.q.old = mutheta.q[1:J]
    rsig.q = rsig.0+N+J
    rtau.q = rtau.0+J
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
    mutheta.q = sig.ratio*sigtheta.q%*%t(vphi)%*%(y-C%*%muphi.q)
    mutheta.q = drop(mutheta.q)
    
    # phi
    sigphi.q = sig.ratio*CtC+diag(rep(c(1/sigbeta2,lam.ratio),c(D+1,n)))
    sigphi.q = solve(sigphi.q)
    muphi.q = sig.ratio*sigphi.q%*%t(C)%*%(y-vphi%*%mutheta.q)
    muphi.q = drop(muphi.q)
    
    # lambda
    slam.q = slam.0+sum((muphi.q^2+diag(sigphi.q))[uidx])
    lam.ratio = rlam.q/slam.q
    
    # sigma
    ssterm = sum((y-C%*%muphi.q-vphi%*%mutheta.q)^2)
    trterm1 = sum(diag((outer(mutheta.q, mutheta.q)+sigtheta.q)%*%DQvec))
    trterm2 = sum(diag(CtC%*%sigphi.q))
    trterm3 = sum(diag(vphitvphi%*%sigtheta.q))
    ssig.q = ssig.0+ssterm+tau.ratio*trterm1+trterm2+trterm3
    sig.ratio = rsig.q/ssig.q
    
    # tau
    stau.q = stau.0+sig.ratio*trterm1
    tau.ratio = rtau.q/stau.q
    
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
    DQsig =  -term1*term5*term3+bindices2/2*term1*term7+term2*term6*term4+bindices2/2*term2*term8    
    pd2sig =  mfactor*pd1sig-0.5*sig.ratio*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*DQsig)
    pd2mu =  mfactor*pd1mu-0.5*sig.ratio*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*DQmu)
    sig2psi.q.old = sig2psi.q
    sig2psi.q =  -0.5/(pd1sig+pd2sig)
    mupsi.q.old = mupsi.q
    
    # convergence
    if(mean((muphi.q.old-muphi.q)^2)<tol & mean((mutheta.q.old-mutheta.q)^2)<tol) break}
  
  post_curve=drop(vphi%*%mutheta.q)
  curve_var=drop(vphi^2%*%diag(sigtheta.q))
  return(list(
    mubeta.q=muphi.q[bidx],sigbeta.q=sigphi.q[bidx,bidx],
    muu.q=muphi.q[uidx],sigu.q=sigphi.q[uidx,uidx],
    mutheta.q=mutheta.q,sigtheta.q=sigtheta.q,
    rsig.q=rsig.q,ssig.q=ssig.q,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve))}