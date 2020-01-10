vbbsasf_exp = function(y,x,w,Z,J,productivity=TRUE,prior=NULL,tol=1e-4,maxiter=500){
  # Pre-calculate objects to be used frequently
  sgn = (-1)^productivity
  n = ncol(Z)
  N = nrow(w)
  D = ncol(w)
  W = cbind(1,w)
  WtW = crossprod(W)
  ZtZ = crossprod(Z)
  ni = apply(Z,2,sum)
  # zero-one standardization of x
  minx = min(x)
  maxx = max(x)
  if(minx < 0 | maxx > 1) x = (x-minx)/(maxx-minx)
  
  # Initialize parameters
  if(is.null(prior)){
    rsig.0 = rtau.0 = rlam.0 = 0.1
    ssig.0 = stau.0 = slam.0 = 10
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
  sb2diag = diag(1/sigbeta2,D+1)
  rlam.q = rsig.0+2*n
  sig.ratio = rsig.0/ssig.0
  lam.ratio = rlam.0/slam.0
  tau.ratio = rtau.0/stau.0
  mubeta.q = rep(0,D+1)
  muu.q = rep(0,n)
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
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
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
    mutheta.q = sig.ratio*sigtheta.q%*%t(vphi)%*%(y-W%*%mubeta.q-sgn*Z%*%muu.q)
    mutheta.q = drop(mutheta.q)
    
    # beta
    sigbeta.q = sig.ratio*WtW+sb2diag
    sigbeta.q = solve(sigbeta.q)
    mubeta.q = sig.ratio*sigbeta.q%*%t(W)%*%(y-sgn*Z%*%muu.q-vphi%*%mutheta.q)
    mubeta.q = drop(mubeta.q)
    
    # u
    sigi2 = 1/(ni*sig.ratio)
    sigi = sqrt(sigi2)
    mui = drop(-lam.ratio+sgn*sig.ratio*(t(Z)%*%(y-W%*%mubeta.q-vphi%*%mutheta.q)))
    mui = sigi2*mui
    musigrat = mui/sigi
    fres = exp(dnorm(musigrat,log=TRUE)-pnorm(musigrat,log.p=TRUE))
    muu.q = mui+sigi*fres
    sigu.q = sigi2*(1-fres*(musigrat+fres))
    
    # lambda
    slam.q = slam.0+2*sum(muu.q)
    lam.ratio = rlam.q/slam.q
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-sgn*Z%*%muu.q-vphi%*%mutheta.q)^2)
    trterm1 = sum(diag((outer(mutheta.q,mutheta.q)+sigtheta.q)%*%DQvec))
    trterm2 = sum(diag(WtW%*%sigbeta.q))
    trterm3 = sum(ZtZ%*%sigu.q)
    trterm4 = sum(diag(vphitvphi%*%sigtheta.q))
    ssig.q = ssig.0+ssterm+tau.ratio*trterm1+trterm2+trterm3+trterm4
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
    mse1 = mean((mubeta.q.old-mubeta.q)^2)
    mse2 = mean((muu.q.old-muu.q)^2)
    mse3 = mean((mutheta.q.old-mutheta.q)^2)
    if(mse1<tol & mse2<tol & mse3<tol) break}
  
  post_curve=drop(vphi%*%mutheta.q)
  curve_var=drop((vphi^2)%*%diag(sigtheta.q))
  return(list(
    mubeta.q=mubeta.q,sigbeta.q=sigbeta.q,
    muu.q=muu.q,sigu.q=sigu.q,
    mutheta.q=mutheta.q,sigtheta.q=sigtheta.q,
    rsig.q=rsig.q,ssig.q=ssig.q,
    rlam.q=rlam.q,slam.q=slam.q,
    rtau.q=rtau.q,stau.q=stau.q,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve))}

# # simulation setting
# source("../misc/make_Z.R")
# set.seed(2)
# nobs = 50
# nrep = 4; D = 10
# Z = make_Z(rep(nrep, nobs))
# betaT = rnorm(D+1); p = length(betaT)
# sigma = 0.5
# lambda = 2
# f = function(x) -(x-1.5)^2
# # data generation
# x = 3*runif(nobs*nrep); ord = order(x)
# w = matrix(rnorm(length(x)*(p-1)), ncol=(p-1))
# u = rexp(nobs,lambda)
# y = drop(cbind(1,w)%*%betaT + f(x) - rep(u, each=nrep) + rnorm(length(x), sd=sigma))
# result = vbbsasf_exp(y,x,w,Z,23)
# res = y-(cbind(1,w)%*%result$mubeta.q-Z%*%result$muu.q)
# par(mar=c(4,4.5,2,1),mfrow=c(1,2))
# plot(betaT, result$mubeta.q, main="Fixed effects",
#      ylab=expression(hat(beta)), xlab=expression(beta))
# lines(-5:5,-5:5)
# plot(u, result$muu.q, main="Random effects",
#      ylab=expression(hat(u)))
# lines(-1:10,-1:10)
# par(mar=c(4,2,2,1),mfrow=c(1,1))
# plot(x,res, main="Fitted mean curve and 95% credible region", ylab="")
# lines(x[ord],result$post_curve[ord],lwd=3,col=3)
# lines(x[ord],result$post_upper[ord],lwd=2,lty=2)
# lines(x[ord],result$post_lower[ord],lwd=2,lty=2)
