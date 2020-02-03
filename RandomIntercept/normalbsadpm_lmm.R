normalbsadpm_lmm = function(y,x,w,Z,J,R,prior=NULL,maxiter=500,tol=1e-5)
{
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  N = ncol(Z)
  T = nrow(Z)
  ti = apply(Z, 2, sum)
  # zero-one standardization of x
  minx = min(x)
  maxx = max(x)
  if(minx < 0 | maxx > 1) x = (x-minx)/(maxx-minx)
  # prior parameters
  # Note : be sure to make base distribution as flat as possible
  if (is.null(prior))
  {
    asig = bsig = 0.1
    atau = btau = 0.1
    aalp = balp = 0.1
    alam = blam = 0.01
    mu0 = 0
    sig02 = 100
    mubeta.0 = rep(0, D+1)
    sigbeta.0 = diag(100, D+1)
    w0 = 2
    mupsi.q.start = 1
  }
  else
  {
    asig = prior$asig
    bsig = prior$bsig
    aalp = prior$aalp
    balp = prior$balp
    atau = prior$atau
    btau = prior$btau
    alam = prior$alam
    blam = prior$blam
    mu0 = prior$mu0
    sig02 = prior$sig02
    mubeta.0 = prior$mubeta.0
    sigbeta.0 = prior$sigbeta.0
    w0 = prior$w0
    mupsi.q.start = prior$mupsi.q.start
  }
  # Define containers
  tau.ratio = atau / btau
  alp.ratio = aalp / balp
  sig.ratio = asig / bsig
  alamtls = rep(alam, R)
  blamtls = rep(blam, R)
  lam.ratio = alamtls / blamtls
  muu.q = rep(0, N)
  sigu.q = rep(sig02, N)
  kappa = matrix(runif(N*R), N, R)
  kappa = kappa / apply(kappa, 1, sum)
  mutls = rep(0, R)
  sigtls = rep(sig02, R)
  gam1s = gam2s = rep(1, R)
  mubeta.q = rep(0, D+1)
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
  # Define predetermined values
  sb0i = solve(sigbeta.0)
  sb0imb0 = sb0i %*% mubeta.0
  asigtl = asig + (N*T/2)
  aalptl = aalp + R - 1
  # Function to use during iteration
  numreplace = function(x, tol=1e-10) {x[x<tol] = tol; x} # function to replace value(s) less than 1e-10 to 1e-10
  revcumsum = function(x) {rev(cumsum(rev(x)))}
  # Begin iteration
  for (iter in 1:maxiter)
  {
    mubeta.q.old = mubeta.q
    muu.q.old = muu.q
    
    # Determine J and corresponding values
    J = min(floor(-15/(-0.5*mupsi.q)),Jfull)    
    mutheta.q.old = mutheta.q[1:J]
    asigtl = asig+0.5*(N*T+J)
    atautl = atau+0.5*J
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
    mutheta.q = sig.ratio*sigtheta.q%*%t(vphi)%*%(y-W%*%mubeta.q-Z%*%muu.q)
    mutheta.q = drop(mutheta.q)
    
    # beta
    sigbeta.q = solve(sb0i + sig.ratio*WtW)
    mubeta.q = sigbeta.q %*% (sb0imb0+sig.ratio*t(W)%*%(y-Z%*%muu.q-vphi%*%mutheta.q))
    mubeta.q = drop(mubeta.q)
    
    # random intercept
    kappa_lam = t(kappa) * lam.ratio
    sigu.q = 1 / (apply(kappa_lam,2,sum) + ti*sig.ratio)
    muu.q = apply(kappa_lam*mutls,2,sum) + sig.ratio*t(Z)%*%(y-W%*%mubeta.q-vphi%*%mutheta.q)
    muu.q = drop(sigu.q * muu.q)
    
    # kappa(prior...)
    logstick = digamma(gam1s) - digamma(gam1s+gam2s)
    log1mstick = cumsum(digamma(gam2s) - digamma(gam1s+gam2s))
    qstick = c(logstick, 0) + c(0, log1mstick)
    S = (t(matrix(muu.q, N, R))-mutls)^2 + sigtls
    S = t(t(S)+sigu.q)*lam.ratio + digamma(alamtls) - log(blamtls)
    S = t(qstick + 0.5*S)
    S = exp(S - apply(S,1,max))
    for (cidx in 1:ncol(S)) S[, cidx] = numreplace(S[, cidx])
    kappa = S / apply(S, 1, sum)
    
    # stick length
    gam1s = 1 + apply(kappa, 2, sum)[-R]
    revcs_kappa = kappa
    for (ridx in 1:nrow(kappa)) revcs_kappa[ridx,] = revcumsum(revcs_kappa[ridx,])
    gam2s = alp.ratio + apply(revcs_kappa[, -1], 2, sum)
    
    # alpha
    balptl = balp - sum(digamma(gam2s)-digamma(gam1s+gam2s))
    alp.ratio = aalptl / balptl
    
    # parameters: mean
    sigtls = 1/(1/sig02 + apply(t(kappa)*lam.ratio, 1, sum))
    mutls = mu0/sig02 + apply(t(kappa*muu.q)*lam.ratio, 1, sum)
    
    # parameters: variance
    alamtls = alam + 0.5*apply(kappa, 1, sum)
    SS = (t(matrix(muu.q, N, R))-mutls)^2 + sigtls
    SS = t(SS) + sigu.q
    SS = kappa*SS
    blamtls = blam + 0.5*apply(SS, 1, sum)
    lam.ratio = alamtls / blamtls
    
    # sigma
    ssterm = sum((y-W%*%mubeta.q-Z%*%muu.q-vphi%*%mutheta.q)^2)
    trterm1 = sum(diag(WtW%*%mubeta.q))
    trterm2 = sum(ti*sigu.q)
    trterm3 = sum(diag(vphitvphi%*%sigtheta.q))
    trterm4 = sum(diag((outer(mutheta.q, mutheta.q)+sigtheta.q)%*%DQvec))
    bsigtl = ssterm + 0.5*(trterm1 + trterm2 + trterm3 + trterm4*tau.ratio)
    sig.ratio = asigtl / bsigtl
    
    # tau
    btautl = btau + 0.5*sig.ratio*trterm4
    tau.ratio = atautl / btautl
    
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
    
    # Convergence
    convergence1 = sum((mubeta.q - mubeta.q.old)^2) < tol
    convergence2 = sum((muu.q - muu.q.old)^2) < tol
    if (convergence1 & convergence2) break
  }
  
  post_curve=drop(vphi%*%mutheta.q)
  curve_var=drop(vphi^2%*%diag(sigtheta.q))
  return(list(
    mubeta.q=mubeta.q, sigbeta.q=sigbeta.q,
    muu.q=muu.q, sigu.q=sigu.q,
    post_lower=qnorm(0.025,post_curve,sqrt(curve_var)),
    post_upper=qnorm(0.975,post_curve,sqrt(curve_var)),
    post_curve=post_curve,
    mutls=mutls, sigtls=sigtls,
    kappa=kappa
  ))
}

# # simulation
# source("../misc/make_Z.R")
# N = 50; T = 4; D = 10 # D = 5 causes error. Seems like moderate unstability is present.
# Z = make_Z(rep(T, N))
# set.seed(10)
# beta = rnorm(D+1)
# w = matrix(rnorm(N*T*D), N*T, D)
# f = function(x) 3*x*sin(pi*x)
# x = 3*runif(N*T); ord = order(x)
# assigner = runif(N)
# mu1 = 2; mu2 = -2
# u = rnorm(N, mu1)
# u[assigner > 0.5] = rnorm(sum(assigner > 0.5), mu2)
# y = cbind(1, w)%*%beta + Z%*%u + f(x) + rnorm(nrow(Z))
# result = normalbsadpm_lmm(y,x,w,Z,23,10)
# par(mfrow=c(1,2))
# plot(beta[-1], result$mubeta.q[-1])
# lines(-5:5,-5:5)
# plot(u, result$muu.q)
# lines(-5:5,-5:5)
# 
# res = y-(cbind(1,w)%*%result$mubeta.q-Z%*%result$muu.q)
# par(mar=c(4,2,2,1),mfrow=c(1,1))
# plot(x,res, main="Fitted mean curve and 95% credible region", ylab="")
# lines(x[ord],result$post_curve[ord],lwd=3,col=3)
# lines(x[ord],result$post_upper[ord],lwd=2,lty=2)
# lines(x[ord],result$post_lower[ord],lwd=2,lty=2)