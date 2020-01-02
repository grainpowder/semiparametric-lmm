wandint = function(f,h,a,b,k,init){
  # Requires numDeriv package while calculating second derivative of h
  # Calculate integral (*) numerically by the method introduced in Appendix B of Wand(2011).
  # 1. Note that h(x) has to be concave function (i.e. h''(x) < 0 for all a < x <b)
  # 2. As this function uses built-in `optim` function to differentiate h, init s.t. a < init < b is required.
  # 3. To prevent numerical underflow, ln{f(x)+k}+h(x) is integrated alternatively.
  #
  # (*) E_X[f(X)] = \int_a^b f(x)exp{h(x)}dx
  #               = \int_a^b {f(x)+k-k}exp{h(x)}dx
  #               = \int_a^b exp[ln{f(x)+k}+h(x)]dx - k\int_a^b exp{h(x)}dx
  #               = I1                              - k*I2
  #
  integrand = function(x) log(f(x)+k)+h(x)
  mu0 = optim(init,integrand,method="L-BFGS-B",lower=a,upper=b)$par
  sig = sqrt(-2/numDeriv::hessian(integrand,mu0)[1])
  h1 = function(x) exp(integrand(mu0+sig*x)-integrand(mu0))
  lnI1 = integrand(mu0)+log(sig)+log(integrate(h1,lower=(a-mu0)/sig,upper=(b-mu0)/sig)$value)
  
  mu0 = optim(init,h,method="L-BFGS-B",lower=a,upper=b)$par
  sig = sqrt(-2/numDeriv::hessian(h,mu0)[1])
  h2 = function(x) exp(h(mu0+sig*x)-h(mu0))
  lnI2 = h(mu0)+log(sig)+log(integrate(h2,lower=(a-mu0)/sig,upper=(b-mu0)/sig)$value)
  
  return(exp(lnI1)-k*exp(lnI2))
}

# C = 1/sqrt(2*pi)
# mom0 = function(x) x-x+1
# mom1 = function(x) x
# mom2 = function(x) x^2
# zkern = function (x) -x^2/2 # Kernel of standard normal distribution
# C*wandint(mom0,zkern,-10,10,50,-9) # P(-10<Z<10)=1
# (C*wandint(mom0,zkern,1,1.5,50,1.1)-(pnorm(1.5)-pnorm(1)))^2 # Should be 0
# C*wandint(mom1,zkern,-10,10,50,-9) # E(Z)=0
# C*wandint(mom2,zkern,-10,10,50,-9) # E(Z^2)=1