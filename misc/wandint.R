wandint = function(f,h,a,b,k,init,log.value=TRUE){
  # Requires numDeriv package while calculating second derivative of h
  # Calculate integral (*) numerically by the method introduced in Appendix B of Wand(2011).
  # 1. Note that h(x) has to be concave function (i.e. h''(x) < 0 for all a < x <b)
  # 2. As this function uses built-in `optim` function to differentiate h, init s.t. a < init < b is required.
  # 3. To prevent numerical underflow, ln{f(x)+k}+h(x) is integrated alternatively.
  #
  # Let I2 = \int_a^b exp{h(x)}dx(i.e. normalizing constant of a distribution whose kernel is exp{h(x)})
  # (*) E_X[f(X)] = \int_a^b f(x)exp{h(x)}dx / I2
  #               = \int_a^b {f(x)+k-k}exp{h(x)}dx / I2
  #               = (\int_a^b exp[ln{f(x)+k}+h(x)]dx - k\int_a^b exp{h(x)}dx) / I2
  #               = (I1                              - k*I2) / I2
  #                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ what this function returns
  integrand = function(x) log(f(x)+k)+h(x)
  mu0 = optim(init,integrand,method="L-BFGS-B",lower=a,upper=b)$par
  sig = sqrt(-2/numDeriv::hessian(integrand,mu0)[1])
  h1 = function(x) exp(integrand(mu0+sig*x)-integrand(mu0))
  lnI1 = integrand(mu0)+log(sig)+log(integrate(h1,lower=(a-mu0)/sig,upper=(b-mu0)/sig)$value)
  
  mu0 = optim(init,h,method="L-BFGS-B",lower=a,upper=b)$par
  sig = sqrt(-2/numDeriv::hessian(h,mu0)[1])
  h2 = function(x) exp(h(mu0+sig*x)-h(mu0))
  lnI2 = h(mu0)+log(sig)+log(integrate(h2,lower=(a-mu0)/sig,upper=(b-mu0)/sig)$value)
  result = exp(lnI1)-k*exp(lnI2)
  if (log.value) return(log(result))
  return(result)
}

# # Test on the accuracy of wandint function
# mom0 = function(x) x-x+1
# mom1 = function(x) x
# mom2 = function(x) x^2
# 
# # On normal distribution
# mu = 0.65; sig = 1.9
# lnC = log(sqrt(2*pi)*sig)
# kern = function (x) -((x-mu)/sig)^2/2 # Kernel of normal distribution
# numprob = exp(wandint(mom0,kern,6,8,50,1.1)-lnC)
# realprob = pnorm(8,mu,sig)-pnorm(6,mu,sig)
# result_mean = exp(wandint(mom1,kern,-10,10,50,-9)-lnC)
# result_var  = exp(wandint(mom2,kern,-10,10,50,-9)-lnC) - result_mean^2
# data.frame(true=c(realprob,mu,sig^2), wandint=c(numprob,result_mean,result_var),
#            row.names = c("some tail probability", "mean", "variance"))
# 
# On gamma distribution
# alpha = 1.7; beta = 4.1 # requires alpha > 1 to ensure log-concavity
# lnC = log(gamma(alpha)) - alpha*log(beta)
# kern = function (x) (alpha-1)*log(x)-beta*x
# numprob = exp(wandint(mom0,kern,2,4,50,1.1)-lnC)
# realprob = pgamma(4,shape=alpha,rate=beta)-pgamma(2,shape=alpha,rate=beta)
# result_mean = exp(wandint(mom1,kern,0.01,10,50,0.01)-lnC)
# result_var = exp(wandint(mom2,kern,0.01,10,50,0.01)-lnC) - result_mean^2
# data.frame(true=c(realprob,alpha/beta,alpha/(beta^2)), wandint=c(numprob,result_mean,result_var),
#            row.names = c("some tail probability", "mean", "variance"))
