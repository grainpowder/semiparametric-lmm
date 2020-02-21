wandint = function(f,h,log.value=TRUE,threshold=1e-8)
{
  # Returns approximate value of \int_0^\infy f(x)exp{h(x)}dx
  integrand = function(x) log(f(x)) + h(x)
  upper_lim = 1.0001
  while (exp(integrand(upper_lim)) > threshold) upper_lim = upper_lim + 1
  
  browser()
  neg_integrand = function(x) -integrand(x)
  mu0 = optimize(neg_integrand, c(0, upper_lim))$minimum
  sig0 = sqrt(-1/numDeriv::hessian(integrand, mu0)[1])
  
  I0_integrand = function(x) exp(integrand(mu0+sig0*x)-integrand(mu0))
  L = -mu0/sig0
  U = (upper_lim - mu0)/sig0
  lnI0 = log(integrate(I0_integrand, lower=L, upper=U)$value)
  
  if (log.value) return(log(sig0) + integrand(mu0) + lnI0)
  return(exp(log(sig0) + integrand(mu0) + lnI0))
}

# Test on the accuracy of wandint function
mom1 = function(x) x
mom2 = function(x) x^2
lnx = function(x) log(x)

# On gamma distribution
alpha = 10; beta = 10 # requires alpha > 1 to ensure log-concavity
kern = function (x) (alpha-1)*log(x)-beta*x
lnC = log(gamma(alpha)) - alpha*log(beta)
result_mean = exp(wandint(mom1,kern)-lnC)
result_var = exp(wandint(mom2,kern)-lnC) - result_mean^2
result_lnx = exp(wandint(lnx,kern)-lnC)
data.frame(true=c(alpha/beta, alpha/(beta^2), -log(beta)+digamma(alpha)), 
           wandint=c(result_mean, result_var, result_lnx),
           row.names = c("mean", "variance", "lnx"))
