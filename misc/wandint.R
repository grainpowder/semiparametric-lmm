wandint = function(f,h,log.value=TRUE,threshold=1e-12,shrink_a=1e-4)
{
  # Returns approximate value of \int_0^\infy f(x)exp{h(x)}dx.
  # Case is divided to ensure quadrature on both f(0) > 0, f(0) <= 0 cases.
  # f(0) > 0(shrink = FALSE) is basic case. f(0) <= 0 is just a replication of this case.
  shrink = f(0) <= 0
  upper_lim = 1
  if (!shrink)
  {
    integrand = function(x) log(f(x)) + h(x)
    while (exp(integrand(upper_lim)) > threshold) upper_lim = upper_lim + 1
    
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
  else
  {
    k = abs(f(shrink_a))
    integrand = function(x) log(f(x)+k) + h(x)
    while (exp(integrand(upper_lim)) > threshold) upper_lim = upper_lim + 1
    neg_integrand = function(x) -integrand(x)
    mu0 = optimize(neg_integrand, c(shrink_a, upper_lim))$minimum
    sig0 = sqrt(-1/numDeriv::hessian(integrand, mu0)[1])
    I1_integrand = function(x) exp(integrand(mu0+sig0*x)-integrand(mu0))
    L = (shrink_a-mu0) / sig0
    U = (upper_lim-mu0) / sig0
    lnI0_1 = log(integrate(I1_integrand, lower=L, upper=U)$value)
    I1 = exp(log(sig0) + integrand(mu0) + lnI0_1)
    
    upper_lim = 1
    while (exp(h(upper_lim)) > threshold) upper_lim = upper_lim + 1
    neg_h = function(x) -h(x)
    mu0 = optimize(neg_h, c(shrink_a, upper_lim))$minimum
    sig0 = sqrt(-1/numDeriv::hessian(h, mu0)[1])
    I2_integrand = function(x) exp(h(mu0+sig0*x)-h(mu0))
    L = (shrink_a-mu0) / sig0
    U = (upper_lim - mu0) / sig0
    lnI0_2 = log(integrate(I2_integrand, lower=L, upper=U)$value)
    I2 = exp(log(sig0) + h(mu0) + lnI0_2)
    
    if (log.value) return(log(I1 - k*I2))
    return(I1 - k*I2)
  }
}

# Test on the accuracy of wandint function
mom0 = function(x) x-x+1
logx = function(x) log(x)
mom1 = function(x) x
mom2 = function(x) x^2

# On gamma distribution
alpha = 1.1; beta = 0.3 # integration w.r.t. p.d.f. returns NA/NaN error when alpha, beta < 1.
kern = function (x) (alpha-1)*log(x)-beta*x
lnC = wandint(mom0,kern)
result_logx = exp(wandint(logx,kern)-lnC)
result_mean = exp(wandint(mom1,kern)-lnC)
result_var = exp(wandint(mom2,kern)-lnC) - result_mean^2
data.frame(true=c(-log(beta)+digamma(alpha), alpha/beta, alpha/(beta^2)), 
           wandint=c(result_int, result_mean, result_var),
           row.names = c("logx", "mean", "variance"))

