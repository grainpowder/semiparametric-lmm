wandint = function(f,h,a,b,inflator,init){
  # Calculate following integral numerically introduced in Appendix B of Wand(2011).
  # Note that h(x) has to be concave function (ex. standard normal(-x^2/2 -> -1))
  #
  # E_X[f(X)] = \int_a^b f(x)exp{h(x)}dx
  #
  # As this function uses built-in `optim` function to differentiate h, init s.t. a < init < b is required.
  # `inflator` is used to prevent numerical underflow. Check the usage below.
  mu0 = optim(init,h,method="L-BFGS-B",lower=a,upper=b)$par
  sig = sqrt(-2/hessian(h,mu0))[1]
  h1 = function(x) exp(h(mu0+sig*x)-h(mu0))
  lnI1 = h(mu0)+log(sig)+log(integrate(h1,lower=(a-mu0)/sig,upper=(b-mu0)/sig)$value)
  return(exp(lnI1))
}
h = function(x) -x^2/2 # kernel of standard normal distribution
wandint(1,h,-1,1,0,-0.9)/sqrt(2*pi) # numerical integration of P(-1<Z<1)
pnorm(1)-pnorm(-1) # Actual value of P(-1<Z<1)
