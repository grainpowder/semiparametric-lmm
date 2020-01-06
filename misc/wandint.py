import numpy as np
from numdifftools import Hessian
from scipy.optimize import minimize
from scipy.integrate import quad


def wandint(f, h, lo, up, k, init, log=True):
    """
        Calculates integral (*) numerically by the method introduced in Appendix B of Wand(2011).
        1. Note that h(x) has to be concave function (i.e. h''(x) < 0 for all a < x <b)
        2. As this function uses minimize function in scipy, integrand is converted into convex function when finding extreme values
        3. To prevent numerical underflow, ln{f(x)+k}+h(x) is integrated alternatively.

        Let I2 = \int_lo^up exp{h(x)}dx(i.e. normalizing constant of a distribution whose kernel is exp{h(x)})
        (*) E_X[f(X)] = \int_a^b f(x)exp{h(x)}dx / I2
                      = \int_a^b {f(x)+k-k}exp{h(x)}dx / I2
                      = (\int_a^b exp[ln{f(x)+k}+h(x)]dx - k\int_a^b exp{h(x)}dx) / I2
                      = (I1                              - k*I2) / I2
                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ what this function returns
        
        Note : Unlike R, this Python implementation occasionally emits RuntimeWarning message when calculating Hessian.
               So suppressing those message using `warning` module is recommended.
    """
    integrand = lambda x: -np.log(f(x)+k) - h(x)
    mu0  = minimize(integrand, init, method="L-BFGS-B", bounds=[(lo,up)]).x[0]
    sig  = np.sqrt(2/Hessian(integrand)(mu0))
    h1   = lambda x: np.exp(integrand(mu0) - integrand(mu0+sig*x))
    lnI1 = -integrand(mu0) + np.log(sig) + np.log(quad(h1, a=(lo-mu0)/sig, b=(up-mu0)/sig)[0])

    neg_h = lambda x: -h(x)
    mu0  = minimize(neg_h, init, method="L-BFGS-B", bounds=[(lo,up)]).x[0]
    sig  = np.sqrt(2/Hessian(neg_h)(mu0))
    h2   = lambda x: np.exp(h(mu0+sig*x) - h(mu0))
    lnI2 = h(mu0) + np.log(sig) + np.log(quad(h2, a=(lo-mu0)/sig, b=(up-mu0)/sig)[0])
    
    result = np.exp(lnI1) - k*np.exp(lnI2)
    if log:
        return(np.log(result))
    else:
        return(result)


# # Test on the accuracy of wandint function
# from scipy.stats import norm, gamma
# from scipy.special import gamma as G
# import warnings
# warnings.filterwarnings("ignore")

# const = lambda x: x-x+1
# x1 = lambda x: x
# x2 = lambda x: x**2

# # On normal distribution
# mu = 0.65; sig = 1.9; low = 6; upp = 8
# true_prob = norm.cdf(upp,mu,sig)-norm.cdf(low,mu,sig)
# true = np.array([mu, sig, true_prob])

# lnC = np.log(np.sqrt(2*np.pi)*sig)
# kern = lambda x: -((x-mu)/sig)**2/2
# num_prob = np.exp(wandint(const, kern, low, upp, 50, 1.1)-lnC)[0]
# num_mean = np.exp(wandint(x1,kern,-10,10,50,-9)-lnC)[0]
# num_var  = (np.exp(wandint(x2,kern,-10,10,50,-9)-lnC) - num_mean**2)[0]
# numerical = np.concatenate((num_mean, np.sqrt(num_var), num_prob), axis=0)

# np.allclose(true, numerical)

# # On gamma distribution
# alpha = 1.7; beta = 4.1; low = 2; upp = 4
# true_prob = gamma.cdf(upp, a=alpha, scale=1/beta)-gamma.cdf(low, a=alpha, scale=1/beta)
# true = np.array([alpha/beta, alpha/beta**2, true_prob])

# lnC = np.log(G(alpha)) - alpha*np.log(beta)
# kern = lambda x: (alpha-1)*np.log(x) - beta*x
# num_prob = np.exp(wandint(const, kern, low, upp, 50, 2.1)-lnC)[0]
# num_mean = np.exp(wandint(x1,kern,0.01,10,50,0.1)-lnC)[0]
# num_var  = (np.exp(wandint(x2,kern,0.01,10,50,0.1)-lnC) - num_mean**2)[0]
# numerical = np.concatenate((num_mean, num_var, num_prob), axis=0)

# true; numerical