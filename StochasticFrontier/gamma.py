from misc import wandint
import numpy as np


const = lambda x: x-x+1
x1 = lambda x: x
x2 = lambda x: x**2

# On gamma distribution
from scipy.stats import norm, gamma
from scipy.special import gamma as G
import warnings
warnings.filterwarnings("ignore")

alpha = 1.7; beta = 4.1; low = 2; upp = 4
true_prob = gamma.cdf(upp, a=alpha, scale=1/beta)-gamma.cdf(low, a=alpha, scale=1/beta)
true = np.array([alpha/beta, alpha/beta**2, true_prob])

lnC = np.log(G(alpha)) - alpha*np.log(beta)
kern = lambda x: (alpha-1)*np.log(x) - beta*x
num_prob = np.exp(wandint(const, kern, low, upp, 50, 2.1)-lnC)[0]
num_mean = np.exp(wandint(x1,kern,0.01,10,50,0.1)-lnC)[0]
num_var  = (np.exp(wandint(x2,kern,0.01,10,50,0.1)-lnC) - num_mean**2)[0]
numerical = np.concatenate((num_mean, num_var, num_prob), axis=0)

true; numerical