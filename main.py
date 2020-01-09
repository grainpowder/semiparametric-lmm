import StochasticFrontier as SF
import numpy as np
from misc import makeZ

# Simulation
gammaSF = SF.Gamma()
n = 20; ti = 5
N = n*ti; D = 10
nrows = np.ones(n, dtype=int)*ti
Z = makeZ(nrows)
sigma = 0.5
theta = 2
lamb = 2

np.random.seed(1)
betaT = np.random.normal(0, sigma, D+1)
uT = np.random.gamma(theta, 1/lamb, size=n)
w = np.random.normal(0, 1, N*D).reshape(N,D)
y = np.hstack((np.ones(N).reshape(N,1),w))@betaT+Z@uT+np.random.normal(0,1,N)

gammaSF.fit(y,w,Z)
