from misc import wandint
from scipy.special import pbdv, gamma, digamma
from numpy.linalg import inv
import numpy as np
import pdb

__all__ = ['Gamma']

class Gamma:


    def __init__(
        self,
        productivity=True,
        sigbeta2=100,
        rsig_0=1.1,
        ssig_0=1.1,
        rthe_0=1.1,
        sthe_0=1.1,
        slam_0=1.1
    ):
        self.sgn = (-1) ** productivity
        self.sigbeta2 = sigbeta2
        self.rsig_0 = rsig_0
        self.ssig_0 = ssig_0
        self.rthe_0 = rthe_0
        self.sthe_0 = sthe_0
        self.slam_0 = slam_0
    

    def moment_u(self, values):
        """
        Takes a tuple of parameters as input and returns m-th moment of u_i.
        A tuple should envelop following parameters in designated order as:
            (theta, m, mu_i, v_i) # `v_i`, not `v_i^2`
        """
        theta, m, mu, var = values
        numerator = gamma(theta+m)*pbdv(-theta-m, mu/var)[0]
        denominator = var**m*gamma(theta)*pbdv(-theta, mu/var)[0]
        return(numerator/denominator)


    def fit(self, y, w, Z, tol=1e-6, maxiter=500):
        """
        Estimates parameters composing VB-SFA with gamma distributed inefficiency contained
        """
        # Arrays/values to be frequently used
        N, D = w.shape
        n = Z.shape[1]
        W = np.hstack((np.ones(N).reshape(N, 1), w))
        WtW = W.T@W
        ZtZ = Z.T@Z
        ti = np.apply_along_axis(np.sum, 0, Z)
        sb2diag = np.diag(np.ones(D+1)*(1/self.sigbeta2))
        rsig_q = self.rsig_0+(N/2)
        rep_n = np.ones(n)
        # Storages to keep updated value
        sig_ratio, lam_ratio, the_ratio = [1]*3
        mubeta_q = np.zeros(D+1)
        muu_q = rep_n*0
        lnu = rep_n*0
        # Predefine function f to be applied into numerical integration
        const = lambda x: x-x+1
        x1 = lambda x: x
        lx = lambda x: np.log(x)
        # Iterate
        for iter in range(maxiter):
            mubeta_q_old = mubeta_q
            muu_q_old = muu_q
            # beta
            sigbeta_q = sig_ratio*WtW+sb2diag
            sigbeta_q = inv(sigbeta_q)
            mubeta_q = sig_ratio*sigbeta_q@W.T@(y-self.sgn*Z@muu_q)
            # u
            mu = -lam_ratio+self.sgn*sig_ratio*Z.T@(y-W@mubeta_q)
            v2 = ti*sig_ratio/2
            theta = rep_n*the_ratio
            args1 = [tup for tup in zip(theta, rep_n, mu, np.sqrt(v2))]
            args2 = [tup for tup in zip(theta, rep_n*2, mu, np.sqrt(v2))]
            muu_q = np.array([result for result in map(self.moment_u, args1)])
            sigu_q = np.array([result for result in map(self.moment_u, args2)])-muu_q**2
            for i in range(n):
                integrand = lambda u: (the_ratio-1)*np.log(u)+mu[i]*u-v2[i]*u**2
                lnC = wandint(const, integrand, 1e-2, 10, 50.1, 0.1, log=False)
                lnu[i] = wandint(lx, integrand, 1e-2, 10, 50.1, 0.1, log=False)/lnC
            # lambda
            rlam_q = (n+1)*the_ratio
            slam_q = self.slam_0+np.sum(muu_q)
            lam_ratio = rlam_q/slam_q
            llamb = -np.log(slam_q)+digamma(rlam_q)
            # theta
            integrand = lambda th: -(n+1)*np.log(gamma(th))+(np.sum(lnu)+(n+1)*llamb+np.log(self.slam_0)-self.sthe_0)*th+(self.rthe_0-1)*np.log(th)
            lnC = wandint(const, integrand, 1e-2, 10, 50.1, 0.1)
            the_ratio = np.exp(wandint(x1, integrand, 1e-2, 10, 50.1, 0.1)-lnC)
            # sigma
            ssterm = np.sum((y-W@mubeta_q-self.sgn*Z@muu_q) ** 2)
            trterm1 = np.sum(np.diag(WtW@sigbeta_q))
            trterm2 = np.sum(ZtZ@sigu_q)
            ssig_q = self.ssig_0+(ssterm+trterm1+trterm2)/2
            sig_ratio = rsig_q/ssig_q
            # Sanity check
            values = np.array(list(mubeta_q)+list(muu_q)+list(lnu))
            if np.sum(np.isnan(values)) > 0:
                print(f">>> While executing iteration after {iter}-th loop,")
                raise ValueError("NaN occured.")
            # Convergence
            rmse1 = np.sqrt(np.mean((mubeta_q-mubeta_q_old)**2))
            rmse2 = np.sqrt(np.mean((muu_q-muu_q_old)**2))
            if rmse1 < tol and rmse2 < tol:
                break
        self.mubeta_q = mubeta_q
        self.sigbeta_q = sigbeta_q
        self.muu_q = muu_q
        self.sigu_q = sigu_q
        self.rsig_q = rsig_q
        self.ssig_q = ssig_q
        self.Elam = lam_ratio
        self.Etheta = the_ratio