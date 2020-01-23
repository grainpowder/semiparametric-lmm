source("onlinevb.R")
source("logdeter.R")

set.seed(1)
# Setting the size of the simulation --------------------------------------

# n = 400
n = 10
T = 5
p = T
sigma2 = 1
# K = 256^2
K = 16

# Generate beta -----------------------------------------------------------

# Beta_i is generated from mixed distribution of two normal distribution of:
#   - N1: N([1.5, 1.5, 1, 2, 2], Simga)
#   - N2: N([-1.5, -1.5, -1, -2, -2], Simga)
# where Sigma is sampled from W(I,10) distribution.
# Resulting betai is a data.frame of n*p dimension, where each row contains time-varying effect of a subject

mu1 = c(1.5, 1.5, 1, 2, 2)
betai = matrix(0, n, p)
vaco = MCMCpack::riwish(10, diag(p))
for (i in 1:n) 
{
  # Generate and assign beta from N1 if runif(1) > .5, from N2 elsewhere
  fromN1 = runif(1) > 0.5
  if (fromN1) 
  {
    betai[i,] = mvtnorm::rmvnorm(1, mu1, vaco)
    next 
  }
  betai[i,] = mvtnorm::rmvnorm(1, -mu1, vaco)
}


# Generate Y and X --------------------------------------------------------

# X_{it} is set to be K*T matrix whose column is an indicator of time.
# ex. X_{43} is indicator matrix of 4th subject whose every element is 0 except third column.
for (i in 1:n)
{
  Yi = matrix(0, K, T)
  etai = rnorm(K)
  for (t in 1:T)
  {
    posi = matrix(0, 1, T)
    posi[1, t] = 1
    Xit = kronecker(rep(1, K), posi)
    Yi[, t] = (Xit%*%betai[i,]) + etai + rnorm(K, 1, sqrt(sigma2))
    filename = paste0("datasets/X_", paste0(i, t), ".csv")
    write.csv(Xit, filename, row.names=FALSE)
  }
  filename = paste0("datasets/Y_", i, ".csv")
  write.csv(Yi, filename, row.names=FALSE)
}
rm("Xit", "Yi")

# Construct proximity matrix and C, Omega accordingly ---------------------

indi = 1:K
rind = matrix(rep(1:sqrt(K), sqrt(K)), ncol=1)
cind = matrix(rep(1:sqrt(K), each=sqrt(K)), ncol=1)
DD = cbind(rind, cind)
rm("rind", "cind")
Omega = rep(0, K)
Cmat = matrix(NA, 0, 3)
for (k in 1:K)
{
  # phi was set to be 2.16 in this simulation
  dis = sqrt((DD[k, 1]-DD[, 1])^2 + (DD[k, 2]-DD[, 2])^2)^(-2.16)
  dis[is.infinite(dis)] = 0
  Omega[k] = 1 / sum(dis)
  dis = dis / sum(dis)
  rind = k * rep(1, sum(dis >= 0.01))
  cind = indi[dis >= 0.01]
  dist = dis[dis >= 0.01]
  Cmat = rbind(Cmat, cbind(rind, cind, dist))
}
write.csv(Cmat, "datasets/Cmat.csv", row.names=FALSE)
write.csv(Omega, "datasets/Omega.csv", row.names=FALSE)
rm("Omega", "Cmat", "dis", "rind", "cind", "DD")

# Parameters for priors ---------------------------------------------------

asigma = bsigma = 1
aalpha = balpha = 1
atau = btau = 1/10000
beta0 = rnorm(p)
capsigma0 = diag(p)
R = 20 # DPM truncation level
npi = 40 # rho related...
pitl = (1/npi) * rep(1, npi) # also rho related...


# Initial values for parameters -------------------------------------------

asigmatls = bsigmatls = aalphatls = balphatls = 1
atautls = btautls = 0.00001
betatls = matrix(0, p, R)
capsigmatls = array(0, c(R, p, p))
for (r in 1:p)
{
  capsigmatls[, ,r][, r] = rep(0.1, R)
}
gama1s = rep(2, R-1)
gama2s = rep(1, R-1)
wbs = rep(0, npi)
pitls = pitl
kappas = runif(R)
kappas = kappas / sum(kappas)
capsis = rep(0.01, K)
rho = seq(0, 0.9, length.out=npi)
# so3nr...logdeter?





