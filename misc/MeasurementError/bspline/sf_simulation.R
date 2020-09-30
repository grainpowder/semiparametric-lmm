library(rstan)

make_Z = function(nrows) 
{
    # nrows: number of replicates of each object(i.e. t_i, i=1,...,n)
    # returns a matrix composed of 0 or 1
    nrows = cumsum(nrows)
    Z = matrix(0,nrows[length(nrows)],length(nrows))
    Z[1:nrows[1],1] = 1
    for(col in 2:ncol(Z)) Z[(nrows[col-1]+1):nrows[col],col] = 1
    return(Z)
}

# Data generation
set.seed(10)
N = 50
T = 7
D = 6
beta = rnorm(D+1)
b = rexp(N)
Z = make_Z(rep(T,N))
w = matrix(rnorm(N*T*D),N*T,D)
y = drop(cbind(1,w)%*%beta + Z%*%b + rnorm(N, sd=0.7))

# MCMC
data = list(
    n_obs = length(y),
    n_regcoef = D,
    n_objects = N,
    sign = 1,
    y = y,
    w = w,
    Z = Z)

# Apply stan code
set.seed(100)
mcmc_result = stan("./stan_repo/sf/code.stan", data=data, iter=1e3)

# Result
extracted_beta = extract(mcmc_result, 'beta')$beta
beta_mean = apply(extracted_beta, 2, mean)
beta_upper = apply(extracted_beta, 2, quantile, 0.975)
beta_lower = apply(extracted_beta, 2, quantile, 0.025)
plot(beta, beta_mean, main='Fixed Effects', xlab='True', ylab='Estimated', pch=19)
lines(-10:10, -10:10)

extracted_b = extract(mcmc_result, 'b')$b
b_mean = apply(extracted_b, 2, mean)
b_upper = apply(extracted_b, 2, quantile, 0.975)
b_lower = apply(extracted_b, 2, quantile, 0.025)
plot(b, b_mean, main='Random Effects', xlab='True', ylab='Estimated', pch=19, cex=0.5)
lines(-10:10, -10:10)
