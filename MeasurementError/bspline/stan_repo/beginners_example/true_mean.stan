data {
    int N;          // Length of the data
    vector[N] y;    // Data whose true mean value is to be estimated
}
transformed data {
    vector[N] ones;
    ones = rep_vector(1, N);
}
parameters {
    real mu;                // True mean of the data
    real<lower=0> sigma;    // Variance of the data
}
model {
    mu ~ normal(0, 100);
    sigma ~ inv_gamma(1e-3, 1e-3);
    for (n in 1:N)
        y[n] ~ normal(mu * ones[n], sigma);
}
