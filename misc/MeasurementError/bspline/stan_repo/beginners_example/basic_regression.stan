data {
    int N;          // Length of the data
    int D;          // Number of the explanatory variables
    matrix[N, D] w; // Matrix of explanatory variable
    vector[N] y;    // Data whose true mean value is to be estimated
}
transformed data {
    vector[N] ones;
    matrix[N, D+1] W;
    ones = rep_vector(1, N);
    W = append_col(ones, w);    // Attached one vector to include intercept
}
parameters {
    vector[D+1] beta;       // Coefficients
    real<lower=0> sigma;    // Variance of the data
}
model {
    sigma ~ inv_gamma(1e-3, 1e-3);
    beta ~ normal(0, 100);
    y ~ normal(W * beta, sigma);
}
