functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // FROM  : https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = (rep_vector(ext_knots[ind+order], size(t)) - to_vector(t)) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int n_obs;                 // number of data points
  int n_intknot;            // number of interior knots
  int n_regcoef;            // number of regression coefficients
  int n_objects;            // number of different index i
  vector[n_intknot] int_knots;  // the sequence of interior knots
  vector[2] boundary;       // two boundary knots
  int order;                // order of the spline basis
  int resolution;           // number of grids to use to report the result
  real y[n_obs];
  real v[n_obs];
  matrix[resolution, n_intknot + order] bsgrid;
  matrix[n_obs, n_regcoef] w;
  matrix[n_obs, n_objects] Z;
}

transformed data {
  int num_basis = n_intknot + order;        // dimension of B-spline basis
  vector[order + n_intknot] ext_knots_temp;
  vector[2*order + n_intknot] ext_knots;    // initialize the set of extended knots
  vector[n_obs] ones;
  matrix[n_obs, n_regcoef+1] W;
  ext_knots_temp = append_row(rep_vector(boundary[1], order), int_knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(boundary[2], order));
  ones = rep_vector(1, n_obs);
  W = append_col(ones, w);    // Attached one vector to include intercept
}

parameters {
  vector[n_objects] b;
  vector[num_basis] u;
  vector[n_regcoef+1] beta;
  real x[n_obs];
  real mu;
  real<lower=0> xi;
  real<lower=0> nu;
  real<lower=0> sigma;
  real<lower=0> phi;
  real<lower=0> upsilon;
}

transformed parameters {
  matrix[n_obs, num_basis] B;
  vector[n_obs] y_hat;
  for (ind in 1:num_basis)
    B[:,ind] = to_vector(build_b_spline(x, to_array_1d(ext_knots), ind, order));
  y_hat = to_vector(W*beta) + to_vector(Z*b) + to_vector(B*u);
}

model {
  // Priors
  b ~ normal(0, phi);
  u ~ normal(0, upsilon);
  beta ~ normal(0, 100);
  x ~ normal(mu, xi);
  mu ~ normal(0, 100);
  phi ~ inv_gamma(1e-3, 1e-3);
  xi ~ inv_gamma(1e-3, 1e-3);
  nu ~ inv_gamma(1e-3, 1e-3);
  sigma ~ inv_gamma(1e-3, 1e-3);
  upsilon ~ inv_gamma(1e-3, 1e-3);

  //Likelihood
  y ~ normal(y_hat, sigma);
  v ~ normal(x, nu);
}

generated quantities {
  vector[resolution] fxGrid;
  fxGrid = bsgrid*u;
}
