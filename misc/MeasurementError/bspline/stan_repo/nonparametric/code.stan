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
  int num_data;             // number of data points
  int n_intknot;            // number of interior knots
  vector[n_intknot] knots;  // the sequence of interior knots
  vector[2] boundary;       // two boundary knots
  int order;                // order of the spline basis
  int resolution;           // number of grids to use
  real Y[num_data];
  real X[num_data];
  matrix[resolution, n_intknot + order] BsGrid;
}

transformed data {
  int num_basis = n_intknot + order;        // dimension of B-spline basis
  matrix[num_data, num_basis] B;            // initialize matrix of B-splines
  vector[order + n_intknot] ext_knots_temp;
  vector[2*order + n_intknot] ext_knots;    // initialize the set of extended knots
  ext_knots_temp = append_row(rep_vector(boundary[1], order), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(boundary[2], order));
  for (ind in 1:num_basis)
    B[:,ind] = to_vector(build_b_spline(X, to_array_1d(ext_knots), ind, order));
}

parameters {
  vector[num_basis] u;
  real beta;
  real<lower=0> sigma;
  real<lower=0> upsilon;
}

transformed parameters {
  vector[num_data] Y_hat;
  Y_hat = beta + to_vector(B*u);
}

model {
  // Priors
  u ~ normal(0, upsilon);
  beta ~ normal(0, 100);
  upsilon ~ inv_gamma(1e-3, 1e-3);
  sigma ~ inv_gamma(1e-3, 1e-3);

  //Likelihood
  Y ~ normal(Y_hat, sigma);
}

generated quantities {
  vector[resolution] fxGrid;
  fxGrid = BsGrid*u;
}
