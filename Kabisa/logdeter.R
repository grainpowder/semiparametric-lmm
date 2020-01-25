logdeter = function(p,m,K,Cmat,rho)
{
  # Computes log(det(I-rho*C)) using method described in Barry and Pace(1999).
  if (!("pracma" %in% row.names(installed.packages()))) 
  {
    stop("Required package 'pracma' is not installed yet.")
  }
  vv = matrix(0, length(rho), p)
  for (cind in 1:p)
  {
    v = matrix(0, length(rho), 1)
    x = rnorm(K, 1)
    H = x
    for (k in 1:m)
    {
      H = pracma::accumarray(Cmat[, 1], H[Cmat[, 2]] * Cmat[, 3])
      v = (K*(rho^k) * sum(H*x)) / k + v
    }
    v = v / sum(x^2)
    vv[, cind] = v
  }
  return(apply(vv, mean, 1))
}