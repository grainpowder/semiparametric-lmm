mer = function(y, w, prior=NULL, maxiter=500, tol=1e-4, n_grid=1e3)
{
  N = length(y)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sb2 = 1e3
    ax = 1e-2
    bx = 1e-2
    au = 1e-2
    bu = 1e-2
    ae = 1e-2
    be = 1e-2
  }
  else
  {
    sb2 = prior$sb2
    ax = prior$ax
    bx = prior$bx
    au = prior$au
    bu = prior$bu
    ae = prior$ae
    be = prior$be
  }
  
  # Define grids for griddy Gibbs routine
  width = max(w) - min(w)
  grid = seq(min(w)-width/10, max(w)+width/10, length.out=n_grid)
  
  # Construct spline bases
  K = min(round(log(length(unique(w)))), 35)
  knots = quantile(grid, seq(0,1,length.out=K+2)[-c(1,K+2)])
  Z = bs(grid, knots)
  
  # Initialize variational parameters
  sigx.ratio = ax/bx
  sigu.ratio = au/bu
  sige.ratio = ae/be
  munu.q = rep(0, 2+n_knots)
  signu.q = diag(rep(sb2, 2+n_knots))
  
}

set.seed(10)
N = 130
RR = 0.8
sig2x = 0.9
sig2v = sig2x/RR-sig2x
mux = 3
x = rnorm(N, mux, sqrt(sig2x))
w = rnorm(N, x, sqrt(sig2v))
f = function(x) x*sin(pi*x)
y = f(x) + rnorm(N)
plot(w,y,main="Noise",pch=19)
plot(w,y,main="Noise and Overshadowed Pattern",pch=19)
points(x,y, col=2,pch=19,cex=0.7)
