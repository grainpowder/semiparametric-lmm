build_ext_knots = function(intKnots, boundary, order=4)
{
  ext_knots = c(rep(boundary[1], order), intKnots, rep(boundary[2], order))
  return(ext_knots)
}

build_b_spline = function(x, ext_knots, ind, order=4)
{
  # Function to construct basis vector of cubic B-spline, B_{i,order}(x)
  b_spline = w1 = w2 = rep(0, length(x))
  
  # Define coefficients of B-spline basis recursively
  if (order == 1) {
    # Base case
    for (i in 1:length(x)) b_spline[i] = (ext_knots[ind] <= x[i]) && (x[i] < ext_knots[ind+1])
  } else {
    # Implementation of sub-problem
    if (ext_knots[ind] != ext_knots[ind+order-1]) {
      w1 = (x - rep(ext_knots[ind],length(x))) / (ext_knots[ind+order-1] - ext_knots[ind])
    }
    if (ext_knots[ind+1] != ext_knots[ind+order]) {
      w2 = (rep(ext_knots[ind+order],length(x)) - x) / (ext_knots[ind+order] - ext_knots[ind+1])
    }
    b_spline = w1*build_b_spline(x, ext_knots, ind, order-1) + w2*build_b_spline(x, ext_knots, ind+1, order-1)
  }
  return(b_spline)
}

# Test on its accuracy
library(splines)
intKnots = c(3,7)
boundary = c(0,10)
x = 1:9

ext_knots = build_ext_knots(intKnots, boundary)
B = matrix(0, nrow=length(x), ncol=(4-1)+length(intKnots)+1)
for (i in 1:ncol(B)) B[,i] = build_b_spline(1:9, ext_knots, i)

norm(bs(x, knots=intKnots, intercept=TRUE, Boundary.knots=boundary), "F") - norm(B, "F")
