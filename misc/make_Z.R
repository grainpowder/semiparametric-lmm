make_Z = function(nrows) {
  # nrows: number of replicates of each object(i.e. t_i, i=1,...,n)
  # returns a matrix composed of 0 or 1
  nrows = cumsum(nrows)
  Z = matrix(0,nrows[length(nrows)],length(nrows))
  Z[1:nrows[1],1] = 1
  for(col in 2:ncol(Z)) Z[(nrows[col-1]+1):nrows[col],col] = 1
  return(Z)}