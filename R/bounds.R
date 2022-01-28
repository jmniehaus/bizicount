#takes in marginal cdf evals (u1 u2) and joint copula eval to impose bounds
frech.bounds = function(u1, u2, p, frech.min = 1e-7) {

  fmax = 1 - frech.min
  low  = pmax(u1 + u2 - 1, frech.min)
  up   = if(frech.min > 0) pmin(u1, u2, fmax) else pmin(u1, u2)

  out = pmin( pmax(p, low), up)

  return(out)
}

bound = function(x, low=1e-7, upper=T){
  if(low >= 1 || low <= 0)
    stop("Invalid lower boundary.")

  up = 1-low

  v = if(!upper)
    pmax(low, x)
  else
    pmin(up, pmax(low, x))

  return(v)
}

