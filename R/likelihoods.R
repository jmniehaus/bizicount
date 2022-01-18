# Bivariate likelihood (univariate is lower down)
cop.lik = function(parms){
  
  ct = lapply(ct.ind, function(x) parms[x])
  

  # Loop over inputs to create list of marginal lambdas
  # ie, lam[1] = f1(x1, z1, w1), lam[2] = f2(x2, z2, w2)
  lam = mapply(
    function(par, design, offset, invlink){
      if(ncol(design) == 1)
        as.vector(invlink(design * par + offset))
      else
        as.vector(invlink(design %*% par + offset))
    },
    par = ct,
    design = X,
    offset = offset.ct,
    invlink = invlink.ct,
    SIMPLIFY = F
  )
  
  # Get zi parameters via subset
  zi = list()
  if (n.zi == 1)
    zi[[zipos]] = parms[zi.ind[[zipos]]]
  if (n.zi == 2){
    zi = lapply(zi.ind, function(x) parms[x])
  }

  # dispersion and dependence params
  disp = list()
  if (n.nb == 1)
    disp[[nbpos]] = exp(tail(parms, 2)[1])
  if (n.nb == 2) {
    disp[[1]] = exp(tail(parms, 3)[1])
    disp[[2]] = exp(tail(parms, 2)[1])
  }
  
  dep = tail(parms, 1)
  if (cop == 'gaus')
    dep = tanh(dep)
  
  # Create list of zero-inflation probability vectors
  psi = list()
  if(n.zi > 0){
    for(i in zipos) 
      psi[[i]] = as.vector(invlink.zi[[i]](Z[[i]] %*% zi[[i]] + offset.zi[[i]]))
  }
  
  
  # get list of arguments to marginal cdfs 
  margin.args = list()
  
  for (i in c(1, 2)) {
    margin.args[[i]] = switch(
      margins[i],
      "pois" = list(
        q      = y[, i],
        lambda = lam[[i]]
      ),
      "nbinom"  = list(
        q     = y[, i],
        mu    = lam[[i]],
        size  = disp[[i]]
      ),
      "zip" = list(
        q      = y[, i],
        lambda = lam[[i]],
        psi    = psi[[i]]
      ),
      "zinb" = list(
        q     = y[, i],
        mu    = lam[[i]],
        size  = disp[[i]],
        psi   = psi[[i]]
      )
    )
  }
  

  # evaluate copula pmf
  d = cop.pmf(
    margin.args = margin.args,
    margins = margins,
    cop = cop,
    dep = dep,
    pmf.min = pmf.min,
    frech.min = frech.min
  )

  # bound likelihoods for logging
  d = bound(d * weights, low = pmf.min) 
  
  #ll
  -sum(log(d))

}



# Univariate Likelihood
zic.ll = function(parms,
                  grad = F) {
  kx = ncol(X)
  kz = ncol(z)
  #setup parameters based on distribution
  if (dist == "pois") {
    ct = head(parms, kx)
    zi = tail(parms, kz)
  } else {
    ct = head(parms, kx)
    zi = parms[(kx + 1):(kx + kz)]
    disp = exp(tail(parms, 1))
  }
  
  #get linear predictor, probability, mean
  eta.ct = X %*% ct + offset.ct
  eta.zi = z %*% zi + offset.zi
  psi = make.link(link.zi)$linkinv(eta.zi)
  lam = make.link(link.ct)$linkinv(eta.ct)
  
  ll = withCallingHandlers(switch(
    dist,
    "pois" = dzip(
      y,
      lambda = lam,
      psi = psi,
      log = T
    ),
    "nbinom"  = dzinb(
      y,
      size = disp,
      psi = psi,
      mu = lam,
      log = T
    )
  ),
  warning = ll.suppress)
  
  ll =  sum(ll * weights)
  
  if (grad)
    attr(ll, "gradient") = get(paste0(dist, ".grad"))(parms = parms)
  
  return(-ll)
}


# Univariate zip gradient
pois.grad = function(parms,
                     grad = F) {
  
  kx = ncol(X)
  kz = ncol(z)
  
  ct = head(parms, kx)
  zi = tail(parms, kz)
  
  eta.ct = (X %*% ct) + offset.ct
  eta.zi = (z %*% zi) + offset.zi
  
  linklist.ct = make.link(link.ct)
  linklist.zi = make.link(link.zi)
  
  lam.ct = linklist.ct$linkinv(eta.ct)
  deta.ct = linklist.ct$mu.eta(eta.ct)
  
  psi = linklist.zi$linkinv(eta.zi)
  deta.zi = linklist.zi$mu.eta(eta.zi)
  
  
  denom.y0 = psi + (1 - psi) * exp(-lam.ct)
  
  d.zi = ifelse(y == 0,
                deta.zi * (1 - exp(-lam.ct)) / denom.y0,
                -deta.zi / (1 - psi))
  
  d.ct = ifelse(y == 0,
                -(1 - psi) * exp(-lam.ct) * deta.ct / denom.y0,
                deta.ct * (y / lam.ct - 1))
  
  d.zi = z * d.zi * weights
  d.ct = X * d.ct * weights
  
  
  return(-unname(colSums(cbind(d.ct, d.zi))))
}


#univariate zinb gradient
nbinom.grad = function(parms,
                       grad = F) {
  kx = ncol(X)
  kz = ncol(z)
  ct = head(parms, kx)
  zi = parms[(kx + 1):(kx + kz)]
  theta = exp(tail(parms, 1))
  
  eta.ct = (X %*% ct) + offset.ct
  eta.zi = (z %*% zi) + offset.zi
  
  linklist.ct = make.link(link.ct)
  linklist.zi = make.link(link.zi)
  
  lam.ct = linklist.ct$linkinv(eta.ct)
  deta.ct = linklist.ct$mu.eta(eta.ct)
  
  psi = linklist.zi$linkinv(eta.zi)
  deta.zi = linklist.zi$mu.eta(eta.zi)
  
  denom.y0 = dzinb(y, size = theta, psi = psi, mu = lam.ct) ^ (-1)
  dnbinom.zi = dnbinom(y, size = theta, mu = lam.ct)
  
  d.zi = ifelse(y == 0,
                denom.y0 * (deta.zi - deta.zi * dnbinom.zi),
                -deta.zi / (1 - psi))
  
  d.ct = ifelse(
    y == 0,
    -denom.y0 * deta.ct * (1 - psi) * (theta / (lam.ct + theta)) ^ (theta + 1),
    deta.ct * theta * (y - lam.ct) / (lam.ct * theta + lam.ct ^ 2)
  )
  
  d.th = ifelse(
    y == 0,
    denom.y0 * (1 - psi) * (theta / (theta + lam.ct)) ^ theta *
      (log(theta) - log(lam.ct + theta) + lam.ct / (lam.ct + theta)),
    digamma(y + theta) - digamma(theta) + 1 + log(theta) -
      log(theta + lam.ct) - (theta + y) / (theta + lam.ct)
  )
  
  d.zi = z * d.zi * weights
  d.ct = X * d.ct * weights
  d.th = theta * d.th * weights # multiply by theta due to exponential transform
  
  return(-unname(colSums(cbind(d.ct, d.zi, d.th))))
  
}


# function to help with getting two offsets instead of one from model frame
# get.offset = function(formula, model.frame){
#   name = grep("offset",
#               unlist(strsplit(as.character(formula), " \\+ ")),
#               value=T)
#
#   if(length(name) == 0) return(0) else return(as.vector(model.frame[, name]))
# }

# function for setting defaults to each optimizer
set.defaults = function(de.list, de.names, de.values) {
  for (i in seq_along(de.names)) {
    if (is.null(de.list[[de.names[i]]]))
      de.list[[de.names[i]]] = de.values[i]
  }
  return(de.list)
}



