cop.cdf = function(margin.args,
                   margins = c("pois", "pois"),
                   dep,
                   cop = "frank",
                   frech.min = 1e-7) {

  # evaluate marginal cdfs using supplied marginal cdf functions and arguments
  u = mapply(function(margin, args)
              do.call(paste0("p", margin), args),
             margin =margins,
             args =margin.args,
             SIMPLIFY = F)

  if (cop == "frank") {
    f = lapply(u, function(x)
      exp(-dep * x) - 1)

    den = (exp(-dep) - 1)

    p = -dep ^ (-1) * log(1 + f[[1]] * f[[2]] / den)

  } else if (cop == "gaus") {

    u = lapply(u, bound, low = 1e-16, upper = T) #avoid inf in qnorm, as this ends up giving NaN in pbivnorm
    qu = lapply(u, qnorm)
    p = pbivnorm(x = qu[[1]], y = qu[[2]], rho = dep)

  }


  p = frech.bounds(
    p = p,
    u1 = u[[1]],
    u2 = u[[2]],
    frech.min = frech.min
  )

  return(p)
}

cop.pmf = function(margin.args,
                   dep,
                   margins = c("pois", "pois"),
                   cop = "frank",
                   pmf.min = 1e-7,
                   frech.min = 1e-7) {

  # gets new quantiles for finite differencing (F(y-1)),
  # puts them into nested list for copula cdf evaluation
  new.args = rep(margin.args, 4)
  for (i in c(3,6,7,8)) {
    new.args[[i]][["q"]] = new.args[[i]][["q"]] - 1
  }

  # takes quantiles after subtracting
  # and puts them into a list for the cdf function
  cdf.args = list()
  iter = 1
  for (i in seq_len(4)) {
    cdf.args[i] = list(new.args[iter:(iter + 1)])
    iter = iter + 2
  }

  # apply cdf to all quantiles, imposing F bounds if indicated
  p = lapply(cdf.args,
             function(w)
               cop.cdf(
                 margin.args = w,
                 dep = dep,
                 cop = cop,
                 margins = margins,
                 frech.min = frech.min
               ))

  d = p[[1]] - p[[2]] - p[[3]] + p[[4]]

  return(d)

}


