#=================================================================
check_uni_args = function(pos = 1) {
  args = as.list(parent.frame(pos))

  with(args, {
    assert(!(any(sapply(
      list(X, y, z), is.null
    )) && is.null(fmla)),
    "Either `fmla` or `X`, `y` and `z` must be defined.",
    pos = -6)

    assert(!any(weights < 0),
           "`weights` cannot be negative.", pos = -6)

    assert(
      is.null(starts) || is.null(dim(starts)) && is.numeric(starts),
      "Argument `starts` must be a numeric vector.",
      pos = -6
    )

    assert(
      is.null(fmla) || (is.null(offset.ct) && is.null(offset.zi)),
      "Offset args are only used when inputting matrices (X, y, Z). To use an offset with formula, add ' + offset(varname)' to formula.",
      pos = -6,
      type = "warning"
    )

  })

  return(invisible(NULL))

}

#=================================================================
check_biv_args = function(pos = 1) {
  args = as.list(parent.frame(pos))


  with(args, {
    assert(
      length(margins) == 2 &&
        all(margins %in% c("nbinom", "pois", "zip", "zinb")),
      "Margins must be length 2, and be one of 'nbinom', 'pois', 'zip', 'zinb' (partial matching not supported for this arg).",
      pos = -6
    )

     assert(
       missing(na.action),
       "`na.action` is deprecated as there is no na.predict method for bizicount. All models implement na.omit.",
       type = "warning",
       pos = -2
     )

    assert(
      is.null(starts) || is.null(dim(starts)) && is.numeric(starts),
      "Argument `starts` must be a numeric vector.",
      pos = -6
    )

    assert(
      frech.min <= 1e-5 & frech.min >= 0,
      "Frechet lower bound must be in [0, .00001].",
      pos = -6
    )

    assert(
      pmf.min <= 1e-5 & pmf.min >= 0,
      "PMF lower bound must be in [0, .00001].",
      pos = -6
    )

    assert(is.logical(keep),
           "`keep` must be logical (T/F).", pos = -6)

    fmla.list = list(as.Formula(fmla1),
                     as.Formula(fmla2))
    fmla.zi  = sapply(fmla.list, function(x) length(x)[2] == 2)
    margin.zi = grepl("zi", margins)
    eq = fmla.zi == margin.zi

    assert(
      all(eq),
      paste0(
        "Structure for formula ",
        paste(which(!eq), collapse = " and "),
        " does not match margin ",
        paste(which(!eq), collapse = " and "),
        ". Either remove zero-inflation components from corresponding formula(s), or specify zero-inflated margin(s)."
      ),
      pos = -6
    )

  })

  return(invisible(NULL))
}

#=================================================================


# For checking inputs to PMF, CDF, Quantile, and rngs
# for univariate zero inflated dists
check_dist_args = function(pos = 1,
                           negbin = F,
                           recycle = F) {
  arglist = as.list(parent.frame(pos))

  on.exit(rm(arglist), add = T)

  if (negbin) {
    assert(
      is.null(arglist$prob) || is.null(arglist$mu),
      "Args `prob` and `mu` cannot be defined simultaneously."
    )
    assert(
      !is.null(arglist$prob) || !is.null(arglist$mu),
      "One of `prob` or `mu` must be specified."
    )
  }

  assert((!is.null(arglist$x) && all(arglist$x >= 0)) || !is.null(arglist$q) && all(arglist$q >=0) || (is.null(arglist$x) & is.null(arglist$q)),
         "Negative values of `q` or `x` supplied.",
         type="warning")

  assert(!missingArg(as.symbol("psi"), envir = parent.frame(pos), eval=T),
         "`psi` must be specified.")

  assert(
    is.logical(arglist$log) ||
      is.logical(arglist$log.p) || is.null(arglist$log),
    "Arg `log` or `log.p` must be logical."
  )
  assert(
    is.logical(arglist$lower.tail) || is.null(arglist$lower.tail),
    "Arg `lower.tail` must be logical."
  )


  # length of non-null args
  arglist.l = sapply(arglist[sapply(arglist, function(x)
    ! is.null(x))],
    length)

  arglist.u = length(unique(arglist.l))



  # make sure all args are same length, or that there are at most two lengths, with one of them being length one
  assert(is.logical(arglist$recycle),
         "`recycle` must be logical.")

  if (!recycle) {
    assert(
      arglist.u == 1 || (arglist.u == 2 && any(arglist.l == 1)),
      "Argument lengths imply ambiguous recycling. Ensure that all args are of
      same length, or that all lengths greater than one are equal."
    )
  } else {
    warning("`recycle=T` may give inaccurate results.")
  }

  return(invisible(NULL))
}

