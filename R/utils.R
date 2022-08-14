# used for print methods
divider = function(symbol, width, prepend=F, append=F){
     if(prepend) cat("\n")
     cat(paste0(rep(symbol, width), collapse=""), "\n")
     if(append) cat("\n")
}

#used for print methods
modlabel = function(name, pad=2, ...){
     cat(paste(name, "|"), "\n")
     divider("-", (nchar(name) + pad), ...)
}


# Wrapper around match.arg to signal more informative errors.
match_arg = function(arg,
                     choices,
                     several.ok = FALSE,
                     max_choice = 2,
                     pos = -1) {

     cal = sys.call(pos)
     def = sys.function(pos)

     call = match.call(def, cal)

     msg = paste("Arg ",
                 paste(sQuote(deparse(substitute(
                      arg
                 ))), "must be one of"),
                 paste(dQuote(choices), collapse = ", "))

     expre = substitute(
          match.arg(arg, choices, several.ok),
          list(
               arg = arg,
               choices = choices,
               several.ok = several.ok,
               call = call,
               msg = msg
          )
     )

     tryCatch(
          eval.parent(expre),
          error = function(e) {
               stop(errorCondition(msg, call = call))
          }
     )

}


#=================================================================
assert = function(cond,
                  expr,
                  type = "error",
                  pos = -2,
                  ...) {
     type = match_arg(type, choices = c("error", "warning", "message"))

     varargs = list(...)
     call = if (!is.null(varargs$call.) &&
                !isTRUE(varargs$call.))
          NULL
     else
          sys.call(pos)

     # for some reason I can't create a custom call object and still set
     # options for the condition (eg, immediate =T, appendLF=T, etc)
     # but this tryCatch block solves the issue
     # See https://stackoverflow.com/questions/63856125/setting-immediate-to-true-with-condition-object-in-r
     if (!cond) {
          out = tryCatch(
               switch(
                    type,
                    "error" = stop(expr, ...),
                    "message" = message(expr, ...),
                    'warning' = warning(expr, ...)
               ),
               error = function(e) {
                    e$call <- call
                    stop(e)
               },
               message = function(m)
                    message(m),
               warning = function(w) {
                    w$call <- call
                    warning(w)
               }
          )
     }

     return(invisible(NULL))

}

# function to help with getting two offsets instead of one from model frame
get.offset = function(formula, model.frame) {
     name = grep("offset",
                 unlist(strsplit(as.character(formula), " \\+ ")),
                 value = T)

     if (length(name) == 0)
          return(rep(0, nrow(model.frame)))
     else
          return(as.vector(model.frame[, name]))
}

# takes vector of all estimated parameters and subets into coefmat
# to extract each part of each equation (eg. count, zi, eq1, eq2)
make.coefmat = function(beta, se, z, pval, index) {
     if (is.null(index))
          return(NULL)
     else
          out = do.call(cbind, lapply(list(beta, se, z, pval), "[", index))

     colnames(out) = c("Estimate", "Std. Err.", "Z value", "Pr(>|z|)")
     return(out)
}

## Environment for regression functions to skip arg checking in below functions
e.check = rlang::env()


# Muffle original NA warning, rewarn/stop with updated call and message from parent function. Can also suppress a pattern
alter.cond = function(expre,
                      suppress = F,
                      pattern = 'NA|NaN',
                      new.message = NULL,
                      new.call = match.call(sys.function(pos), sys.call(pos)),
                      pos = -1){

     # if(is.null(new.call)){
     #
     # }

     withCallingHandlers(
          eval.parent(substitute(expre, list(expre = expre))),
          warning = function(w) {

               w$call = new.call
               old.message = w$message
               if(grepl(pattern, old.message) && !is.null(new.message)){
                    w$message = new.message
               }

               if(suppress && grepl(pattern, old.message)){
                    invokeRestart("muffleWarning")
               }
               else {
                    warning(w)
                    invokeRestart("muffleWarning")
               }
          },
          error = function(e){

               e$call = new.call

               if(!is.null(new.message))
                    e$message = new.message

               stop(e)
          }
     )

}

# Function to suppress warnings about NA and NaN values in likelihood search
# but keep other warnings
ll.suppress = function(warn) {
     if (grepl("NA|NaN", warn))
          invokeRestart("muffleWarning")
}


rep_len.custom = function(n, var){
     if(is.null(var))
          return(rep(.5, n))
     else if(length(var) > 1)
          rep(var, length.out = n)
     else
          return(var)
}


# function for setting defaults to each optimizer
set.defaults = function(de.list, de.names, de.values) {
     for (i in seq_along(de.names)) {
          if (is.null(de.list[[de.names[i]]]))
               de.list[[de.names[i]]] = de.values[i]
     }
     return(de.list)
}





scaler = function(df, scaling = 1){

  to_scale = !(sapply(df, function(col) length(unique(col)) %in% c(1,2) || is.character(col) || is.factor(col))
               | grepl("^offset\\(|^\\(weights\\)$", names(df)))

  to_scale[c(1,2)] = c(F,F) # first two are always outcome vectors, never scale them.

  if(sum(to_scale) == 0)
       return(df)

  if(scaling == 3) {
    mn = apply(df[, to_scale, drop=F], 2, min)
    mx = apply(df[, to_scale, drop=F], 2, max)
  }

  df[,to_scale] = switch(
    scaling,
     scale(df[,to_scale], center = T, scale = T),
     scale(df[,to_scale], center = T, scale = T)/2,
     scale(df[,to_scale], center = mn, scale = (mx - mn))
  )


  attr(df, "scaled") = names(to_scale)[to_scale]
  return(df)
}




vprob = function(x, eq = 'both'){
     eq = match_arg(eq, choices = c("both", "left", "right"))

     out = if(!is.null(x)){
          switch(eq,
                 "both" = x <= 1 & x >= 0,
                 "left" = x < 1 & x >= 0,
                 "right" = x <= 1 & x > 0)
     }
     else{
          TRUE
     }


     return(out)

}
