#### Define mpfr methods for Summary  group functions
####			     =======

### "Math" are done in ./Math.R , "Ops", "Arith", "Logic", "Compare" in ./Arith.R

###--> "max"   "min"   "range" "prod"  "sum"   "any"   "all"
.Summary.codes <-
    c("max" = 1, "min" = 2, "range" = 3, "prod" = 4, "sum" = 5,
      "any" = 6, "all" = 7)
storage.mode(.Summary.codes) <- "integer"

setMethod("Summary", "mpfr",
	  function(x, ..., na.rm=FALSE) {
	      iop <- .Summary.codes[.Generic]
	      r <- .Call(Summary_mpfr, if(length(x)) c(x, ...) else x, na.rm, iop)
	      if(iop <= 5)
		  new("mpfr", r)
	      else ## any, all :
		  r
	  })

## "mean": based on sum() :
setMethod("mean", "mpfr", function(x, ...) sum(x, ...)/length(x))

## FIXME: can do this considerably faster in C:
setMethod("which.max", "mpfr", function(x) which.max(x == max(x)))
setMethod("which.min", "mpfr", function(x) which.max(x == min(x)))
