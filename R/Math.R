#### Define mpfr methods for Math  and Math2  group functions
####                        ======     =====

### "Arith", "Compare",..., are in ./Arith.R
###  ----                            ~~~~~~~

## [1] "abs"    "sign"    "sqrt"    "ceiling" "floor" "trunc" "cummax"
## [8] "cummin" "cumprod" "cumsum"  "exp"     "expm1" "log"   "log10"
##[15] "log2"   "log1p"   "cos"     "cosh"    "sin"   "sinh"  "tan"
##[22] "tanh"   "acos"    "acosh"   "asin"    "asinh" "atan"  "atanh"
##[29] "gamma"  "lgamma"  "digamma" "trigamma"

if(FALSE) ## here are the individual function
    dput(getGroupMembers("Math"))

## Uniform interface to C:
##
## Pass integer code to call and do the rest in C
## Codes from ~/R/D/r-devel/R/src/main/names.c :
.Math.codes <-
    c(
      "floor" =     1,
      "ceiling" =   2,
      "sqrt" =      3,
      "sign" =      4,
      "exp" =      10,
      "expm1" =    11,
      "log1p" =    12,
      "cos" =      20,
      "sin" =      21,
      "tan" =      22,
      "acos" =     23,
      "asin" =     24,
      "cosh" =     30,
      "sinh" =     31,
      "tanh" =     32,
      "acosh" =    33,
      "asinh" =    34,
      "atanh" =    35,
      "lgamma" =   40,
      "gamma" =    41,
      "digamma" =  42,
      "trigamma" = 43)

.Math.gen <- getGroupMembers("Math")

## Those "Math" group generics that are not in the do_math1 table above

.Math.codes <-
    c(.Math.codes,
      "trunc" = 0, "atan" = 25, # "abs" has own method!
      "log" = 13, "log2" = 14, "log10" = 15,
      "cummax" = 51, "cummin" = 52, "cumprod" = 53, "cumsum" = 54,
      ## These are *NOT* in R's  Math group, but 1-argument math functions
      ## available in the mpfr - library:
      "erf" = 101, "erfc" = 102, "zeta" = 104, "Eint" = 106, "Li2" = 107,
      "j0" = 111, "j1" = 112, "y0" = 113, "y1" = 114)
storage.mode(.Math.codes) <- "integer"

if(FALSE)
.Math.gen[!(.Math.gen %in% names(.Math.codes))]
## "abs" -- only one left

## A few ones have a very simple method:
setMethod("sign", "mpfr",
	  function(x) sapply(x, function(e) e@sign))

setMethod("abs", "mpfr",
	  function(x) {
	      for(i in seq_along(x)) x[[i]]@sign <- 1L
	      x
	  })

## Note that  factorial() and lfactorial() automagically work through  [l]gamma()
## but for the sake of "exact for integer"
setMethod("factorial", "mpfr",
	  function(x) {
	      r <- gamma(x + 1)
	      if(mpfr.is.integer(x)) round(r) else r
	  })

## "log" is still special with its 'base' :
setMethod("log", signature(x = "mpfr"),
	  function(x, base) {
	      if(!missing(base) && base != exp(1))
		  stop("base != exp(1) is not yet implemented")
	      x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["log"],
				 PACKAGE="Rmpfr")
	      x
	  })

setMethod("Math", signature(x = "mpfr"),
	  function(x) {
	      x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[.Generic],
			   PACKAGE="Rmpfr")
	      x
	  })

setMethod("Math2", signature(x = "mpfr"),
	  function(x, digits) {
	      ## NOTA BENE: vectorized in  'x'
	      if(any(ret.x <- !is.finite(x) | mpfr.is.0(x))) {
		  if(any(ok <- !ret.x))
		      x[ok] <- callGeneric(x[ok], digits=digits)
		  return(x)
	      }
              if(!missing(digits)) {
                  digits <- as.integer(round(digits))
                  if(is.na(digits)) return(x + digits)
              } ## else: default *depends* on the generic

	      ## now: both x and digits are finite
	      pow10 <- function(d) mpfr(rep.int(10., length(d)),
					precBits = log2(10)*as.numeric(d))^ d
	      rint <- function(x) { ## have x >= 0 here
		  sml.x <- (x < .Machine$integer.max)
		  r <- x
		  if(any(sml.x)) {
		      x.5 <- x[sml.x] + 0.5
		      ix <- as.integer(x.5)
		      ## implement "round to even" :
		      if(any(doDec <- (abs(x.5 - ix) < 10*.Machine$double.eps & (ix %% 2))))
			  ix[doDec] <- ix[doDec] - 1L
		      r[sml.x] <- ix
		  }
		  if(!all(sml.x)) { ## large x - no longer care for round to even
		      r[!sml.x] <- floor(x[!sml.x] + 0.5)
		  }
		  r
	      }
	      neg.x <- x < 0
	      x[neg.x] <- - x[neg.x]
	      sgn <- ifelse(neg.x, -1, +1)
	      switch(.Generic,
		     "round" = { ## following ~/R/D/r-devel/R/src/nmath/fround.c :
			 if(missing(digits) || digits == 0)
			     sgn * rint(x)
			 else if(digits > 0) {
			     p10 <- pow10(digits)
			     intx <- floor(x)
			     sgn * (intx + rint((x-intx) * p10) / p10)
			 }
			 else { ## digits < 0
			     p10 <- pow10(-digits)
			     sgn * rint(x/p10) * p10
			 }
		     },
		     "signif" = { ## following ~/R/D/r-devel/R/src/nmath/fprec.c :
                         if(missing(digits)) digits <- 6L
			 if(digits > max(getPrec(x)) * log10(2))
			     return(x)
			 if(digits < 1) digits <- 1L
			 l10 <- log10(x)
			 e10 <- digits - 1L - floor(l10)
			 r <- x
			 pos.e <- (e10 > 0) ##* 10 ^ e, with e >= 1 : exactly representable
			 if(any(pos.e)) {
			     p10 <- pow10(e10[pos.e])
			     r[pos.e] <- sgn[pos.e]* rint(x[pos.e]*p10) / p10
			 }
			 if(any(neg.e <- !pos.e)) {
			     p10 <- pow10(-e10[neg.e])
			     r[neg.e] <- sgn[neg.e]* rint(x[neg.e]/p10) * p10
			 }
			 r
		     },
		     stop(gettextf("Non-Math2 group generic '%s' -- should not happen",
				   .Generic)))
	  })

##---- mpfrArray / mpfrMatrix --- methods -----------------

## not many needed: "mpfrArray" contain "mpfr",
## i.e., if the above methods are written "general enough", they apply directly

setMethod("sign", "mpfrArray",
	  function(x) structure(sapply(x, function(e) e@sign),
				dim = dim(x),
				dimnames = dimnames(x)))
