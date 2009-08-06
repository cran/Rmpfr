## erf(), erfc()

erf <- function(x) {
    if(is.numeric(x)) 2 * pnorm(x * sqrt(2)) - 1
    else if(is(x, "mpfr")) { # maybe also mpfrMatrix
        ##new("mpfr", .Call("Math_mpfr", x, .Math.codes["erf"], PACKAGE="Rmpfr"))
        x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["erf"], PACKAGE="Rmpfr")
        x
    }
    else stop("invalid class(x): ", class(x))
}
##    pnorm(x* sqrt(2)) = (1 + erf(x))/2
##==> pnorm(x.)  = (1 + erf(x./sqrt(2)))/2

pnorm <- function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
{
    if(is.numeric(q) && is.numeric(mean) && is.numeric(sd))
        .Internal(pnorm(q, mean, sd, lower.tail, log.p))
    else if(is(q, "mpfr") || is(mean, "mpfr") || is(sd, "mpfr")) {
        stopifnot(length(lower.tail) == 1, length(log.p) == 1)
        q <- as(q, "mpfr")
        prec.q <- max(sapply(q, slot, "prec"))
        rt2 <- sqrt(mpfr(2, prec.q))
        if(lower.tail) {
            if(log.p && all(mean == 0))
                ## log(sd/2 * (1 + erf(q/rt2)))
                log(sd/2) + log1p(erf(q/rt2))
            else {
                r <- mean + sd/2 * (1 + erf(q/rt2))
                if(log.p) log(r) else r
            }
        } else { ## upper.tail
            if(log.p && all(mean == 0))
                ## log(sd/2 * erfc(q/rt2))
                log(sd/2) + log(erfc(q/rt2))
            else {
                r <- mean + sd/2 * erfc(q/rt2)
                if(log.p) log(r) else r
            }
        }

    } else stop("invalid arguments (q,mean,sd)")
}


erfc <- function(x) {
    if(is.numeric(x)) 2 * pnorm(x * sqrt(2), lower = FALSE)
    else if(is(x, "mpfr")) {
        x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["erfc"], PACKAGE="Rmpfr")
        x
    }
    else stop("invalid class(x): ", class(x))
}
##    pnorm(x* sqrt(2), lower=FALSE) = erfc(x))/2
##==> pnorm(x., lower=TRUE)  = erfc(x./sqrt(2))/2


## zeta()
zeta <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["zeta"], PACKAGE="Rmpfr")
    x
}

Bernoulli <- function(k, precBits = 128) {
    ## Purpose: Bernoulli Numbers (in high precision)
    ## -----------------------------------------------------------
    ## Arguments: k: positive integer vector
    ## -----------------------------------------------------------
    ## Author: Martin Maechler, Date: 12 Dec 2008, 11:35
    stopifnot(all(k > 0), k == as.integer(k))
    - k * zeta(if(is(k, "mpfr")) 1 - k else mpfr(1 - k, precBits=precBits))
}


## eint() "Exponential integral"
Ei <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["Eint"], PACKAGE="Rmpfr")
    x
}

### ------------- Bessel: ---------
## j0, j1, jn
## y0, y1, yn
j0 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["j0"], PACKAGE="Rmpfr")
    x
}
j1 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["j1"], PACKAGE="Rmpfr")
    x
}
y0 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["y0"], PACKAGE="Rmpfr")
    x
}
y1 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes["y1"], PACKAGE="Rmpfr")
    x
}

jn <- function(n, x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("R_mpfr_jn", x, as.integer(n), PACKAGE="Rmpfr")
    x
}
yn <- function(n, x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("R_mpfr_yn", x, as.integer(n), PACKAGE="Rmpfr")
    x
}

###-------- 2-argument cases -------

## atan2()
setMethod("atan2", signature(y = "mpfr", x = "mpfr"),
	  function(y, x)
	      new("mpfr", .Call("R_mpfr_atan2", y, x, PACKAGE="Rmpfr")))

setMethod("atan2", signature(y = "mpfr", x = "ANY"),
	  function(y, x) {
	      new("mpfr",
		  .Call("R_mpfr_atan2", y, as(x, "mpfr"), PACKAGE="Rmpfr"))
	  })
setMethod("atan2", signature(y = "ANY", x = "mpfr"),
	  function(y, x) {
	      new("mpfr",
		  .Call("R_mpfr_atan2", as(y, "mpfr"), x, PACKAGE="Rmpfr"))
	  })

setMethod("atan2", signature(y = "mpfrArray", x = "mpfrArray"),
	  function(y, x) {
	      if(dim(x) != dim(y))
		  stop("array dimensions differ")
	      x@.Data[] <- .Call("R_mpfr_atan2", y, x, PACKAGE="Rmpfr")
	      x
	  })

setMethod("atan2", signature(y = "mpfrArray", x = "ANY"),
	  function(y, x) {
	      if(length(y) %% length(x) != 0)
		  stop("length of first argument (array) is not multiple of the second argument's one")
	      y@.Data[] <- .Call("R_mpfr_atan2", y, as(x, "mpfr"),
				 PACKAGE="Rmpfr")
	      y
	  })

setMethod("atan2", signature(y = "ANY", x = "mpfrArray"),
	  function(y, x) {
	      if(length(x) %% length(y) != 0)
		  stop("length of second argument (array) is not multiple of the first argument's one")
	      x@.Data[] <- .Call("R_mpfr_atan2", as(y, "mpfr"), x,
				 PACKAGE="Rmpfr")
	      x
	  })

## hypot()
hypot <- function(x,y) {
    if(is(x, "mpfrArray") || is.array(x)) {
	if(is.array(x)) x <- mpfrArray(x, 128L, dim=dim(x), dimnames(x))
	if(is.array(y)) y <- mpfrArray(y, 128L, dim=dim(y), dimnames(y))
	if(is(y, "mpfrArray")) {
	    if(dim(x) != dim(y))
		stop("array dimensions differ")
	    x@.Data[] <- .Call("R_mpfr_hypot", x, y, PACKAGE="Rmpfr")
	    x
	} else { ## y is not (mpfr)Array
	    if(length(x) %% length(y) != 0)
		stop("length of first argument (array) is not multiple of the second argument's one")
	    x@.Data[] <- .Call("R_mpfr_hypot", x, as(y, "mpfr"), PACKAGE="Rmpfr")
	    x
	}
    } else if(is(y, "mpfrArray")) {
	if(length(y) %% length(x) != 0)
	    stop("length of second argument (array) is not multiple of the first argument's one")
	y@.Data[] <- .Call("R_mpfr_hypot", as(x, "mpfr"), y, PACKAGE="Rmpfr")
	y
    }
    else
	new("mpfr", .Call("R_mpfr_hypot", as(x, "mpfr"), as(y, "mpfr"),
			  PACKAGE="Rmpfr"))
}
