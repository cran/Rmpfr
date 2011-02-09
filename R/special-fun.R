## erf(), erfc()

erf <- function(x) {
    if(is.numeric(x)) 2 * pnorm(x * sqrt(2)) - 1
    else if(is(x, "mpfr")) { # maybe also mpfrMatrix
        ##new("mpfr", .Call("Math_mpfr", x, .Math.codes[["erf"]], PACKAGE="Rmpfr"))
        x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["erf"]], PACKAGE="Rmpfr")
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
        x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["erfc"]], PACKAGE="Rmpfr")
        x
    }
    else stop("invalid class(x): ", class(x))
}
##    pnorm(x* sqrt(2), lower=FALSE) = erfc(x))/2
##==> pnorm(x., lower=TRUE)  = erfc(x./sqrt(2))/2


## zeta()
zeta <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["zeta"]], PACKAGE="Rmpfr")
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
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["Eint"]], PACKAGE="Rmpfr")
    x
}

## Li_2() the dilogarithm
Li2 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["Li2"]], PACKAGE="Rmpfr")
    x
}


### ------------- Bessel: ---------
## j0, j1, jn
## y0, y1, yn
j0 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["j0"]], PACKAGE="Rmpfr")
    x
}
j1 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["j1"]], PACKAGE="Rmpfr")
    x
}
y0 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["y0"]], PACKAGE="Rmpfr")
    x
}
y1 <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["y1"]], PACKAGE="Rmpfr")
    x
}

Ai <- function(x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("Math_mpfr", x, .Math.codes[["Ai"]], PACKAGE="Rmpfr")
    x
}

jn <- function(n, x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("R_mpfr_jn", x, n, PACKAGE="Rmpfr")
    x
}
yn <- function(n, x) {
    x <- as(x, "mpfr")
    x@.Data[] <- .Call("R_mpfr_yn", x, n, PACKAGE="Rmpfr")
    x
}


###-------- 2-argument cases -------

## We want to automatically construct the methods needed:
## But atan2() as argument list and  signature  (y, x)
## where  beta() has  (a,b)

mpfrMath2setMeth.y.x <- function(fname, Csub) {
    stopifnot(existsFunction(fname),
	      is.character(fname), length(fname) == 1,
	      is.character(Csub ), length(Csub ) == 1)

    setMethod(fname, signature(y = "mpfr", x = "mpfr"),
	      function(y, x)
	      new("mpfr", .Call(Csub, y, x, PACKAGE="Rmpfr")))
    setMethod(fname, signature(y = "mpfr", x = "ANY"),
	      function(y, x)
	      new("mpfr", .Call(Csub, y, as(x, "mpfr"), PACKAGE="Rmpfr")))
    setMethod(fname, signature(y = "ANY", x = "mpfr"),
	      function(y, x)
	      new("mpfr", .Call(Csub, as(y, "mpfr"), x, PACKAGE="Rmpfr")))


    setMethod(fname, signature(y = "mpfrArray", x = "mpfrArray"),
              function(y, x) {
                  if(dim(x) != dim(y))
                      stop("array dimensions differ")
                  x@.Data[] <- .Call(Csub, y, x, PACKAGE="Rmpfr")
                  x
              })
    setMethod(fname, signature(y = "mpfrArray", x = "ANY"),
              function(y, x) {
                  if(length(y) %% length(x) != 0)
                      stop("length of first argument (array) is not multiple of the second argument's one")
                  y@.Data[] <- .Call(Csub, y, as(x, "mpfr"), PACKAGE="Rmpfr")
                  y
              })
    setMethod(fname, signature(y = "ANY", x = "mpfrArray"),
              function(y, x) {
                  if(length(x) %% length(y) != 0)
                      stop("length of second argument (array) is not multiple of the first argument's one")
                  x@.Data[] <- .Call(Csub, as(y, "mpfr"), x, PACKAGE="Rmpfr")
                  x
              })

} ## end{mpfrMath2setMeth.y.x}

## atan2():
mpfrMath2setMeth.y.x("atan2", "R_mpfr_atan2")

mpfrMath2setMeth.a.b <- function(fname, Csub) {
    stopifnot(existsFunction(fname),
	      is.character(fname), length(fname) == 1,
	      is.character(Csub ), length(Csub ) == 1)

    setMethod(fname, signature(a = "mpfr", b = "mpfr"),
	      function(a, b)
	      new("mpfr", .Call(Csub, a, b, PACKAGE="Rmpfr")))
    setMethod(fname, signature(a = "mpfr", b = "ANY"),
	      function(a, b)
	      new("mpfr", .Call(Csub, a, as(b, "mpfr"), PACKAGE="Rmpfr")))
    setMethod(fname, signature(a = "ANY", b = "mpfr"),
	      function(a, b)
	      new("mpfr", .Call(Csub, as(a, "mpfr"), b, PACKAGE="Rmpfr")))


    setMethod(fname, signature(a = "mpfrArray", b = "mpfrArray"),
              function(a, b) {
                  if(dim(b) != dim(a))
                      stop("array dimensions differ")
                  b@.Data[] <- .Call(Csub, a, b, PACKAGE="Rmpfr")
                  b
              })
    setMethod(fname, signature(a = "mpfrArray", b = "ANY"),
              function(a, b) {
                  if(length(a) %% length(b) != 0)
                      stop("length of first argument (array) is not multiple of the second argument's one")
                  a@.Data[] <- .Call(Csub, a, as(b, "mpfr"), PACKAGE="Rmpfr")
                  a
              })
    setMethod(fname, signature(a = "ANY", b = "mpfrArray"),
              function(a, b) {
                  if(length(b) %% length(a) != 0)
                      stop("length of second argument (array) is not multiple of the first argument's one")
                  b@.Data[] <- .Call(Csub, as(a, "mpfr"), b, PACKAGE="Rmpfr")
                  b
              })

} ## end{mpfrMath2setMeth.a.b}

mpfrMath2setMeth.a.b("beta",  "R_mpfr_beta")
mpfrMath2setMeth.a.b("lbeta", "R_mpfr_lbeta")

rm(mpfrMath2setMeth.y.x,
   mpfrMath2setMeth.a.b)


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
