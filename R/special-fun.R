## erf(), erfc()

erf <- function(x) {
    if(is.numeric(x)) 2 * pnorm(x * sqrt(2)) - 1
    else if(is.mpfr(x)) { # maybe also mpfrMatrix
	##new("mpfr", .Call(Math_mpfr, x, .Math.codes[["erf"]]))
	x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["erf"]])
	x
    }
    else stop("invalid class(x): ", class(x))
}
##    pnorm(x* sqrt(2)) = (1 + erf(x))/2
##==> pnorm(x.)	 = (1 + erf(x./sqrt(2)))/2

##    pnorm(x* sqrt(2), lower=FALSE) = erfc(x)/2
##==> pnorm(x., lower=TRUE)  = erfc(x./sqrt(2))/2
erfc <- function(x) {
    if(is.numeric(x)) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
    else if(is.mpfr(x)) {
	x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["erfc"]])
	x
    }
    else stop("invalid class(x): ", class(x))
}

pnorm <- function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
{
    if(is.numeric(q) && is.numeric(mean) && is.numeric(sd))
	stats__pnorm(q, mean, sd, lower.tail=lower.tail, log.p=log.p)
    else if((q.mp <- is.mpfr(q)) || is.mpfr(mean) || is.mpfr(sd)) {
	stopifnot(length(lower.tail) == 1L, length(log.p) == 1L)
	rr <- q <- ((if(q.mp) q else as(q, "mpfr")) - mean) / sd
	if(any(neg <- (q < 0))) ## swap those:	Phi(-z) = 1 - Phi(z)
	    rr[neg] <- pnorm(-q[neg], lower.tail = !lower.tail, log.p=log.p)
	if(any(pos <- !neg)) {
	    q <- q[pos] #==> now  q >= 0
	    prec.q <- max(.getPrec(q))
	    two <- mpfr(2, prec.q)
	    rt2 <- sqrt(two)
	    rr[pos] <- if(lower.tail) {
		if(log.p) {
		    r <- q
		    sml <- q < 0.67448975
		    if(any(sml)) {
			eq2 <- erf(q[sml]/rt2) ## |eq2| < 1/2 <==> |q/rt2| < 0.47693627620447
			##                        <==>  sml   <==>   |q|   < 0.67448975019608
			r[ sml] <- log1p(eq2) - log(two)
		    }
		    if(any(!sml)) {
			ec2 <- erfc(q[!sml]/rt2) ## ==> ec2 = 1-eq2 <= 1 - 1/2 = 1/2
			r[!sml] <- log1p(-0.5*ec2)
		    }
		    r
		}
		else ## !log.p
		    (1 + erf(q/rt2))/2
	    } else { ## upper.tail
		r <- erfc(q/rt2) / 2
		if(log.p) log(r) else r
	    }
	}
	rr
    } else stop("(q,mean,sd) must be numeric or \"mpfr\"")
}#{pnorm}

dnorm <- function (x, mean = 0, sd = 1, log = FALSE) {
    if(is.numeric(x) && is.numeric(mean) && is.numeric(sd))
	stats__dnorm(x, mean, sd, log=log)
    else if((x.mp <- is.mpfr(x)) || is.mpfr(mean) || is.mpfr(sd)) {
	## stopifnot(length(log) == 1)
	prec <- pmax(53, getPrec(x), getPrec(mean), getPrec(sd))
	if(!x.mp) x <- mpfr(x, prec)
	x <- (x - mean) / sd
	twopi <- 2*Const("pi", prec)
	## f(x) =  1/(sigma*sqrt(2pi)) * exp(-1/2 x^2)
	if(log) ## log( f(x) ) = -[ log(sigma) + log(2pi)/2 + x^2 / 2]
	    -(log(if(is.mpfr(sd)) sd else mpfr(sd, prec))  + (log(twopi) + x*x)/2)
	else exp(-x^2/2) / (sd*sqrt(twopi))
    } else stop("invalid arguments (x,mean,sd)")
}

dpois <- function (x, lambda, log = FALSE,
                   useLog = { ## MPFR overflow:
                       ln2 <- log(2)
                       any(lambda >= -.mpfr_erange("Emin")*ln2) ||
                           any(x*log(lambda) >= .mpfr_erange("Emax")*ln2)
                   })
{
    if(is.numeric(x) && is.numeric(lambda)) ## standard R
	stats__dpois(x, lambda, log=log)
    else if((l.mp <- is.mpfr(lambda)) | (x.mp <- is.mpfr(x))) {
	prec <- pmax(53, getPrec(lambda), getPrec(x))
	if(!l.mp) lambda <- mpfr(lambda, prec)
	if(!x.mp) x <- mpfr(x, prec)
        if(log || useLog) {
            r <-  -lambda  + x*log(lambda) - lfactorial(x)
            if(log) r else exp(r)
        }
	else exp(-lambda) * lambda^x /  factorial(x)
    } else
	stop("(x,lambda) must be numeric or \"mpfr\"")
}

dbinom <- function(x, size, prob, log = FALSE,
                   useLog = any(abs(x) > 1e6) || ## MPFR overflow:
                       max(x*log(prob), (size-x)*log1p(-prob)) >=
                         .mpfr_erange("Emax")*log(2))
{
    if(is.numeric(x) && is.numeric(size) && is.numeric(prob)) ## standard R
	stats__dbinom(x, size, prob, log=log)
    else if((s.mp <- is.mpfr(size)) |
	    (p.mp <- is.mpfr(prob)) || is.mpfr(x)) {
	stopifnot(is.whole(x)) # R's dbinom() gives NaN's with a warning..
        if(!is.integer(x)) x <- as.integer(x) # chooseMpfr() needs
	prec <- pmax(53, getPrec(size), getPrec(prob), getPrec(x))
	if(!s.mp) size <- mpfr(size, prec)
	if(!p.mp) prob <- mpfr(prob, prec)
	## n:= size, p:= prob,	compute	 P(x) = choose(n, x) p^x (1-p)^(n-x)
        if(useLog) { # do *not* use chooseMpfr() {which is O(x^2)}
            lC.nx <- ## lchoose(size, x), but that is *not* yet available for "mpfr" __FIXME?__
                lfactorial(size) - (lfactorial(x) + lfactorial(size-x))
        } else {
            C.nx <- chooseMpfr(size, x)
            lC.nx <- log(C.nx)
        }
        if(log || useLog) {
            r <- lC.nx + x*log(prob) + (size-x)*log1p(-prob)
            if(log) r else exp(r)
        }
	else C.nx * prob^x * (1-prob)^(size-x)
    } else
	stop("(x,size, prob) must be numeric or \"mpfr\"")
}## {dbinom}

dnbinom <- function (x, size, prob, mu, log = FALSE, useLog = any(x > 1e6)) {
    if(!missing(mu)) {
        if (!missing(prob))
            stop("'prob' and 'mu' both specified")
        ## Using argument 'mu' instead of 'prob'
        if (size == Inf)
            return( dpois(x, lambda=mu, log) )
        else
            prob <- size/(size+mu)  #  and continue :
    }
    if(is.numeric(x) && is.numeric(size) && is.numeric(prob)) { ## standard R
        if(!missing(mu))
            stats__dnbinom(x, size, mu=mu,     log=log)
        else
            stats__dnbinom(x, size, prob=prob, log=log)
    } else if((s.mp <- is.mpfr(size)) |
	    (p.mp <- is.mpfr(prob)) || is.mpfr(x)) {
	stopifnot(is.whole(x)) # R's dbinom() gives NaN's with a warning..
        if(!is.integer(x) && !useLog)
            x <- as.integer(x) # chooseMpfr() needs it
	prec <- pmax(53, getPrec(size), getPrec(prob), getPrec(x))
	if(!s.mp) size <- mpfr(size, prec)
	if(!p.mp) prob <- mpfr(prob, prec)
	## n:= size, p:= prob,	compute	 P(x) = choose(n+x-1, x) * p^n * (1-p)^x
        if(!useLog) {
            C.nx <- chooseMpfr(size+x-1, x)
            if(log || ## MPFR overflow:
               max(size*log(prob), x*log1p(-prob)) >= .mpfr_erange("Emax")*log(2))
            {
                r <- log(C.nx) + size*log(prob) + x*log1p(-prob)
                if(log) r else exp(r)
            }
            else C.nx * prob^size * (1-prob)^x
        } else { # x not integer, typically  |x| > .Machine$integer.max (= 2147483647 = 2^31 - 1)
            ## => x is large but  size >= x is even larger ... so everything is large
            ## FIXME (?) suffering from cancellation (when ?) !
            logC.nx <- lgamma(size+x) - lgamma(size) - lgamma(x+1)
            if(log)  logC.nx + size*log(prob) + x*log1p(-prob)
            else exp(logC.nx + size*log(prob) + x*log1p(-prob))
        }
    } else
	stop("(x,size, prob | mu) must be numeric or \"mpfr\"")
}## {dnbinom}


dgamma <- function(x, shape, rate = 1, scale = 1/rate, log = FALSE) {
    missR <- missing(rate)
    missS <- missing(scale)
    if (!missR && !missS) { ## as stats::dgamma()
        if (abs(rate * scale - 1) < 1e-15)
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
    ## and now use 'scale' only
    if(is.numeric(x) && is.numeric(shape) && is.numeric(scale))
        stats__dgamma(x, shape, scale=scale, log=log)
    else if((sh.mp <- is.mpfr(shape)) |
	    (sc.mp <- is.mpfr(scale)) || is.mpfr(x)) {
        ##     f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)  ; a=shape, s=scale
        ## log f(x) = -a*log(s) - lgamma(a) + (a-1)*log(x) - (x/s)
	if(!sh.mp || !sc.mp) {
	    prec <- pmax(53, getPrec(shape), getPrec(scale), getPrec(x))
	    if(!sh.mp)
		shape <- mpfr(shape, prec)
	    else ## !sc.mp :
		scale <- mpfr(scale, prec)
	}
	## for now, "cheap", relying on "mpfr" arithmetic to be smart
	## "TODO":  Use C.Loader's formulae via dpois_raw() ,  bd0() etc
	## lgam.sh <- lgamma(shape)
	## ldgamma <- function(x, shp, s) -shp*log(s) -lgam.sh + (shp-1)*log(x) - (x/s)
        ldgamma <- function(x, shp, s) -shp*log(s) -lgamma(shp) + (shp-1)*log(x) - (x/s)
	if(log)
	    ldgamma(x, shape, scale)
	else { ## use direct [non - log-scale] formula when applicable
	    ## ok <- lgam.sh < log(2) * Rmpfr:::.mpfr.erange("Emax") & ## finite gamma(shape) := exp(lgam.sh)
	    ##       is.finite(xsh1 <- x^(shape-1)) &
	    ##       !is.na(r <- xsh1 * exp(-(x/scale)) / (scale^shape * exp(lgam.sh)))
	    ## r[!ok] <- exp(ldgamma(x[!ok], shape[!ok], scale[!ok]))
	    ## r
	    exp(ldgamma(x, shape, scale))
	}
    } else
	stop("(x, shape, scale) must be numeric or \"mpfr\"")
}## {dgamma}



## zeta()
zeta <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr") # keep "mpfrArray"
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["zeta"]])
    x
}

## "FIXME" -- rather use 'bigq' in gmp and the "sumBin" algorithm from copula!
Bernoulli <- function(k, precBits = 128) {
    ## Purpose: Bernoulli Numbers (in high precision)
    ## -----------------------------------------------------------
    ## Arguments: k: non-negative integer vector
    ## -----------------------------------------------------------
    ## Author: Martin Maechler, Date: 12 Dec 2008, 11:35
    stopifnot(all(k >= 0), k == as.integer(k))
    r <- - k * zeta(if(is.mpfr(k)) 1 - k else mpfr(1 - k, precBits=precBits))
    if(any(k0 <- k == 0)) r[k0] <- mpfr(1, precBits=precBits)
    r
}


## eint() "Exponential integral"
Ei <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr") # keep "mpfrArray"
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["Eint"]])
    x
}

## Li_2() the dilogarithm
Li2 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr") # keep "mpfrArray"
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["Li2"]])
    x
}


### ------------- Bessel: ---------
## j0, j1, jn
## y0, y1, yn
j0 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr") # keep "mpfrArray"
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["j0"]])
    x
}
j1 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["j1"]])
    x
}
y0 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["y0"]])
    x
}
y1 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["y1"]])
    x
}

Ai <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["Ai"]])
    x
}

jn <- function(n, x, rnd.mode = c('N','D','U','Z','A')) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(R_mpfr_jn, x, n, match.arg(rnd.mode))
    x
}
yn <- function(n, x, rnd.mode = c('N','D','U','Z','A')) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(R_mpfr_yn, x, n, match.arg(rnd.mode))
    x
}


###-------- 2-argument cases -------

## We want to automatically construct the methods needed:
## But atan2() as argument list and  signature	(y, x)
## where  beta() and lbeta()  have  (a,b) --> cannot treat them identically;
## and treat atan2() speparately

## NB: atan2(), beta() and lbeta() all have implicitGeneric()s in methods with no '...'
## ==  ---> canNOT have 3rd argument : rnd.mode = c('N','D','U','Z','A')
##     ---> using   "N"  instead of    match.arg(rnd.mode)
setMethod("atan2", signature(y = "mpfr", x = "mpfr"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, y, x, "N")))
setMethod("atan2", signature(y = "mpfr", x = "numeric"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, y, .mpfr(x, 128L), "N")))
setMethod("atan2", signature(y = "numeric", x = "mpfr"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, .mpfr(y, 128L), x, "N")))
setMethod("atan2", signature(y = "mpfr", x = "ANY"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, y, as(x, "mpfr"), "N")))
setMethod("atan2", signature(y = "ANY", x = "mpfr"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, as(y, "mpfr"), x, "N")))

setMethod("atan2", signature(y = "mpfrArray", x = "mpfrArray"),
          function(y, x) {
              if(dim(x) != dim(y))
                  stop("array dimensions differ")
              x@.Data[] <- .Call(R_mpfr_atan2, y, x, "N")
              x
          })
setMethod("atan2", signature(y = "mpfrArray", x = "ANY"),
          function(y, x) {
              if(length(y) %% length(x) != 0)
                  stop("length of first argument (array) is not multiple of the second argument's one")
              y@.Data[] <- .Call(R_mpfr_atan2, y, if(is.numeric(x))
                  .mpfr(x, 128L) else as(x, "mpfr"), "N")
              y
          })
setMethod("atan2", signature(y = "ANY", x = "mpfrArray"),
          function(y, x) {
              if(length(x) %% length(y) != 0)
                  stop("length of second argument (array) is not multiple of the first argument's one")
              x@.Data[] <- .Call(R_mpfr_atan2, if(is.numeric(y))
                  .mpfr(y, 128L) else as(y, "mpfr"), x, "N")
              x
          })

## Using  "macro"  {instead of previous aux. function  mpfrMath2setMeth.a.b() :
for (ff in list(c("beta",  "R_mpfr_beta"),
                c("lbeta", "R_mpfr_lbeta"))) eval(substitute(
{
    setMethod(fname, signature(a = "mpfr", b = "mpfr"),
	      function(a, b) new("mpfr", .Call(Csub, a, b, "N")))
    setMethod(fname, signature(a = "mpfr", b = "numeric"),
	      function(a, b) new("mpfr", .Call(Csub, a, .mpfr(b, 128L), "N")))
    setMethod(fname, signature(a = "numeric", b = "mpfr"),
	      function(a, b) new("mpfr", .Call(Csub, .mpfr(a, 128L), b, "N")))
    setMethod(fname, signature(a = "mpfr", b = "ANY"),
	      function(a, b) new("mpfr", .Call(Csub, a, as(b, "mpfr"), "N")))
    setMethod(fname, signature(a = "ANY", b = "mpfr"),
	      function(a, b) new("mpfr", .Call(Csub, as(a, "mpfr"), b, "N")))

    setMethod(fname, signature(a = "mpfrArray", b = "mpfrArray"),
	      function(a, b) {
		  if(dim(b) != dim(a))
		      stop("array dimensions differ")
		  b@.Data[] <- .Call(Csub, a, b, "N")
		  b
	      })
    setMethod(fname, signature(a = "mpfrArray", b = "ANY"),
	      function(a, b) {
		  if(length(a) %% length(b) != 0)
		      stop("length of first argument (array) is not multiple of the second argument's one")
		  a@.Data[] <- .Call(Csub, a, if(is.numeric(b))
				     .mpfr(b, 128L) else as(b, "mpfr"), "N")
		  a
	      })
    setMethod(fname, signature(a = "ANY", b = "mpfrArray"),
	      function(a, b) {
		  if(length(b) %% length(a) != 0)
		      stop("length of second argument (array) is not multiple of the first argument's one")
		  b@.Data[] <- .Call(Csub, if(is.numeric(a))
				     .mpfr(a, 128L) else as(a, "mpfr"), b, "N")
		  b
	      })

}, list(fname = ff[[1]], Csub = as.symbol(ff[[2]]))))


## hypot()
hypot <- function(x,y, rnd.mode = c('N','D','U','Z','A')) {
    if(is(x, "mpfrArray") || is.array(x)) {
	if(is.array(x)) x <- mpfrArray(x, 128L, dim=dim(x), dimnames(x))
	if(is.array(y)) y <- mpfrArray(y, 128L, dim=dim(y), dimnames(y))
	if(is(y, "mpfrArray")) {
	    if(dim(x) != dim(y))
		stop("array dimensions differ")
	    x@.Data[] <- .Call(R_mpfr_hypot, x, y, match.arg(rnd.mode))
	    x
	} else { ## y is not (mpfr)Array
	    if(length(x) %% length(y) != 0)
		stop("length of first argument (array) is not multiple of the second argument's one")
	    x@.Data[] <- .Call(R_mpfr_hypot, x, as(y, "mpfr"), match.arg(rnd.mode))
	    x
	}
    } else if(is(y, "mpfrArray")) {
	if(length(y) %% length(x) != 0)
	    stop("length of second argument (array) is not multiple of the first argument's one")
	y@.Data[] <- .Call(R_mpfr_hypot, as(x, "mpfr"), y, match.arg(rnd.mode))
	y
    }
    else
	new("mpfr", .Call(R_mpfr_hypot, as(x, "mpfr"), as(y, "mpfr"), match.arg(rnd.mode)))
}

## The Beta(a,b)  Cumulative Probabilities are exactly computable for *integer* a,b:
pbetaI <- function(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
		   precBits = NULL,
                   useRational = !log.p && !is.mpfr(q) && is.null(precBits),
                   rnd.mode = c('N','D','U','Z','A'))
{
    stopifnot(length(shape1) == 1, length(shape2) == 1,
	      is.whole(shape1), is.whole(shape2),
	      shape1 >= 0, shape2 >= 0,
	      length(lower.tail) == 1, length(log.p) == 1,
	      0 <= q, q <= 1, ncp == 0,
	      is.null(precBits) ||
	      (is.numeric(precBits) && is.whole(precBits) && precBits >= 2))
    ## Care for too large (a,b) and "integer overflow".
    ## NB:  below have 0:(b - 1) or 0:(a - 1)
    max.ab <- 2^20
    if(is.na(a <- as.integer(shape1)) || (!lower.tail && a > max.ab))
        stop("a = shape1 is too large for 'lower.tail=FALSE' and the current algorithm")
    if(is.na(b <- as.integer(shape2)) || (lower.tail && b > max.ab))
        stop("b = shape2 is too large for 'lower.tail=TRUE' and the current algorithm")
    n <- a + b - 1L
    if(!useRational) {
        pr.x <- getPrec(q, bigq. = 256L)
        if(is.null(precBits)) {
            aq <- abs(as.numeric(q))
            mq <- if(any(po <- aq > 0)) min(aq[po]) else 1 # ==> log = 0
            ## -n*log(|x|): such that 1 - |x|^n does not completely cancel
            precBits <- max(128L, pr.x, -as.numeric(n)*log(mq))
        }
        if(pr.x < precBits || !is.mpfr(q))
            q <- mpfr(q, precBits=precBits)
        mpfr1 <- list(.Call(const_asMpfr, 1, 16L, "N")) # as prototype for vapply()
    }
    F <- if(log.p) log else identity
    ## FIXME: logspace add sum   lsum(.) should be more accurate for large n ==> could use larger a,b

    if(lower.tail) {
	## The prob. is	  P[ X <= x ] = \sum_{k=a}^ n    (n \\ k) x^k (1-x)^(n-k)
        ## but we want to sum from 0 {smallest --> largest} as well:
        ##                P[ X <= x ] = \sum_{k=0}^{b-1} (n \\ k) (1-x)^k x^(n-k)
	k <- 0:(b - 1L)
        FUN.x <- function(x) sum(n.choose.k * (1-x)^k * x^(n-k))
    } else { ## upper tail
	## Prob. P[ X > q ] =  1 - P[X <= q ] = \sum_{k=0}^{a-1} (n \\ k) x^k (1-x)^(n-k)
	k <- 0:(a - 1L)
        FUN.x <- function(x) sum(n.choose.k * x^k * (1-x)^(n-k))
    }
    n.choose.k <- chooseZ(n, k)
    if(useRational) {
        q <- as.bigq(q)
        if(length(q) == 1L) FUN.x(q) else c_bigq(lapply(q, FUN.x))
    } else { # mpfr
        roundMpfr(F(
            ## "vapply() for "mpfr"
            new("mpfr",
                vapply(q, FUN.x, mpfr1))),
            ## reduce the precision, in order to not "claim wrongly":
            precBits=precBits, match.arg(rnd.mode))
    }
}

### MPFR version >= 3.2.0 :
   "https://www.mpfr.org/mpfr-current/mpfr.html#index-mpfr_005fgamma_005finc"
##
## >>> Note: the current implementation of mpfr_gamma_inc is slow for large values of rop or op,
## >>> ====  in which case some internal overflow might also occur.
##
##  mpfr_gamma_inc(a,x) =: igamma(a,x)    where
##
## igamma(a,x) = "upper" incomplete gamma  Γ(a,x) :=: Γ(a) - γ(a,x);
##	γ(a,x) = "lower" incomplete gamma  γ(a,x) := ₀∫ˣ tᵃ⁻¹ e⁻ᵗ dt, and
## R's  pgamma(x, a) :==  γ(a,x) / Γ(a)
##
## >>> ../man/igamma.Rd <<<
igamma <- function(a,x, rnd.mode = c('N','D','U','Z','A')) {
    if(mpfrVersion() < "3.2.0")
	stop("igamma()'s MPFR equivalent needs mpfr version >= 3.2.0, but mpfrVersion()=",
	     mpfrVersion())
    if(is(a, "mpfrArray") || is.array(a)) {
	if(is.array(a)) a <- mpfrArray(a, 128L, dim=dim(a), dimnames(a))
	if(is.array(x)) x <- mpfrArray(x, 128L, dim=dim(x), dimnames(x))
	if(is(x, "mpfrArray")) {
	    if(dim(a) != dim(x))
		stop("array dimensions differ")
	    a@.Data[] <- .Call(R_mpfr_igamma, a, x, match.arg(rnd.mode))
	    a
	} else { ## x is not (mpfr)Array
	    if(length(a) %% length(x) != 0)
		stop("length of first argument (array) is not multiple of the second argument's one")
	    a@.Data[] <- .Call(R_mpfr_igamma, a, as(x, "mpfr"), match.arg(rnd.mode))
	    a
	}
    } else if(is(x, "mpfrArray")) {
	if(length(x) %% length(a) != 0)
	    stop("length of second argument (array) is not multiple of the first argument's one")
	x@.Data[] <- .Call(R_mpfr_igamma, as(a, "mpfr"), x, match.arg(rnd.mode))
	x
    }
    else
	new("mpfr", .Call(R_mpfr_igamma, as(a, "mpfr"), as(x, "mpfr"), match.arg(rnd.mode)))
}

## only as long as we still may have   mpfrVersion() < "3.2.0", e.g. in Fedora 30 (2019f)

## mpfrVersion() cannot be called at package build time  (underlying C entry point not ready):
## if(mpfrVersion() < "3.2.0")
## dummy .. to pacify "R CMD check"
## R_mpfr_igamma <- quote(dummy) # gives NOTE  ‘R_mpfr_igamma’ is of class "name"


## These are identical from package copuula/R/special-func.R -- where MM authored the function also:
## We want to export these, but cannot easily import from copula which "weekly depends" on Rmpfr

##' @title Compute  f(a) = log(1 - exp(-a))  stably
##' @param a numeric vector of positive values
##' @param cutoff  log(2) is optimal, see  Maechler (201x) .....
##' @return f(a) == log(1 - exp(-a)) == log1p(-exp(-a)) == log(-expm1(-a))
##' @author Martin Maechler, May 2002 .. Aug. 2011
##' @references Maechler(2012)
##' Accurately Computing log(1 - exp(-|a|)) Assessed by the Rmpfr package.
##' http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
## MM: ~/R/Pkgs/Rmpfr/inst/doc/log1mexp-note.Rnw
##--> ../man/log1mexp.Rd
log1mexp <- function(a, cutoff = log(2)) ## << log(2) is optimal >>
{
    if(has.na <- any(ina <- is.na(a))) {
	y <- a
	a <- a[ok <- !ina]
    }
    if(any(a < 0))## a == 0  -->  -Inf	(in both cases)
	warning("'a' >= 0 needed")
    tst <- a <= cutoff
    r <- a
    r[ tst] <- log(-expm1(-a[ tst]))
    r[!tst] <- log1p(-exp(-a[!tst]))
    if(has.na) { y[ok] <- r ; y } else r
}

##' @title Compute  f(x) = log(1 + exp(x))  stably and quickly
##--> ../man/log1mexp.Rd
log1pexp <- function(x, c0 = -37, c1 = 18, c2 = 33.3)
{
    if(has.na <- any(ina <- is.na(x))) {
	y <- x
	x <- x[ok <- !ina]
    }
    r <- exp(x)
    if(any(i <- c0 < x & (i1 <- x <= c1)))
	r[i] <- log1p(r[i])
    if(any(i <- !i1 & (i2 <- x <= c2)))
	r[i] <- x[i] + 1/r[i] # 1/exp(x) = exp(-x)
    if(any(i3 <- !i2))
	r[i3] <- x[i3]
    if(has.na) { y[ok] <- r ; y } else r
}
