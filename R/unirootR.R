### This is a translation of R_zeroin2 in ~/R/D/r-devel/R/src/appl/zeroin.c
### from  C to R  by John Nash,
### ---> file rootoned/R/zeroin.R of the new (2011-08-18) R-forge package rootoned
###
### Where John Nash calls it zeroin(), I call it unirootR()

##' Simple modification of uniroot()  which should work with mpfr-numbers
##' MM: uniroot() is in ~/R/D/r-devel/R/src/library/stats/R/nlm.R
##'
unirootR <- function(f, interval, ...,
		     lower = min(interval), upper = max(interval),
		     f.lower = f(lower, ...), f.upper = f(upper, ...),
                     extendInt = c("no", "yes", "downX", "upX"),
                     trace = 0,
		     verbose = as.logical(trace), verbDigits = max(3, min(20, -log10(tol)/2)),
		     tol = .Machine$double.eps^0.25, maxiter = 1000L,
                     check.conv = FALSE,
                     ## Rmpfr-only:
                     warn.no.convergence = !check.conv,
		     epsC = NULL)
{
    if(!missing(interval) && length(interval) != 2L)
	stop("'interval' must be a vector of length 2")
    ## For many "quick things", we will use as.numeric(.) but we do *NOT* assume that
    ## lower and upper are numeric!
    .N <- as.numeric
    if(lower >= upper) # (may be mpfr-numbers *outside* double.xmax)
	stop("lower < upper  is not fulfilled")
    if(is.na(.N(f.lower))) stop("f.lower = f(lower) is NA")
    if(is.na(.N(f.upper))) stop("f.upper = f(upper) is NA")

    form  <- function(x, digits = verbDigits) format(x, digits=digits, drop0trailing=TRUE)
    formI <- function(x, di = getOption("digits")) format(x, digits=di, drop0trailing=TRUE)
    Sig <- switch(match.arg(extendInt),
		  "yes" = NULL,
		  "downX"= -1,
		  "no"   =  0,
		  "upX"  =  1,
		  stop("invalid 'extendInt'; please report"))
    ## protect against later   0 * Inf  |--> NaN  and Inf * -Inf.
    truncate <- function(x) { ## NA are already excluded; deal with +/- Inf
        if(is.numeric(x))
            pmax.int(pmin(x, .Machine$double.xmax), -.Machine$double.xmax)
        else if(inherits(x, "mpfr") && is.infinite(x)) # use maximal/minimal mpfr-number instead:
            sign(x) * mpfr(2, .getPrec(x))^((1 - 2^-52)*.mpfr_erange("Emax"))
        else x
    }
    f.low. <- truncate(f.lower)
    f.upp. <- truncate(f.upper)
    doX <- (   is.null(Sig) && f.low. * f.upp. > 0 ||
	    is.numeric(Sig) && (Sig*f.low. > 0 || Sig*f.upp. < 0))
    if(doX) { ## extend the interval = [lower, upper]
	if(trace)
	    cat(sprintf("search {extendInt=\"%s\", Sig=%s} in [%s,%s]%s",
                        extendInt, formI(Sig), formI(lower), formI(upper),
			if(trace >= 2)"\n" else " ... "))
	Delta <- function(u) 0.01* pmax(1e-4, abs(u)) ## <-- FIXME? [= R's uniroot() for double]
        it <- 0L
	## Two cases:
	if(is.null(Sig)) {
	    ## case 1)	'Sig' unspecified --> extend (lower, upper) at the same time
	    delta <- Delta(c(lower,upper))
	    while(isTRUE(f.lower*f.upper > 0) &&
                  any(iF <- is.finite(c(lower,upper)))) {
		if((it <- it + 1L) > maxiter)
		    stop(gettextf("no sign change found in %d iterations", it-1),
			 domain=NA)
		if(iF[1]) {
		    ol <- lower; of <- f.lower
		    if(is.na(f.lower <- f(lower <- lower - delta[1], ...))) {
			lower <- ol; f.lower <- of; delta[1] <- delta[1]/4
		    }
		}
		if(iF[2]) {
		    ol <- upper; of <- f.upper
		    if(is.na(f.upper <- f(upper <- upper + delta[2], ...))) {
			upper <- ol; f.upper <- of; delta[2] <- delta[2]/4
		    }
		}
		if(trace >= 2)
		    cat(sprintf(" .. modified lower,upper: (%15g,%15g)\n",
				.N(lower), .N(upper)))
		delta <- 2 * delta
	    }
	} else {
	    ## case 2) 'Sig' specified --> typically change only *one* of lower, upper
	    ## make sure we have Sig*f(lower) <= 0 and Sig*f(upper) >= 0:
	    delta <- Delta(lower)
	    while(isTRUE(Sig*f.lower > 0)) {
		if((it <- it + 1L) > maxiter)
		    stop(gettextf("no sign change found in %d iterations", it-1),
			 domain=NA)
		f.lower <- f(lower <- lower - delta, ...)
		if(trace >= 2) cat(sprintf(" .. modified lower: %s, f(.)=%s\n",
                                           formI(lower), formI(f.lower)))
		delta <- 2 * delta
	    }
	    delta <- Delta(upper)
	    while(isTRUE(Sig*f.upper < 0)) {
		if((it <- it + 1L) > maxiter)
		    stop(gettextf("no sign change found in %d iterations", it-1),
			 domain=NA)
		f.upper <- f(upper <- upper + delta, ...)
		if(trace >= 2) cat(sprintf(" .. modified upper: %s, f(.)=%s\n",
                                           formI(upper), formI(f.upper)))
		delta <- 2 * delta
	    }
	}
	if(trace && trace < 2)
            cat(sprintf("extended to [%s, %s] in %d steps\n", formI(lower), formI(upper), it))
    }
    if(!isTRUE(as.vector(sign(f.lower) * sign(f.upper) <= 0)))
	stop(if(doX)
	"did not succeed extending the interval endpoints for f(lower) * f(upper) <= 0"
	     else sprintf("f() values at end points = (%s, %s) not of opposite sign",
                          formI(f.lower), formI(f.upper)))

    if(is.null(epsC) || is.na(epsC)) {
	## determine 'epsC'  ``the achievable Machine precision''
	## -- given the class of f.lower, f.upper
        ff <- f.lower * f.upper
	if(is.double(ff)) epsC <- .Machine$double.eps
	else if(is(ff, "mpfr"))
	    epsC <- 2^-min(getPrec(f.lower), getPrec(f.upper))
	else { ## another number class -- try to see if getPrec() is defined..
	    ## if not, there's not much we can do
	    if(is(prec <- tryCatch(min(getPrec(f.lower), getPrec(f.upper)),
				   error = function(e)e),
		  "error")) {
		warning("no valid getPrec() for the number class(es) ",
			paste(unique(class(f.lower),class(f.upper)), collapse=", "),
			".\n Using double precision .Machine$double.eps.")
		epsC <- .Machine$double.eps
	    } else {
		epsC <- 2^-prec
		message("using epsC = %s ..", format(epsC))
	    }
	}
    }
    if(tol < epsC / 8) # "8 fudge factor" (otherwise happens too often)
        warning(sprintf("tol (%g) < epsC (%g)  is rarely sensical,
 and the resulting precision is probably not better than epsC",
                        tol, epsC))

    ## Instead of the call to C code, now "do it in R" :
    ## val <- .Internal(zeroin2(function(arg) as.numeric(f(arg, ...)),
    ##			     lower, upper, f.lower, f.upper,
    ##			     tol, as.integer(maxiter)))
    a <- lower # interval[1]
    b <- upper # interval[2]
    fa <- f.lower # f(ax, ...)
    fb <- f.upper # f(bx, ...)
    if (verbose)
	cat(sprintf("==> Start zeroin: f(%g)= %g;  f(%g)= %g\n", .N(a), .N(fa),  .N(b), .N(fb)))
    c <- a
    fc <- fa
    ## First test if we have found a root at an endpoint
    maxit <- maxiter + 2L # count evaluations as maxiter-maxit
    converged <- FALSE

    while(!converged && maxit > 0) { ##---- Main iteration loop ------------------------------
	if (verbose) cat("Iteration >>>", maxiter+3L-maxit, "<<< ;")
	d.prev <- b-a
	## Distance from the last but one to the last approximation	*/
	##double tol.2;		       ## Actual tolerance		*/
	##double p;			## Interpolation step is calcu- */
	##double q;			## lated in the form p/q; divi-
	##				 * sion operations is delayed
	##				 * until the last moment	*/
	##double d.new;		## Step at this iteration	*/

	if(abs(fc) < abs(fb)) { ## Swap data for b to be the smaller
	    if (verbose) cat(sprintf("fc (=%s) smaller than fb\n", form(fa)))
	    a <- b
	    b <- c
	    c <- a	## best approximation
	    fa <- fb
	    fb <- fc
	    fc <- fa
	}
	tol.2 <- 2*epsC*abs(b) + tol/2
	d.new <- (c-b)/2 # bisection
	if (verbose) cat("tol.2(epsC,b) = ",.N(tol.2), "; d.new= ",.N(d.new),"\n", sep="")

	## converged <- (abs(d.new) <= tol.2 && is.finite(fb))  ||  fb == 0
	converged <- (abs(d.new) <= tol.2)  ||  fb == 0
	if(converged) {
	    if (verbose) cat("DONE! -- small d.new or fb=0\n")
	    ## Acceptable approx. is found :
	    val <- list(root=b, froot=fb, rtol = abs(c-b), maxit=maxiter-maxit)
	}
	else {
	    ## Decide if the interpolation can be tried	*/
	    if( (abs(d.prev) >= tol.2) ## If d.prev was large enough*/
	       && (abs(fa) > abs(fb)) ) { ## and was in true direction,
		## Interpolation may be tried	*/
		##    register double t1,cb,t2;
		if (verbose) cat("d.prev larger than tol.2 and fa bigger than fb --> ")

		cb <- c-b
		if (a == c) { ## If we have only two distinct points, linear interpolation
		    ## can only be applied
		    t1 <- fb/fa
		    p <- cb*t1
		    q <- 1 - t1
		    if (verbose) cat("a == c: ")
		}
		else { ## Quadric inverse interpolation*/
		    if (verbose) cat("a != c: ")
		    q <- fa/fc
		    t1 <- fb/fc
		    t2 <- fb/fa
		    p <- t2 * ( cb*q*(q-t1) - (b-a)*(t1-1) )
		    q <- (q-1) * (t1-1) * (t2-1)
		}
		if(p > 0) { ## p was calculated with the */
		    if (verbose) cat(" p > 0; ")
		    q <- -q ## opposite sign; make p positive */
		} else {    ## and assign possible minus to	*/
		    if (verbose) cat(" p <= 0; ")
		    p <- -p ## q				*/
		}
		if (p < 0.75*cb*q - abs(tol.2*q)/2 ## If b+p/q falls in [b,c]*/
		    && p < abs(d.prev*q/2)) {	## and isn't too large	*/
		    if (verbose) cat("p satisfies conditions for changing d.new\n")
		    d.new <- p/q ## it is accepted
		} else if(verbose) cat("\n")
		## If p/q is too large, then the
		## bisection procedure can reduce [b,c] range to more extent
	    }

	    if( abs(d.new) < tol.2) { ## Adjust the step to be not less than tolerance
		if (verbose) cat("d.new smaller than tol.2, adjusted to it.\n")
                d.new <- if(d.new > 0) tol.2 else -tol.2
	    }
	    a <- b
	    fa <- fb ## Save the previous approx. */
	    b <- b + d.new
	    fb <- f(b, ...)
            if (verbose) cat(sprintf("new f(b=%s) = %s;\n", form(b), form(fb)))
	    maxit <- maxit-1
	    ## Do step to a new approxim. */
	    if( ((fb > 0) && (fc > 0)) || ((fb < 0) && (fc < 0)) ) {
		if (verbose) cat(sprintf("  make c:=a=%s to have sign opposite to b: f(c)=%s\n",
                                         form(a), form(fa)))
		## Adjust c for it to have a sign opposite to that of b */
		c <- a
		fc <- fa
            }

	}## else not converged

    } ## end{ while(maxit > 0) } --------------------------------------------

    if(converged) {
	iter <- val[["maxit"]]
	if(!is.na(fb) &&  abs(fb) > 0.5*max(abs(f.lower), abs(f.upper)))# from John Nash:
	    warning("Final function magnitude seems large -- maybe converged to sign-changing 'pole' location?")
    } else { ## (!converged) : failed!
        if(check.conv)
            stop("no convergence in zero finding in ", iter, " iterations")
        ## else
	val <- list(root= b, rtol = abs(c-b))
	iter <- maxiter
	if(warn.no.convergence)
	    warning("_NOT_ converged in ", iter, " iterations")
    }
    list(root = val[["root"]], f.root = f(val[["root"]], ...),
	 iter = iter, estim.prec = .N(val[["rtol"]]), converged = converged)
} ## {unirootR}


### qnorm() via inversion of pnorm()  by  unirootR() ==========================
"
qnorm(p) == q   <==>  pnorm(q) == p
"

##====== TODO: Use good i.e. *tight* inequalities for  Phi(x) ..  which are invertible
##       ===== to (still tight) inequalities for Phi^{-1}(x) == qnorm(x)
## ===> get good *start* interval for unirootR()

## qnorm  (p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
qnormI <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE,
                   trace = 0, verbose = as.logical(trace), # <- identical defaults as unirootR()
                   tol, # for "base" = .Machine$double.eps^0.25, but here use getPrec()
                   useMpfr = any(prec > 53), # if true, use mpfr
                   give.full = FALSE,
                   ...) # <- arguments to pass to unirootR()
{
    ## The function  whose "zeros" aka "roots" we want to find:
    zFun <- function(q) pnorm(q, mean=mean, sd=sd, lower.tail=lower.tail, log.p=log.p) - p.
    if(missing(tol) || !is.finite(tol)) {
        prec <- max(getPrec(if(missing(mean) && missing(sd)) p else p+(mean+sd)))
        ## not max(getPrec(c(p, mean, sd)))  as it's >= 128 by default
        tol <- 2^-(prec + 2) # 2^(-prec - 1) gives less accurate
    } else
        prec <- as.integer(ceiling(1 - log2(tol)))
    ##
    if(verbose) cat(sprintf("prec=%d ==> useMpfr=%s:\n", prec, format(useMpfr)))
    verbDigs <- if(any(is.vd <- "verbDigits" == ...names()))
                 ...elt(which(is.vd))
             else max(3, min(20, -log10(tol)/2))
    if(useMpfr) { ## This **IS** important here:
        old_eranges <- .mpfr_erange() # typically -/+ 2^30
        myERng <- (1-2^-52) * .mpfr_erange(c("min.emin","max.emax"))
        if(!isTRUE(all.equal(myERng, old_eranges))) {
            .mpfr_erange_set(value = myERng)
            on.exit( .mpfr_erange_set( , old_eranges) )
        }
    }
    sgn <- if(lower.tail) -1 else 1
    INf <- if(lower.tail) Inf else -Inf
    ## "start"-interval for unirootR() --- see 'TODO' .. *tight* .. above <<<<<<<<
    ## Notably for the relevant (p <<< -1, log.p=TRUE) case !
    qnInt <-
        if(log.p) {
            Pi <- if(useMpfr) Const("pi", prec) else pi
            function(p) { s2 <- -2*p
                if(useMpfr && !inherits(s2, "mpfr")) s2 <- mpfr(s2, precBits = prec)
                xs1 <- s2 - log(2*Pi*s2)
                if(p < -54) {
                    qn <- sqrt(s2 - log(2*Pi*xs1)); e <- 1e-4
                } else if(p <= -15) {
                    qn <- sqrt(s2 - log(2*Pi*xs1) - 1/(2 + xs1)); e <- 1e-3
                } else { # p >= -15
                    qn <- stats__qnorm(asNumeric(p), log.p=TRUE)
                    ## FIXME: not good enough, e.g., for  p = mpfr(-1e-5, 128)
                    e <- 1e-2
                }
                qn*(sgn + c(-e,e))
            }
        } else { ## log.p is FALSE
            Id <- if(useMpfr) function(.) mpfr(., precBits = prec) else identity
            function(p) { q <- stats__qnorm(asNumeric(p), lower.tail=lower.tail)
                Id(if(abs(q) < 1e-3) c(-2,2)*1e-3 else c(.99, 1.01) * q)
            }
        }
    ## Deal with prob in {0, 1} which correspond to quantiles -/+ Inf :
    ## idea from {DPQ}'s .D_0 and .D_1 :
    .p_1 <- as.integer(!log.p)
    .p_0 <- if (log.p) -Inf else 0
    r <- if(give.full) vector("list", length(p))
         else if(useMpfr) mpfr(p, precBits=prec)
         else p
    for(ip in seq_along(p)) {
        p. <- p[ip]
        if(verbose) cat(sprintf("p. = p[ip=%d] = %s:\n", ip,
                                format(p., digits = verbDigs, drop0trailing=TRUE)))
        ri <-
            if (p. == .p_1)
                INf
            else if (p. == .p_0)
                -INf
            else if(is.na(p.) || p. > .p_1 || p. < .p_0)
                NaN
            else { ## zFun() uses p.
                ur <- unirootR(zFun, interval = qnInt(p.),
                               extendInt = if(lower.tail) "upX" else "downX",
                               trace=trace, verbose=verbose, tol=tol,
                               ...) # verbDigits, maxiter, check.conv, warn.no.convergence, epsC
                if(give.full) ur else ur$root
            }
        if(give.full) r[[ip]] <- ri else r[ip] <- ri
    }
    r
}
