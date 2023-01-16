#### Check tools -- notably for Rmpfr
#### ================================
## to be used as
##	source(system.file("check-tools.R", package="Rmpfr"), keep.source=FALSE)

## Some of the Matrix package check tools {MM = ~/R/Pkgs/Matrix/inst/test-tools-1.R}
##                                                       ~~~~~~~~~~~~~~~~~~~~~~~~~~
assert.EQ <- function(target, current, tol = if(showOnly) 0 else 1e-15,
                      giveRE = FALSE, showOnly = FALSE, ...) {
    ## Purpose: check equality *and* show non-equality
    ## ----------------------------------------------------------------------
    ## showOnly: if TRUE, return (and hence typically print) all.equal(...)
    T <- isTRUE(ae <- all.equal(target, current, tolerance = tol, ...))
    if(showOnly) return(ae) else if(giveRE && T) { ## don't show if stop() later:
	ae0 <- if(tol == 0) ae else all.equal(target, current, tolerance = 0, ...)
	if(!isTRUE(ae0)) writeLines(ae0)
    }
    if(!T) stop("all.equal() |-> ", paste(ae, collapse=sprintf("%-19s","\n")))
    else if(giveRE) invisible(ae0)
}
##' a version with other "useful" defaults (tol, giveRE, check.attr..)
assert.EQ. <- function(target, current,
		       tol = if(showOnly) 0 else .Machine$double.eps^0.5,
		       giveRE = TRUE, showOnly = FALSE, ...) {
    assert.EQ(target, current, tol=tol, giveRE=giveRE, showOnly=showOnly,
	      check.attributes=FALSE, ...)
}

showSys.time <- function(expr, ...) {
    ## prepend 'Time' for R CMD Rdiff
    st <- system.time(expr, ...)
    writeLines(paste("Time", capture.output(print(st))))
    invisible(st)
}



### ------- Part I --  do not need 'Rmpfr'

`%=N=%` <- function(x,y) (x == y) | (is.na(x) & is.na(y))

all.eq.finite <- function(x,y, ...) {
    ## x = 'target'   y = 'current'
    if(any(is.finite(y[!(fx <- is.finite(x))])))
	return("current has finite values where target has not")
    if(any(is.finite(x[!(fy <- is.finite(y))])))
	return("target has finite values where current has not")
    ## now they have finite values at the same locations
    all.equal(x[fx], y[fy], ...)
}

### ------- Part II --  do not make sense or work outside of 'Rmpfr' :

all.EQ <- function(x,y, tolerance = 2^-98, ...) # very small tol. for MPFR
    all.equal.finite(x, y, tolerance=tolerance, ...)

