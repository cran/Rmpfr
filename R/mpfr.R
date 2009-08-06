#### All methods for  "mpfr" (and "mpfr1") class
#### apart from coercions and the group methods

setMethod("is.finite", "mpfr",
          function(x) .Call("R_mpfr_is_finite", x, PACKAGE="Rmpfr"))
setMethod("is.infinite", "mpfr",
          function(x) .Call("R_mpfr_is_infinite", x, PACKAGE="Rmpfr"))
## MPFR has only "NaN" ( == "NA"  -- hence these two are identical :
setMethod("is.na", "mpfr",
          function(x) .Call("R_mpfr_is_na", x, PACKAGE="Rmpfr"))
setMethod("is.nan", "mpfr",
          function(x) .Call("R_mpfr_is_na", x, PACKAGE="Rmpfr"))

mpfr.is.0 <- function(x) .Call("R_mpfr_is_zero", x, PACKAGE="Rmpfr")
    ## sapply(x, function(.) .@exp == - .Machine$integer.max)

mpfr_default_prec <- function(prec) {
    if(missing(prec) || is.null(prec))
	.Call("R_mpfr_get_default_prec", PACKAGE="Rmpfr")
    else {
	stopifnot((prec <- as.integer(prec[1])) > 0)
	.Call("R_mpfr_set_default_prec", prec, PACKAGE="Rmpfr")
    }
}

print.mpfr1 <- function(x, digits = NULL, ...) {
    stopifnot(is(x, "mpfr1"), is.null(digits) || digits >= 2)
    cat("'mpfr1' ", format(as(x, "mpfr"), digits=digits),"\n", sep="")
    invisible(x)
}

setMethod(show, "mpfr1", function(object) print.mpfr1(object))

## For testing, debugging etc
if(.Platform$OS.type != "windows") {## No R_Outputfile (in C) on Windows

.print.mpfr <- function(x, digits = NA, ...) {
    stopifnot(is(x, "mpfr"), is.na(digits) || digits >= 2)
    ## digits = NA --> the inherent precision of x will be used
    if(length(x) >= 1)
	.Call("print_mpfr", x, as.integer(digits), PACKAGE="Rmpfr")
    invisible(x)
}
}# non-Windows only

## Get or Set the C-global  'R_mpfr_debug_' variable:
.mpfr.debug <- function(i = NA)
    .Call("R_mpfr_set_debug", as.integer(i), PACKAGE="Rmpfr")

print.mpfr <- function(x, digits = NULL, ...) {
    stopifnot(is(x, "mpfr"), is.null(digits) || digits >= 2)
    ## digits = NA --> the inherent precision of x will be used
    n <- length(x)
    ch.prec <-
	if(n >= 1) {
	    rpr <- range(sapply(x, slot, "prec"))
	    paste("of precision ", rpr[1],
		   if(rpr[1] != rpr[2]) paste("..",rpr[2]), " bits")
	}
    cat(n, "'mpfr'", if(n == 1) "number" else "numbers", ch.prec, "\n")
    if(n >= 1)
	print(format(x, digits=digits), ..., quote = FALSE)
    ## .Call("print_mpfr", x, as.integer(digits), PACKAGE="Rmpfr")
    invisible(x)
}
setMethod(show, "mpfr", function(object) print.mpfr(object))


## "[" which also keeps names ... JMC says that names are not support(ed|able)
## ---  for such objects..
setMethod("[", signature(x = "mpfr", i = "ANY", j = "missing", drop = "missing"),
          function(x,i,j, ..., drop) {
              nA <- nargs()
	      if(nA == 2) { ## x[i] etc -- vector case -- to be fast, need C! --
		  x@.Data <- structure(x@.Data[i], names=names(x)[i])
		  x
	      } else if(nA == 3 && !is.null(d <- dim(x))) { ## matrix indexing(!)
                  ## not keeping dimnames though ...
                  message("nargs() == 3  'mpfr' array indexing ... ")
                  new("mpfr", structure(x@.Data[i,j,...,drop=drop], dim = d))
## keeping dimnames: maybe try
##                   D <- x@.Data; dim(D) <- d
##                   if(!is.null(dn <- dimnames(x))) dimnames(D) <- dn
##                   D <- D[i,,drop=drop]
##                   new("mpfr", D)

              }
              else
                  stop(sprintf("invalid 'mpfr' subsetting (nargs = %d)",nA))
          })

## "[<-" :
setReplaceMethod("[", signature(x = "mpfr", i = "ANY", j = "missing",
				value = "ANY"),
	  function(x,i,value) { x[i] <- as(value, "mpfr"); x })

setReplaceMethod("[", signature(x = "mpfr", i = "ANY", j = "missing",
				value = "mpfr"),
	  function(x,i, ..., value) {
	      if(length(list(...)))# should no longer happen:
		  stop("extra replacement arguments",
		       deparse(list(...))," not dealt with")
	      x@.Data[i] <- value
	      x })



## I don't see how I could use setMethod("c", ...)
## but this works "magically"  when the first argument is an mpfr :
c.mpfr <- function(...) new("mpfr", unlist(lapply(list(...), as, Class = "mpfr")))


setMethod("unique", signature(x="mpfr", incomparables="missing"),
	  function(x, incomparables = FALSE, ...)
	  new("mpfr", unique(x@.Data, incomparables, ...)))

## -> duplicated() now work

## sort() works too  (but could be made faster via faster
## ------  xtfrm() method !  [ TODO ]


setGeneric("pmin", signature = "...")# -> message about override ...
setGeneric("pmax", signature = "...")

setMethod("pmin", "Mnumber",
	  function(..., na.rm = FALSE) {
	      args <- list(...)
	      if(all(sapply(args, is.atomic)))
		  return( base::pmin(..., na.rm = na.rm) )
	      ## else: at least one is "mpfr(Matrix/Array)"
              is.m <- sapply(args, is, "mpfr")
              if(!any(is.m))
                  stop("no \"mpfr\" argument -- wrong method chosen")

	      N <- max(lengths <- sapply(args, length))
              ## precision needed -- FIXME: should be *vector*
              mPrec <- max(unlist(lapply(args[is.m], getPrec)),
                           if(any(sapply(args[!is.m], is.double)))
                           .Machine$double.digits)
              ## to be the result :
              r <- mpfr(rep.int(Inf, N), precBits = mPrec)

              ## modified from ~/R/D/r-devel/R/src/library/base/R/pmax.R
              has.na <- FALSE
              for(i in seq_along(args)) {
                  x <- args[[i]]
                  if((n.i <- lengths[i]) != N)
                      x <- x[rep(seq_len(n.i), length.out = N)]
                  nas <- cbind(is.na(r), is.na(x))
                  if(!is.m[i]) x <- mpfr(x, prec=mPrec)
                  if(has.na || (has.na <- any(nas))) {
                      r[nas[, 1L]] <- x[nas[, 1L]]
                      x[nas[, 2L]] <- r[nas[, 2L]]
                  }
                  change <- r > x
                  change <- change & !is.na(change)
                  r[change] <- x[change]
                  if (has.na && !na.rm)
                      r[nas[, 1L] | nas[, 2L]] <- NA
              }

              mostattributes(r) <- attributes(args[[1L]])
              r
          })

setMethod("pmax", "Mnumber",
	  function(..., na.rm = FALSE) {
	      args <- list(...)
	      if(all(sapply(args, is.atomic)))
		  return( base::pmax(..., na.rm = na.rm) )
	      ## else: at least one is "mpfr(Matrix/Array)"
              is.m <- sapply(args, is, "mpfr")
              if(!any(is.m))
                  stop("no \"mpfr\" argument -- wrong method chosen")

	      N <- max(lengths <- sapply(args, length))
              ## precision needed -- FIXME: should be *vector*
              mPrec <- max(unlist(lapply(args[is.m], getPrec)),
                           if(any(sapply(args[!is.m], is.double)))
                           .Machine$double.digits)
              ## to be the result :
              r <- mpfr(rep.int(-Inf, N), precBits = mPrec)

              ## modified from ~/R/D/r-devel/R/src/library/base/R/pmax.R
              has.na <- FALSE
              for(i in seq_along(args)) {
                  x <- args[[i]]
                  if((n.i <- lengths[i]) != N)
                      x <- x[rep(seq_len(n.i), length.out = N)]
                  nas <- cbind(is.na(r), is.na(x))
                  if(!is.m[i]) x <- mpfr(x, prec=mPrec)
                  if(has.na || (has.na <- any(nas))) {
                      r[nas[, 1L]] <- x[nas[, 1L]]
                      x[nas[, 2L]] <- r[nas[, 2L]]
                  }
                  change <- r < x
                  change <- change & !is.na(change)
                  r[change] <- x[change]
                  if (has.na && !na.rm)
                      r[nas[, 1L] | nas[, 2L]] <- NA
              }

              mostattributes(r) <- attributes(args[[1L]])
              r
          })


### seq() :

## seq.default()  and  seq.Date()  as examples :
## ~/R/D/r-devel/R/src/library/base/R/seq.R    and
## ~/R/D/r-devel/R/src/library/base/R/dates.R

seqMpfr <- function(from = 1, to = 1, by = ((to - from)/(length.out - 1)),
                    length.out = NULL, along.with = NULL, ...)
{

    if(missing(from)) stop("'from' must be specified")
    if (!is(from, "mpfr")) from <- as(from, "mpfr")
    if(length(from) != 1) stop("'from' must be of length 1")
    if(!missing(to)) {
        if (!is(to, "mpfr")) to <- as(to, "mpfr")
        if (length(to) != 1) stop("'to' must be of length 1")
    }
    if (!missing(along.with)) {
        length.out <- length(along.with)
    } else if (!is.null(length.out)) {
        if (length(length.out) != 1) stop("'length.out' must be of length 1")
        length.out <- ceiling(length.out)
    }
##     status <- c(!missing(to), !missing(by), !is.null(length.out))
##     if(sum(status) != 2)
## ## stop("exactly two of 'to', 'by' and 'length.out' / 'along.with' must be specified")
##         warning("not exactly two of 'to', 'by' and 'length.out' / 'along.with' have been specified")

    if(is.null(length.out) && missing(by))
        by <- mpfr(1, from[[1]]@prec)

    if (!is(by, "mpfr")) by <- as(by, "mpfr")
    if (length(by) != 1) stop("'by' must be of length 1")

    ## ---- This is  cut n paste  from seq.default() :
    ## ---- It should work, since "arithmetic works for mpfr :
    if(is.null(length.out)) {
        del <- to - from
        if(del == 0 && to == 0) return(to)
        n <- del/by
        if(!(length(n) && is.finite(n))) {
            if(length(by) && by == 0 && length(del) && del == 0)
                return(from)
            stop("invalid (to - from)/by in seq(.)")
        }
        if(n < 0)
            stop("wrong sign in 'by' argument")
        if(n > .Machine$integer.max)
            stop("'by' argument is much too small")

        dd <- abs(del)/max(abs(to), abs(from))
        if (dd < 100*.Machine$double.eps) return(from)
        n <- as.integer(n + 1e-7)
        x <- from + (0:n) * by
        ## correct for overshot because of fuzz
        if(by > 0) pmin(x, to) else pmax(x, to)
    }
    else if(!is.finite(length.out) || length.out < 0)
	stop("length must be non-negative number")
    else if(length.out == 0)
	from[FALSE] # of same precision
    ## else if (One) 1:length.out
    else if(missing(by)) {
	# if(from == to || length.out < 2) by <- 1
	if(missing(to))
	    to <- from + length.out - 1
	if(missing(from))
	    from <- to - length.out + 1
	if(length.out > 2)
	    if(from == to)
		rep.int(from, length.out)
	    else as.vector(c(from, from + (1:(length.out - 2)) * by, to))
	else as.vector(c(from, to))[1:length.out]
    }
    else if(missing(to))
	from + (0:(length.out - 1)) * by
    else if(missing(from))
	to - ((length.out - 1):0) * by
    else stop("too many arguments")
}

if(FALSE) ## fails: seq(1, length.out=3)
setGeneric("seq", function(from, to, by, ...) standardGeneric("seq"),
           useAsDefault = function(from, to, by, ...)
           base::seq(from, to, by, ...))
if(FALSE) ## fails: seq(1, length.out=3)
setGeneric("seq", function(from, to, by, ...) standardGeneric("seq"),
           useAsDefault =
           function(from=1, to=1, by=((to-from)/(length.out-1)), ...)
           base::seq(from, to, by, ...))

if(FALSE) { ##-- but this also fails: afterwards  seq(1, length.out=3)

setGeneric("seq", function (from, to, by, length.out, along.with, ...)
           standardGeneric("seq"),
           signature = c("from", "to", "by"),
           useAsDefault = {
               function(from=1, to=1, by=((to-from)/(length.out-1)),
                        length.out=NULL, along.with=NULL, ...)
                   base::seq(from, to, by,
                             length.out=length.out, along.with=along.with, ...)
           })


setMethod("seq", c(from="mpfr", to="ANY", by = "ANY"), seqMpfr)
setMethod("seq", c(from="ANY", to="mpfr", by = "ANY"), seqMpfr)
setMethod("seq", c(from="ANY", to="ANY", by = "mpfr"), seqMpfr)

}#not yet

getPrec <- function(x) sapply(x, slot, "prec")

### all.equal()

## TODO ?? <<<<<<<<<<<
## ====
## 2) instead of  as(., "mpfr")  use  mpfr(., precBits = <smart>)

## For two "mpfr"s, use a  "smart" default tolerance :
setMethod("all.equal", signature(target = "mpfr", current = "mpfr"),
	  function (target, current,
		    tolerance =
		    2^-(0.5 * min(mean(getPrec(target)),
				  mean(getPrec(current)))), ...)
      {
	  ## to use "our" mean() :
	  environment(all.equal.numeric) <- environment()
	  all.equal.numeric(target, current, tolerance=tolerance, ...)
      })

setMethod("all.equal", signature(target = "mpfr", current = "ANY"),
	  function (target, current,
		    tolerance = .Machine$double.eps^0.5, ...) {
	      ## to use "our" mean() :
	      environment(all.equal.numeric) <- environment()
	      all.equal.numeric(target, as(current, "mpfr"),
				tolerance=tolerance, ...)
	  })

setMethod("all.equal", signature(target = "ANY", current = "mpfr"),
	  function (target, current,
		    tolerance = .Machine$double.eps^0.5, ...) {
	      ## to use "our" mean() :
	      environment(all.equal.numeric) <- environment()
	      all.equal.numeric(as(target, "mpfr"), current,
				tolerance=tolerance, ...)
	  })
