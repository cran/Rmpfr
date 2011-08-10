#### All methods for  "mpfr" (and "mpfr1") class
#### apart from coercions and the group methods

setMethod("is.finite", "mpfr",
          function(x) .Call(R_mpfr_is_finite, x))
setMethod("is.infinite", "mpfr",
          function(x) .Call(R_mpfr_is_infinite, x))
## MPFR has only "NaN" ( == "NA"  -- hence these two are identical :
setMethod("is.na", "mpfr",
          function(x) .Call(R_mpfr_is_na, x))
setMethod("is.nan", "mpfr",
          function(x) .Call(R_mpfr_is_na, x))

mpfr.is.0 <- function(x) .Call(R_mpfr_is_zero, x)
    ## sapply(x, function(.) .@exp == - .Machine$integer.max)

mpfr.is.integer <- function(x)
    .Call(R_mpfr_is_integer, x)

is.whole <- function(x) {
    if(is.integer(x) || is.logical(x)) rep.int(TRUE, length(x))
    else if(is.numeric(x)) x == floor(x)
    else if(is.complex(x)) x == round(x)
    else if(is(x,"mpfr")) mpfr.is.integer(x)
    else rep.int(FALSE, length(x))
}

mpfr_default_prec <- function(prec) {
    if(missing(prec) || is.null(prec))
	.Call(R_mpfr_get_default_prec)
    else {
	stopifnot((prec <- as.integer(prec[1])) > 0)
	.Call(R_mpfr_set_default_prec, prec)
    }
}

.mpfrVersion <- function() .Call(R_mpfr_get_version)
mpfrVersion <- function()
    numeric_version(sub("^([0-9]+\\.[0-9]+\\.[0-9]+).*","\\1", .mpfrVersion()))

print.mpfr1 <- function(x, digits = NULL, drop0trailing = TRUE, ...) {
    stopifnot(is(x, "mpfr1"), is.null(digits) || digits >= 2)
    cat("'mpfr1' ",
	format(as(x, "mpfr"), digits=digits, drop0trailing=drop0trailing),
	"\n", sep="")
    invisible(x)
}

setMethod(show, "mpfr1", function(object) print.mpfr1(object))

## For testing, debugging etc
if(.Platform$OS.type != "windows") {## No R_Outputfile (in C) on Windows

.print.mpfr <- function(x, digits = NA, ...) {
    stopifnot(is(x, "mpfr"), is.na(digits) || digits >= 2)
    ## digits = NA --> the inherent precision of x will be used
    if(length(x) >= 1)
	.Call(print_mpfr, x, as.integer(digits))
    invisible(x)
}
}# non-Windows only

## a faster version of getDataPart(.) - as we *KNOW* we have a list
## !! If ever the internal representation of such S4 objects changes, this can break !!
getD <- function(x) { attributes(x) <- NULL; x }

## Get or Set the C-global  'R_mpfr_debug_' variable:
.mpfr.debug <- function(i = NA)
    .Call(R_mpfr_set_debug, as.integer(i))

print.mpfr <- function(x, digits = NULL, drop0trailing = TRUE, ...) {
    stopifnot(is(x, "mpfr"), is.null(digits) || digits >= 2)
    ## digits = NA --> the inherent precision of x will be used
    n <- length(x)
    ch.prec <-
	if(n >= 1) {
	    rpr <- range(.getPrec(x))
	    paste("of precision ", rpr[1],
		   if(rpr[1] != rpr[2]) paste("..",rpr[2]), " bits")
	}
    cat(n, "'mpfr'", if(n == 1) "number" else "numbers", ch.prec, "\n")
    if(n >= 1)
	print(format(x, digits=digits, drop0trailing=drop0trailing), ...,
	      quote = FALSE)
    ## .Call(print_mpfr, x, as.integer(digits))
    invisible(x)
}
setMethod(show, "mpfr", function(object) print.mpfr(object))


## "[" which also keeps names ... JMC says that names are not support(ed|able)
## ---	for such objects..
.mpfr.subset <- function(x,i,j, ..., drop) {
    nA <- nargs()
    if(nA == 2) { ## x[i] etc -- vector case -- to be fast, need C! --
        xd <- structure(getD(x)[i], names=names(x)[i])
        if(any(iN <- vapply(xd, is.null, NA))) # e.g. i > length(x)
            xd[iN] <- mpfr(NA, precBits = 2L)
        ## faster than  { x@.Data <- xd ; x }:
        setDataPart(x, xd, check=FALSE)
    } else if(nA == 3 && !is.null(d <- dim(x))) { ## matrix indexing(!)
        ## not keeping dimnames though ...
        message("nargs() == 3	 'mpfr' array indexing ... ")
        new("mpfr", structure(getD(x)[i,j,...,drop=drop], dim = d))
        ## keeping dimnames: maybe try
        ##		     D <- getD(x); dim(D) <- d
        ##		     if(!is.null(dn <- dimnames(x))) dimnames(D) <- dn
        ##		     D <- D[i,,drop=drop]
        ##		     new("mpfr", D)
    }
    else
        stop(sprintf("invalid 'mpfr' subsetting (nargs = %d)",nA))
}
setMethod("[", signature(x = "mpfr", i = "ANY", j = "missing", drop = "missing"),
          .mpfr.subset)

setMethod("[[", signature(x = "mpfr", i = "ANY"),
	  function(x,i) {
	      if(length(i) > 1L) # give better error message than x@.Data[[i]] would:
		  stop("attempt to select more than one element")
	      xd <- getD(x)[[i]] # also gives error when i is "not ok"
              ## faster than { x@.Data <- list(xd) ; x }
              setDataPart(x, list(xd), check=FALSE)
	  })

## "[<-" :
.mpfr.repl <- function(x, i, ..., value, check=TRUE) {
    if(length(list(...))) ## should no longer happen:
	stop("extra replacement arguments ", deparse(list(...)),
	     " not dealt with")
    n <- length(xD <- getD(x))
    xD[i] <- value
    if((nn <- length(xD)) > n+1)
	## must "fill" the newly created NULL entries
	xD[setdiff((n+1):(nn-1), i)] <- mpfr(NA, precBits = 2L)
    setDataPart(x, xD, check=check)
}
setReplaceMethod("[", signature(x = "mpfr", i = "ANY", j = "missing",
				value = "mpfr"),
                 function(x, i, j, ..., value) .mpfr.repl(x, i, ..., value=value))
## for non-"mpfr", i.e. "ANY" 'value', coerce to mpfr with correct prec:
setReplaceMethod("[", signature(x = "mpfr", i = "missing", j = "missing",
				value = "ANY"),
	  function(x,i,j, ..., value)
		 .mpfr.repl(x, , value = mpfr(value, precBits =
				 pmax(getPrec(value), .getPrec(x)))))
setReplaceMethod("[", signature(x = "mpfr", i = "ANY", j = "missing",
				value = "ANY"),
	  function(x,i,j, ..., value)
		 .mpfr.repl(x, i, value = mpfr(value, precBits =
				  pmax(getPrec(value), .getPrec(x[i])))))


## I don't see how I could use setMethod("c", ...)
## but this works "magically"  when the first argument is an mpfr :
c.mpfr <- function(...) new("mpfr", unlist(lapply(list(...), as, Class = "mpfr")))


setMethod("unique", signature(x="mpfr", incomparables="missing"),
	  function(x, incomparables = FALSE, ...)
	  new("mpfr", unique(getD(x), incomparables, ...)))

## -> duplicated() now work

## sort() works too  (but could be made faster via faster
## ------  xtfrm() method !  [ TODO ]


setGeneric("pmin", signature = "...")# -> message about override ...
setGeneric("pmax", signature = "...")

setMethod("pmin", "Mnumber",
	  function(..., na.rm = FALSE) {
	      args <- list(...)
	      if(all(vapply(args, is.atomic, NA)))
		  return( base::pmin(..., na.rm = na.rm) )
	      ## else: at least one is "mpfr(Matrix/Array)"
	      is.m <- vapply(args, is, NA, "mpfr")
	      if(!any(is.m))
		  stop("no \"mpfr\" argument -- wrong method chosen")

	      N <- max(lengths <- vapply(args, length, 1L))
	      ## precision needed -- FIXME: should be *vector*
	      mPrec <- max(unlist(lapply(args[is.m], .getPrec)),# not vapply
			   if(any(vapply(args[!is.m], is.double, NA)))
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
	      if(all(vapply(args, is.atomic, NA)))
		  return( base::pmax(..., na.rm = na.rm) )
	      ## else: at least one is "mpfr(Matrix/Array)"
	      is.m <- vapply(args, is, NA, "mpfr")
	      if(!any(is.m))
		  stop("no \"mpfr\" argument -- wrong method chosen")

	      N <- max(lengths <- vapply(args, length, 1L))
	      ## precision needed -- FIXME: should be *vector*
	      mPrec <- max(unlist(lapply(args[is.m], .getPrec)),# not vapply
			   if(any(vapply(args[!is.m], is.double, NA)))
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
    if(h.from <- !missing(from)) {
	lf <- length(from)
	if(lf != 1) stop("'from' must be of length 1")
    }
    if ((One <- nargs() == 1L) && h.from) {
	if(is.numeric(from) || is(from,"mpfr")) {
	    to <- from; from <- mpfr(1, getPrec(from))
	} else stop("'from' is neither numeric nor \"mpfr\"")
    }
    ## else if (!is(from, "mpfr")) from <- as(from, "mpfr")

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
##	   warning("not exactly two of 'to', 'by' and 'length.out' / 'along.with' have been specified")

    if(is.null(length.out)) {
	if(!is(to,   "mpfr")) to   <- as(to,   "mpfr")
	if(!is(from, "mpfr")) from <- as(from, "mpfr")# need it again
	del <- to - from
	if(del == 0 && to == 0) return(to)
	if(missing(by)) {
	    by <- mpfr(sign(del), getD(from)[[1]]@prec)
	}
    }
    if (!is(by, "mpfr")) by <- as(by, "mpfr")
    if (length(by) != 1) stop("'by' must be of length 1")

    ## ---- This is  cut n paste  from seq.default() :
    ## ---- It should work, since "arithmetic works for mpfr :
    if(is.null(length.out)) {
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
	as(from,"mpfr")[FALSE] # of same precision
    ## else if (One) 1:length.out
    else if(missing(by)) {
	# if(from == to || length.out < 2) by <- 1
        length.out <- as.integer(length.out)
	if(missing(to))
	    to <- as(from,"mpfr") + length.out - 1
	if(missing(from))
	    from <- to - length.out + 1
	if(length.out > 2)
	    if(from == to)
		rep.int(as(from,"mpfr"), length.out)
	    else { f <- as(from,"mpfr")
		   as.vector(c(f, f + (1:(length.out - 2)) * by, to)) }
	else as.vector(c(as(from,"mpfr"), to))[seq_len(length.out)]
    }
    else if(missing(to))
	as(from,"mpfr") + (0:(as.integer(length.out) - 1L)) * by
    else if(missing(from))
	to - ((as.integer(length.out) - 1L):0) * by
    else stop("too many arguments")
}

if(FALSE) { ##-- --- I don't see *any* way  to define  seq() {S4} methods
    ## 1. Currently  need a  setGeneric() :
    ## ---- just calling setMethod("seq",...) as below fails directly {signature problem}

    ## 2. Trying three different variations --- all of them render the
    ##    *default method invalid :
    ###   --->    seq(1, length.out=3)  # afterwards fails with   " missing 'by' "
setGeneric("seq", function(from, to, by, ...) standardGeneric("seq"),
	   useAsDefault = function(from, to, by, ...)
	   base::seq(from, to, by, ...))

setGeneric("seq", function(from, to, by, ...) standardGeneric("seq"),
	   useAsDefault =
	   function(from=1, to=1, by=((to-from)/(length.out-1)), ...)
	   base::seq(from, to, by, ...))

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

}##--not yet-- defining seq() methods -- as it fails

## the fast mpfr-only version - should not return NULL
.getPrec <- function(x) {
    if(length(x)) vapply(getD(x), slot, 1L, "prec")
    else mpfr_default_prec()
}
## the user version
getPrec <- function(x, base = 10, doNumeric = TRUE, is.mpfr = NA) {
    ## if(!length(x)) ## NULL (from sapply(.) below) is not ok
    ##     return(mpfr_default_prec())
    if(isTRUE(is.mpfr) || is(x,"mpfr")) vapply(getD(x), slot, 1L, "prec")
    else if(is.character(x)) ## number of digits --> number of bits
	ceiling(log2(base) * nchar(gsub("[-.]", '', x)))
    else if(is.logical(x))
	2L # even 1 would suffice - but need 2 (in C ?)
    else if(is.raw(x))
	8L
    else {
	if(!doNumeric) stop("must specify 'precBits' for numeric 'x'")
	## else
	if(is.integer(x)) 32L
	else if(is.double(x)) 53L
	else stop(sprintf("cannot determine 'precBits' for x of type '%s'",
			  typeof(x)))
    }
}


### all.equal()

## TODO ?? <<<<<<<<<<<
## ====
## 2) instead of  as(., "mpfr")	 use  mpfr(., precBits = <smart>)

## For two "mpfr"s, use a  "smart" default tolerance :
setMethod("all.equal", signature(target = "mpfr", current = "mpfr"),
	  function (target, current,
		    tolerance =
		    2^-(0.5 * min(mean(.getPrec(target)),
				  mean(.getPrec(current)))), ...)
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
