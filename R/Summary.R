#### Define mpfr methods for Summary  group functions
####			     =======

### "Math" are done in ./Math.R , "Ops", "Arith", "Logic", "Compare" in ./Arith.R

.Summary.codes <-
    c("max" = 1, "min" = 2, "range" = 3, "prod" = 4, "sum" = 5,
      "any" = 10, "all" = 11)
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


stats__quantile.default <- stats:::quantile.default

setMethod("quantile", "mpfr", stats__quantile.default)
## FIXME: is *slow*  *and* uses double precision epsilon internally .Machine$double.epsilon
##
## Not perfect: has the "0%" "25%" "50%" ... names but not as slot ... hmm ...
## 'mpfr' numbers do not have 'names' slot ... (etc) -- but "work" with names
	  ## function(x, ...) {
	  ##     if((match("names", names(list(...)), nomatch = 0L)) == 0L)
	  ##         stats__quantile.default(x, ..., names=FALSE)
	  ##     else ## ... contains 'names = ..'
	  ##         stats__quantile.default(x, ...)
	  ## })

setMethod("mean", "mpfr", function(x, trim = 0, na.rm = FALSE, ...) {
    if(trim == 0) ## based on sum() :
	sum(x, na.rm=na.rm, ...) / length(x)
    else {
	## cut'n'paste from  mean.default() :
	if (!is.numeric(trim) || length(trim) != 1L || trim < 0)
	    stop("'trim' must be numeric of length one, in  [0, 1/2]")
	if (na.rm)
	    x <- x[!is.na(x)]
	n <- length(x)
	if (anyNA(x))
	    mpfr(NA)
	else if (trim >= 0.5)
	    quantile(x, probs = 0.5, na.rm = FALSE, names = FALSE)
	else {
	    lo <- floor(n * trim) + 1
	    hi <- n + 1 - lo
	    mean(sort(x, partial = unique(c(lo, hi)))[lo:hi], na.rm = FALSE)
	}
    }
})

setMethod("median", "mpfr",
	  function(x, na.rm=FALSE, ...)
	      quantile(x, probs = 0.5, na.rm=na.rm, names = FALSE))


setMethod("summary", "mpfr", function (object, ..., digits, quantile.type = 7) {
    ## Should work almost like  summary(asNumeric(object, ..., digits=digits))
    ## *but* w/o underflow and overflow:
    nas <- is.na(object)
    object <- object[!nas]
    qq <- quantile(object, names=FALSE, type = quantile.type)
    qq <- c(qq[1L:3L], mean(object), qq[4L:5L])
    ## names(qq) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    if (!missing(digits))
        qq <- signif(qq, digits)
    if(any(nas)) # names() updatingn works for "mpfr"
        qq  <- c(qq, "NA's" = sum(nas))
    ## loses names: as(qq, "summaryMpfr")
    ## workaround :
    new("summaryMpfr", qq,
        names = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))
})

setClass("summaryMpfr", contains = "mpfr",
         slots = c(names = "character"))

print.summaryMpfr <- function (x, digits=max(3L, getOption("digits") - 3L), ...)
{
    xx <- x
    names(xx) <- NULL # will be lost anyway
    if(getRversion() >= "3.5.2") { ## for zapsmall() to work
        finite <- is.finite(xx)
        xx[finite] <- zapsmall(xx[finite])
    }
    m <- match("NA's", names(xx), nomatch = 0L)
    xx <-
        if(m)
            c(format(xx[-m], digits = digits), `NA's` = as.character(xx[m]))
        else
            format(xx, digits = digits)
    names(xx) <- names(x)
    print.table(xx, digits = digits, ...)
    invisible(x)
}

setMethod(show, "summaryMpfr", function(object) print.summaryMpfr(object))


## FIXME: can do this considerably faster in C: [which.max(): loc.(first TRUE)]
setMethod("which.max", "mpfr", function(x) which.max(x == max(x)))
setMethod("which.min", "mpfr", function(x) which.max(x == min(x)))
