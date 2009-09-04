#### All  coercion methods for the  "Rmpfr" classes


## was  toMpfr <-
mpfr <- function(x, precBits, base = 10)
{
    if(is.character(x))
	stopifnot(length(base) == 1, 2 <= base, base <= 36)
    if(missing(precBits)) {
	if(is.character(x)) ## number of digits --> number of bits
	    precBits <- ceiling(log2(base) * nchar(gsub("[-.]", '', x)))
	else if(is.logical(x))
	    precBits <- 2L # even 1 would suffice - but need 2 (in C ?)
	else if(is.raw(x))
	    precBits <- 8L
	else stop("must specify 'precBits' for numeric 'x'")
    }
    ## libmpfr would exit (after good error message) for precBits == 1
    stopifnot(precBits >= 2)
    if(is.numeric(x) || is.logical(x) || is.raw(x)) {
	new("mpfr", .Call("d2mpfr1_list", x, precBits, PACKAGE="Rmpfr"))
    } else if(is.character(x)) {
	new("mpfr", .Call("str2mpfr1_list", x, precBits, base, PACKAGE="Rmpfr"))
    }
    else stop("invalid 'x'. Must be numeric (logical, raw) or character")
}

setAs("numeric", "mpfr1", ## use default precision of 128 bits
      function(from) .Call("d2mpfr1", from, 128L, PACKAGE="Rmpfr"))
setAs("numeric", "mpfr", function(from) mpfr(from, 128L))
setAs("integer", "mpfr", function(from) mpfr(from,  32L))
setAs("raw",     "mpfr", function(from) mpfr(from,   8L))
setAs("logical", "mpfr", function(from) mpfr(from,   2L))
## TODO?  base=16 for "0x" or "0X" prefix -- but base must have length 1 ..
setAs("character", "mpfr", function(from) mpfr(from))

setAs("mpfr", "numeric", function(from) .Call("mpfr2d", from, PACKAGE="Rmpfr"))
setAs("mpfr", "integer", function(from) .Call("mpfr2i", from, PACKAGE="Rmpfr"))
setMethod("as.numeric", "mpfr", function(x) .Call("mpfr2d", x, PACKAGE="Rmpfr"))
setMethod("as.integer", "mpfr", function(x) .Call("mpfr2i", x, PACKAGE="Rmpfr"))

setAs("mpfr1", "numeric",  ## just for user-de-confusion :
      function(from) {
	  warning("coercing \"mpfr1\" via \"mpfr\" (inefficient)")
	  as(new("mpfr", list(from)), "numeric") })

setAs("mpfr1", "mpfr", function(from) new("mpfr", list(from)))
setAs("mpfr", "mpfr1", function(from) {
    if(length(from) == 1) from[[1]] else
    stop("only \"mpfr\" objects of length 1 can be coerced to \"mpfr1\"")
})


.mpfr2str <- function(x, digits = NULL) {
    stopifnot(is.null(digits) ||
	      (is.numeric(digits) && digits >= 0))
    ##	digits = NULL : use as many digits "as needed"
    .Call("mpfr2str", x, digits, PACKAGE="Rmpfr")
}

formatMpfr <-
    function(x, digits = NULL, trim = FALSE, scientific = NA,
             big.mark = "", big.interval = 3L,
             small.mark = "", small.interval = 5L, decimal.mark = ".",
	     zero.print = NULL, drop0trailing = FALSE, ...)
{
    stopifnot(is.null(digits) ||
	      (is.numeric(digits) && digits >= 0))
    ##	digits = NULL : use as many digits "as needed"

    ff <- .mpfr2str(x, digits)
    hasMinus <- sign(x) == -1
    isNum <- ff$finite	## ff$finite == is.finite(x)
    i0 <- ff$is.0	## == mpfr.is.0(x)
    ex <- ff$exp ## the *decimal* exp : one too large *unless* x == 0
    r  <- ff$str
    if(is.null(digits)) digits <- nchar(r)

    if(any(i0)) {
	Ex <- ex
	Ex[!i0] <- ex[!i0] - 1L
    } else Ex <- ex - 1L

    if(!all(isNum))
        r[!isNum] <- gsub("@", '', r[!isNum], fixed=TRUE)

    ## (maybe) add decimal point
    patch <- function(str, k)
	paste(substr   (str, 1, k),
	      substring(str, k+1), sep = decimal.mark)

    if(is.na(scientific))
        scientific <- as.numeric(getOption("scipen"))
    ## This very much depends on the desired format.
    ## if(scientific) --> all get a final "e<exp>"; otherwise, we
    ## adopt the following simple scheme :
    hasE <- { if(is.logical(scientific)) scientific else
              isNum & (Ex < -4 + scientific | Ex >= digits) }

    if(any(hasE)) {
	i. <- 1+hasMinus
	ii <- isNum & hasE
	r[ii] <- patch(r[ii], i.[ii])
	## FIXME : if this is correct, make it simpler :
	## r[hasE] <- paste(r[hasE], as.character(Ex[hasE]), sep = "e")
	r[ii] <- paste(r[ii], as.character(Ex[ii]), sep = "e")
    }
    if(!all(hasE)) { ## "non-scientific" i.e. with out final  e<nn> :
	ii <- isNum & !hasE
	iNeg <- ex <= 0	 & ii ## i.e., ex	 in {0,-1,-2,-3}
	iPos <- ex > 0	 & ii ## i.e., ex	 in {1,2..., digits}

	nZeros <- function(n)
	    sapply(n, function(k) paste(rep.int("0", k), collapse = ""))
	if(any(iNeg)) { ## "0.00..." : be careful with minus sign
	    if(any(isMin <- hasMinus[iNeg])) {
		rr <- r[iNeg]
		rr[isMin] <- substring(rr[isMin], 2)
		r[iNeg] <- paste(c("","-")[1+isMin], "0.",
				 nZeros(-ex[iNeg]), rr, sep="")
	    }
	    else {
		r[iNeg] <- paste("0.", nZeros(-ex[iNeg]), r[iNeg], sep="")
	    }
	}
	if(any(iPos)) ## "xy.nnnn" :
	    r[iPos] <- patch(r[iPos], (hasMinus + ex)[iPos])
    }
    prettyNum(r, big.mark = big.mark, big.interval = big.interval,
              small.mark = small.mark,
              small.interval = small.interval,
              decimal.mark = decimal.mark,
              zero.print = zero.print, drop0trailing = drop0trailing,
              preserve.width = if (trim) "individual" else "common")
}
setMethod("format", "mpfr", formatMpfr)

setAs("mpfr", "character", function(from)
      format(from, digits=NULL, drop0trailing = TRUE))

setAs("character", "mpfr", function(from) mpfr(from))
