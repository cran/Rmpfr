#### All  coercion methods for the  "Rmpfr" classes

if(getRversion() < "3.5") {
    isFALSE <- function (x) is.logical(x) && length(x) == 1L && !is.na(x) && !x
    isTRUE  <- function (x) is.logical(x) && length(x) == 1L && !is.na(x) && x

if(getRversion() < "3.3")
    strrep <- function (x, times) { ## (x, times) must be "recycled"
	if((lx <- length(x)) < (lt <- length(times)))
	    x <- rep_len(x, lt)
	else if(lt < lx)
	    times <- rep_len(times, lx)
	vapply(seq_along(x),
	       function(i) paste(rep.int(x[i], times[i]), collapse = ""), "")
    }
if(getRversion() < "3.2")
    lengths <- function(x, use.names = TRUE) vapply(x, length, 1L, USE.NAMES = use.names)
}

##' fast pre-test (for numeric, bigz, bigq, ..):
is.mpfr <- function(x) isS4(x) && is(x, "mpfr")

mpfr <- function(x, precBits, ...) UseMethod("mpfr")

mpfr.mpfr <- function(x, precBits, rnd.mode = c('N','D','U','Z','A'), ...)
    roundMpfr(x, precBits=precBits, rnd.mode=rnd.mode)

mpfr.bigz <- function(x, precBits, ...) {
    if(missing(precBits)) precBits <- max(2L, frexpZ(x)$exp)
    if(getOption("verbose"))
	warning("mpfr(<bigz>) --> .bigz2mpfr() [not efficiently via character]")
    ..bigz2mpfr(x, precBits)
}

mpfr.bigq <- function(x, precBits, ...) {
    if(missing(precBits)) precBits <- getPrec(x)#-> warning
    if(getOption("verbose"))
	warning("mpfr(<bigq>) --> .bigq2mpfr() [not efficiently via character]")
    ..bigq2mpfr(x, precBits)
}

mpfr.NULL <- function(x, ...) mpfr(logical(), ...)

mpfr.default <- function(x, precBits, base = 10, rnd.mode = c('N','D','U','Z','A'),
                         scientific = NA, ...)
{
    if(is.ch <- is.character(x))
	stopifnot(length(base) == 1, 2 <= base, base <= 62)
    stopifnot(is.character(rnd.mode <- toupper(rnd.mode)))
    rnd.mode <- match.arg(rnd.mode)

    if(is.raw(x)) { # is.raw() is faster
	stopifnot(missing(precBits) || precBits >= 2)
	## else warning("unrecognized raw 'x'") # <- ?? {see use in ../tests/create.R }
        ## {but 'raw' is treated below}
    } else { ## typically the result of Vectorize() or similar on "mpfr"
        if(is.list(x) && all(lengths(lc <- lapply(x, class)) == 1L) &&
           all(unlist(lc) == "mpfr1"))
        return(new("mpfr", x))
    } ## else
    if(missing(precBits)) {
	precBits <- getPrec(x, base = base, doNumeric = FALSE)
    }
    ## libmpfr would exit (after good error message) for precBits == 1
    stopifnot(precBits >= 2)

    ml <-
	if(is.numeric(x) || is.logical(x) || is.raw(x))
	    .Call(d2mpfr1_list, x, precBits, rnd.mode)
	else if(is.ch)
	    .Call(str2mpfr1_list,x, precBits, base, rnd.mode)
	else stop("invalid 'x'. Must be numeric (logical, raw) or character")
    if(is.array(x)) {
	dim <- dim(x) ; dn <- dimnames(x)
	new(if(length(dim) == 2) "mpfrMatrix" else "mpfrArray",
	    ml, Dim = dim,
	    Dimnames = if(is.null(dn)) vector("list", length(dim)) else dn)
    }
    else new("mpfr", ml)
} ## mpfr.default()

.mpfr <- function(x, precBits)
    new("mpfr", .Call(d2mpfr1_list, x, precBits, "N"))
.mpfr. <- function(x, precBits, rnd.mode)
    new("mpfr", .Call(d2mpfr1_list, x, precBits, rnd.mode))


##' to be used in our own low-level R programming
.d2mpfr1 <- function(x, precBits) .Call(d2mpfr1, x, precBits, "N")
setAs("numeric", "mpfr1", ## use default precision of 128 bits
      function(from) .Call(d2mpfr1, from, 128L, "N"))# <- round to [N]earest
setAs("numeric", "mpfr", function(from) .mpfr(from, 128L))
setAs("integer", "mpfr", function(from) .mpfr(from,  32L))
setAs("raw",     "mpfr", function(from) .mpfr(from,   8L))
setAs("logical", "mpfr", function(from) .mpfr(from,   2L))
## TODO?  base=16 for "0x" or "0X" prefix -- but base must have length 1 ..
setAs("character", "mpfr", function(from) mpfr(from))

setAs("mpfr", "numeric", function(from) .Call(mpfr2d, from, rnd.mode="N"))
setAs("mpfr", "integer", function(from) .Call(mpfr2i, from, rnd.mode="N"))
setMethod("as.numeric", "mpfr", function(x, rnd.mode="N") .Call(mpfr2d, x, rnd.mode))
## "Z": round towards [Z]ero -- crucial for as.integer() :
setMethod("as.integer", "mpfr", function(x, rnd.mode="Z") .Call(mpfr2i, x, rnd.mode))

## FIXME (in gmp!!): asNumeric() should get "..." argument
setMethod("asNumeric", "mpfr",      function(x) .Call(mpfr2d, x, rnd.mode="N"))
setMethod("asNumeric", "mpfrArray", function(x) toNum(x, rnd.mode="N"))

setAs("mpfr1", "numeric",  ## just for user-de-confusion :
      function(from) {
	  warning("coercing \"mpfr1\" via \"mpfr\" (inefficient)")
	  as(new("mpfr", list(from)), "numeric") })

setAs("mpfr1", "mpfr", function(from) new("mpfr", list(from)))
setAs("mpfr", "mpfr1", function(from) {
    if(length(from) == 1) getD(from)[[1]] else
    stop("only \"mpfr\" objects of length 1 can be coerced to \"mpfr1\"")
})

.mpfr1tolist <- function(x)
    sapply(.slotNames(x), slot, object=x, simplify=FALSE)
.mpfr2list <- function(x, names=FALSE) {
    if(isTRUE(names)) names <- format(x)
    x <- lapply(getD(x), .mpfr1tolist)
    if(is.character(names))
	names(x) <- names
    x
}


## Breaks the working of vapply(q, FUN.x) in pbetaI() in ./special-fun.R :
## as.list.mpfr1 <- function(x, ...) .mpfr1tolist(x)
## as.list.mpfr  <- function(x, ...) .mpfr2list(x)

## and then
mpfrXport <- function(x, names=FALSE) {
    if(!is.mpfr(x)) stop("argument is not a \"mpfr\" object")
    structure(class = "mpfrXport",
	      list(gmp.numb.bits = .mpfr_gmp_numbbits(),
		   ## currently unused, but in case:
		   mpfr.version	 = .mpfrVersion(),
		   Machine  = .Machine[grepl("sizeof",names(.Machine))],
		   Sys.info = Sys.info()[c("sysname", "machine")],
		   mpfr = .mpfr2list(x, names=names)))
}

mpfrImport <- function(mxp) {
    if(!inherits(mxp, "mpfrXport")) stop("need an \"mpfrXport\" object")
    nbits <- .mpfr_gmp_numbbits()
    if(!identical(nbits, mxp$gmp.numb.bits))
	stop("GMP bits not matching: 'x' has ", mxp$gmp.numb.bits,
	     "; the loaded 'Rmpfr' package has ", nbits)
    m1 <- lapply(mxp$mpfr, function(o) do.call(new, c("mpfr1", o)))
    new("mpfr", m1)
}

.mpfr2str <- function(x, digits = NULL, maybe.full = !is.null(digits), base = 10L) {
    ## digits = NULL : use as many digits "as needed" for the precision
    stopifnot(is.null(digits) || (is.numeric(digits) && length(digits) == 1 && digits >= 0),
              is.logical(maybe.full), length(maybe.full) == 1L, !is.na(maybe.full),
	      is.numeric(base),       length(base)       == 1L, base == as.integer(base),
	      2 <= base, base <= 62)
    if(!is.null(digits) && digits == 1 && base %in% 2L^(1:5)) {
	## MPFR mpfr_get_str(): "N must be >= 2"; we found that N = 1 is ok unless
	##      for these bases where it aborts (in C). ==> prevent that:
	digits <- 2L
	message(gettextf("base = %d, digits = 1 is increased to digits = 2", base))
    }
    .Call(mpfr2str, x, digits, maybe.full, base) # -> ../src/convert.c
}

##' very low level version, not exported :
..mpfr2str <- function(x, digits = NULL, maybe.full = !is.null(digits), base = 10L)
    .Call(mpfr2str, x, digits, maybe.full, base) # -> ../src/convert.c

##' more efficient, just getting the (exp, finite, is0) list, 'exp'  wrt base = 2
.mpfr_formatinfo <- function(x) .Call(R_mpfr_formatinfo, x)

##' getting the 'exp' (wrt base = 2) only  [also for extended erange!]
.mpfr2exp <- function(x) .Call(R_mpfr_2exp, x)

formatMpfr <-
    function(x, digits = NULL, trim = FALSE, scientific = NA,
	     maybe.full = !is.null(digits) && is.na(scientific),
             base = 10, showNeg0 = TRUE, max.digits = Inf,
	     big.mark = "", big.interval = 3L,
	     small.mark = "", small.interval = 5L,
             decimal.mark = ".",
             exponent.char = if(base <= 14) "e" else if(base <= 36) "E" else "|e",
             exponent.plus = TRUE,
	     zero.print = NULL, drop0trailing = FALSE, ...)
{
    ## digits = NULL : use as many digits "as needed"
    ff <- .mpfr2str(x, digits, maybe.full=maybe.full, base=base) # (checks its args!)
    ## FIXME/TODO: If have very large numbers, but not high precision, should detect it
    ## ==========  and use  maybe.full = FALSE also for the default scientific = NA
    ## digs.x <- ceiling(.getPrec(x) / log2(base))

    stopifnot(length(scientific) == 1L)
### max.digits "doomed":  scientific := number, (~= getOption("scipen")) should replace it
    stopifnot(is.numeric(max.digits), max.digits > 0)
    if(is.numeric(digits)) stopifnot(digits <= max.digits)

    isNum <- ff$finite	## ff$finite == is.finite(x)
    i0 <- ff$is.0	## == mpfrIs0(x)
    ex <- ff$exp ## the *decimal* exp (wrt given 'base' !): one too large *unless* x == 0
    r  <- ff$str
    r.dig <- nchar(r) # (in both cases, digits NULL or not)
    ## Note that r.dig[] entries may vary, notably for digits NULL when .getPrec(x) is non-constant
    if(any(Lrg <- r.dig > max.digits)) { ## now "cut down", e.g. in print() when max.digits < Inf
	r    [Lrg] <- substr(r[Lrg], 1L, max.digits)
	r.dig[Lrg] <- max.digits
    }
    if(any(i0)) {
	## sign(x) == -1 "fails" for '-0'
	hasMinus <- substr(ff$str, 1L,1L) == "-"
	if(!showNeg0 && any(iN0 <- hasMinus & i0)) {
	    ## get rid of "-" for "negative zero"
	    r[iN0] <- substring(r[iN0], 2)
	    hasMinus[iN0] <- FALSE
	}
	Ex <- ex
	Ex[!i0] <- ex[!i0] - 1L
    } else {
	Ex <- ex - 1L
	hasMinus <- sign(x) == -1
    }

    if(!all(isNum)) ## "@Inf@", "@NaN@", ...
	r[!isNum] <- gsub("@", '', r[!isNum], fixed=TRUE)

    ##' (maybe) add decimal point after position  k
    patch <- function(str, k)
	paste(substr   (str, 1L, k),
	      substring(str, k+1L), sep = decimal.mark)

    ## scipen := penalty for using "scientific", i.e., exponential format
    scipen <-
        if(is.na(scientific))
            as.numeric(getOption("scipen"))
        else if(!(is.logical(scientific) ||
                  (is.numeric(scientific) && round(scientific) == scientific)))
            stop("'scientific' must be logical or a whole number")
        else if(is.logical(scientific)) {
            if(scientific) -32L
            else max(Ex) + 64 # << penalize much
        } else ## is.numeric(scientific)  and a whole number
            scientific

    ## This very much depends on the desired format.
    ## if(scientific) --> all get a final "e<exp>"; otherwise, we
    ## adopt the following simple scheme :
### TODO: new argument   jointly = (NA | TRUE | FALSE) or just (T | F)
### ---- if(jointly) use scalar ("global") hasE and have things *align*
    ## 'hasE' is *vector* (along 'x') :
    hasE <- {
        if(isTRUE(scientific)) TRUE
	## hasE := (wF <= wE + scipen) , where (in R's format.default, which has  jointly = TRUE ):
        ##          ~~~~~~~~~~~~~~~~
        ##      wE = neg + (d > 0) + d + 4 + e (width for E format); d = mxns - 1, mxns = max_i{nsig_i}
        ##                                               e = #{digits of exponent} -1 (= 1 or 2 in R)
        ##      wF = mxsl + rgt + (rgt != 0); rgt := max_i{ (digits right of ".")_i }
        ##	     mxsl := max_i{sleft_i}; sleft_i = sign_i + (digits left of ".")_i
        else { ## scientific = (FALSE | NA | number) --- for now :
            if(is.na(scientific))
                scientific <- scipen
            isNum & (Ex < -4 + scientific | Ex > r.dig)
        }
    }
    if(aE <- any(ii <- isNum & hasE)) {
        ii <- which(ii)
	i. <- 1L + hasMinus
	r[ii] <- patch(r[ii], i.[ii])
	if(drop0trailing)
	    ## drop 0's only after decimal mark (and drop it, if immediately there)
	    r[ii] <- sub(paste0("\\", decimal.mark, "?0+$"), "", r[ii])
        chE <- if(exponent.plus)
                   sprintf("%+.0f", Ex[ii]) # "%..f": also when Ex is outside integer range!
               else  as.character(Ex[ii])
	r[ii] <- paste(r[ii], chE, sep = exponent.char)
    }
    use.prettyN <- (base <= 14 && (!aE || exponent.char == "e"))
    if(non.sci <- !all(hasE)) { ## "non-scientific" i.e. without final  e[+-]?<n>+ :
	ii <- isNum & !hasE
	## iNeg <- ex <= 0 & ii ## i.e., ex	 in {0,-1,-2,-3}
	## iPos <- ex >  0 & ii ## i.e., ex	 in {1,2..., digits}
	iNeg <- Ex <  0	 & ii ## i.e., ex	 in {0,-1,-2,-3}
	iPos <- Ex >= 0	 & ii ## i.e., ex	 in {1,2..., digits}

	if(any(eq <- (Ex == r.dig))) {
	    r[eq] <- paste0(r[eq], "0")
	    Ex[eq] <- Ex[eq] + 1L
	}
	if(any(iNeg)) { ## "0.00..." : be careful with minus sign
	    if(any(isMin <- hasMinus[iNeg])) {
		rr <- r[iNeg]
		rr[isMin] <- substring(rr[isMin], 2)
		r[iNeg] <- paste0(c("","-")[1+isMin], "0.",
				  strrep("0", -ex[iNeg]), rr)
	    }
	    else {
		r[iNeg] <- paste0("0.", strrep("0", -ex[iNeg]), r[iNeg])
	    }
	}
	if(any(iPos)) ## "xy.nnnn" :
	    r[iPos] <- patch(r[iPos], (hasMinus + Ex+1L)[iPos])
    }
    if(use.prettyN)
        r <- prettyNum(r, big.mark = big.mark, big.interval = big.interval,
                       small.mark = small.mark,
                       small.interval = small.interval,
                       decimal.mark = decimal.mark,
                       zero.print = zero.print, drop0trailing = drop0trailing,
                       preserve.width = if (trim) "individual" else "common")
    else {
	if(non.sci && drop0trailing)
	    ## drop 0's only *after* (and together with!) decimal mark:
	    r <- sub(paste0(decimal.mark, "0+$"), "", r)
	if(!missing(big.mark) || !missing(big.interval) || !missing(small.interval) ||
	    !missing(small.mark) || !missing(big.interval) || !missing(zero.print))
	    warning("with base >= 15 or 'exponent.char != \"e\", cannot use prettyNum()")
    }
    if(is.null(d <- dim(x))) r
    else array(r, dim=d, dimnames = dimnames(x))
}
setMethod("format", "mpfr", formatMpfr)

formatN.mpfr <- function(x, drop0trailing = TRUE, ...) {
    paste0(formatMpfr(x, drop0trailing=drop0trailing, ...),"_M")
}

setAs("mpfr", "character", function(from)
      format(from, digits=NULL, drop0trailing = TRUE))

setAs("character", "mpfr", function(from) mpfr(from))
