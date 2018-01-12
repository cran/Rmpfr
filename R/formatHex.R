## sprintf("%+13.13a", x) ## hex digits after the hex point = 13

## precBits: double precision = 53 = 1 + 13*4


## conversion from Hex digits to binary sequences of digits
HextoBin <- c(
 "0"="0000",
 "1"="0001",
 "2"="0010",
 "3"="0011",
 "4"="0100",
 "5"="0101",
 "6"="0110",
 "7"="0111",
 "8"="1000",
 "9"="1001",
 "A"="1010",
 "B"="1011",
 "C"="1100",
 "D"="1101",
 "E"="1110",
 "F"="1111",
 "a"="1010",
 "b"="1011",
 "c"="1100",
 "d"="1101",
 "e"="1110",
 "f"="1111")

if(FALSE) {
## the code isn't using either of these inverses.
BintoHex <- names( HextoBin[1:16])
names(BintoHex) <- HextoBin[1:16]

Bintohex <- tolower(BintoHex)
}

## RMH mentioned that  sprintfMpfr() is "parallel" to formatMpfr()
## and agreed that in principle everything should rather be based on formatMpfr(), hence
## sprintMfpr() should become unneeded  (or be *based* on formatMpfr() and renamed as basic formatFOO()
## utility for formatHex() style format -- which differs from what the MPFR lib provides (<--> our .mpfr2str())

##' @title sprintf("%a", *)-like formatting of mpfr numbers
##' @param x mpfr-number vector
##' @param bits integer (scalar) specifing the desired number of bits ("binary digits")
##' @param style 1-character string specifying
##' @return character vector of same length as \code{x}
##' @author Martin Maechler
sprintfMpfr <- function(x, bits, style = "+", expAlign=TRUE, showNeg0 = TRUE) {
    stopifnot(length(style <- as.character(style)) == 1, nchar(style) == 1,
	      style %in% c("+", " "),
	      length(bits) == 1, bits %% 1 == 0)
    hexdigits <- 1L + (bits-1L) %/% 4L ## common to both branches
### TODO: For consistency, no longer use sprintf() for bits <= 52
### ----  currently "fails", e.g., in  mpfr(formatBin(mpfr(2, 60)))
    if(bits > 52) { # <== precBits > 53
	neg <- sign(x) == -1
	ff <- .mpfr2str(x, hexdigits + 1L, maybe.full=FALSE, ## ????
                                                             base = 16)  ## need +1
        if(!showNeg0) {
            negzero <- substr(ff$str, 1L, 2L) == "-0"
            ff$str[negzero] <- substr(ff$str[negzero], 2L, 1000000L)
            ## force "-0" to "0". neg is already consistent.
        }
	isNum <- ff$finite	## ff$finite == is.finite(x)
	i0 <- ff$is.0	## == mpfrIs0(x)
        FirstDigit <- substr(ff$str, 1L, 1L)
        FirstDigit[neg] <- substr(ff$str[neg], 2L, 2L)
        BinPlace <- c("0"=0,
                      "1"=0,
                      "2"=1, "3"=1,
                      "4"=2, "5"=2, "6"=2, "7"=2,
                      "8"=3, "9"=3, "a"=3, "b"=3, "c"=3, "d"=3, "e"=3, "f"=3)
        bitMod4 <- 2^BinPlace[FirstDigit]
        x[isNum] <- x[isNum] / bitMod4[isNum] ## reduce mantissa by 2^BinPlace
        ff <- .mpfr2str(x, hexdigits + 1L, base = 16)  ## revised input value
        if(!showNeg0) # force "-0" to "0"
            ff$str[negzero] <- substr(ff$str[negzero], 2L, 1000000L)
	ex <- ff$exp ## the *decimal* value of base-2 exp : one too large *unless* x == 0
	r  <- ff$str # the mantissa, including "-" if negative
	Ex <- ex - 1L
	if(any(i0)) Ex[i0] <- ex[i0]
	if(!all(isNum)) ## "@Inf@", "@NaN@", ...
	    r[!isNum] <- gsub("@", '', r[!isNum], fixed=TRUE)
	if(any(i <- neg & isNum))
            ## r[i] <- sub("^-", "-0x", r[i])  wrongly gives e.g. "-0x.18"; want "-0x1.8"
	    r[i] <- paste0("-0x", substr(r[i], 2L, 2L), ".",
			   substring(r[i], 3L), "p")
	if(any(i <- !neg & isNum))
	    r[i] <- paste0(style, "0x", substr(r[i], 1L, 1L), ".",
			   substring(r[i], 2L), "p")
	## r[isNum] <- paste0(r[isNum], c("", "+")[1+ (isNum & (Ex >= 0))], 4*Ex)
	Exp <- 4*Ex
	Exp[!i0] <- Exp[!i0] + BinPlace[FirstDigit[!i0]]  ## increase exponent by BinPlace
	if (expAlign) {
	    Exp.format <- c("%1.1i", "%2.2i", "%3.3i")[max(1, ceiling(log10(max(abs(Exp[isNum])))))]
	    Exp[isNum] <- sprintf(Exp.format, Exp[isNum])
	}
	r[isNum] <- paste0(r[isNum], ## add "+" for positive exponents:
			   c("", "+")[1+(isNum & (Ex >= 0))][isNum], Exp[isNum])
	r
    }
    else {
	nX <- as.character(hexdigits)
	if(!showNeg0) {
	    negzero <- substr(format(x), 1L, 2L) == "-0"
	    x[negzero] <- 0
	}
	result <- sprintf(paste0("%", style, nX, ".", nX, "a"), x)
	if(any(pInf <- is.infinite(x) & x > 0))
	    result[pInf] <- sub("+", " ", result[pInf], fixed=TRUE)
	result
    }
}

##___ ../man/formatHex.Rd ___
##           ~~~~~~~~~~~~
formatHex <- function(x, precBits = min(getPrec(x)), style = "+", expAlign=TRUE) {
    if (is.numeric(x)) {
	precBits <- getPrec(x)
	x <- mpfr(x, precBits)
    }
    precB <- as.integer(precBits)
    structure(sprintfMpfr(x, bits=precB-1L, style=style, expAlign=expAlign),
	      ##---------
	      dim = dim(x), dimnames = dimnames(x),
	      base = 16L, precBits = precB, class = c("Ncharacter", "character"))
}

formatBin <- function(x, precBits = min(getPrec(x)), scientific = TRUE,
		      left.pad = "_", right.pad = left.pad, style = "+", expAlign=TRUE)
{
    H <- formatHex(x, precBits=precBits, style=style, expAlign=expAlign)
    ## bindigits is number of binary digits after the precision point
    bindigits <- attr(H, "precBits") - 1L
    ## hexdigits is the number of hex digits after the precision point
    hexdigits <- 1L + ((bindigits-1L) %/% 4L)# *must* be correct = #{pure digits between "." and "p"}
    attributes(H) <- NULL
    finite <- is.finite(x)
    H <- H[finite]
    S <- substr(H, 1L, 1L) # sign
    A <- substr(H, 4L, 4L)
    B <- substr(H, 6L, 6L+hexdigits-1L)
    ## assumes *always* an exponent "p" which is correct
    pow <- substr(H, 6L+hexdigits+1L, 1000000L)
    sB <- strsplit(B, "")
    rsB <- do.call(rbind, sB)
    hrsB <- HextoBin[rsB]
    dim(hrsB) <- dim(rsB)
    hrsBa <- apply(hrsB, 1, paste, collapse="")
    hrsBb <- substr(hrsBa, 1, bindigits)
    ## While this is a truncation,
    ## the mpfr conversion assures that
    ## only zero characters are truncated.
    if (!scientific) {
	powers <- as.integer(pow)
        Left <- -powers + max(powers, 2-precBits)
        Right <- powers - min(powers, precBits-1)
	D <- cbind(S, "0b", strrep(left.pad, Left),
		   A, hrsBb, strrep(right.pad, Right))
	D2 <- apply(D, 1, function(x) do.call(paste, list(x, collapse="")))
	ilft <- as.integer(max(Left) + min(powers)) + 4L
	res <- paste0(substr(D2,      1L,   ilft  ), ".",
		      substr(D2, ilft+1L, 1000000L))
    }
    else {
	res <- cbind(S, "0b", A, ".", hrsBb, "p", pow)
	res <- apply(res, 1, function(x) do.call(paste, list(x, collapse="")))
    }
    result <- rep("", length(x))
    result[finite] <- res
    result[!finite] <- as.numeric(x[!finite])
    structure(result, dim = dim(x), dimnames = dimnames(x),
              base = 2L, precBits = precBits, class = c("Ncharacter", "character"))
}

print.Ncharacter <- function(x, ...) {
    y <- unclass(x)
    attr(y,"base") <- NULL
    attr(y,"precBits") <- NULL
    myR <- attr(x,"base") != 10L ## formatDec() currently left-aligns [yes, this is a hack]
    ## print(y, quote=FALSE, right = myR, ...)  # protecting against multiple 'quote' and 'right'
    ## ensuring  'quote=*' and 'right=*' in  '...'  take precedence :
    pa <- c(list(...), list(quote=FALSE, right = myR))
    do.call(print, c(list(y), pa[unique(names(pa))]))
    invisible(x)
}


## RMH 2017-05-23, ~/R/MM/Pkg-ex/Rmpfr/formatDec-revised2.R :
formatDec <- function(x, precBits = min(getPrec(x)), digits=decdigits,
                      nsmall=NULL, scientific=FALSE,
                      style="+", decimalPointAlign = TRUE,
                      ...) {
    if (is.character(x)) x <- as.numeric(x)
    if (is.numeric(x)) x <- mpfr(x, precBits)
    else if (is.complex(x)) stop("complex 'x' are not supported in \"Rmpfr\" (yet)")
    decdigits <- ceiling(log(2^precBits, 10)) + 1
    chx <- format(x, digits=max(digits, decdigits), nsmall=nsmall,
                  scientific=scientific, style=style, ...)
    if (decimalPointAlign) {
	fin.x <- is.finite(x)
	chx[fin.x] <- formatAlign(chx[fin.x], ...)
    }
    structure(chx,
	      dim = dim(x),
	      dimnames = dimnames(x),
              base = 10L, precBits = precBits, class = c("Ncharacter", "character"))
}

##' Non exported utility  currently only used in  formatDec();
##' NB:   '...' here, so we can pass '...' above which may have arguments not for here
formatAlign <- function(x, leftpad=" ", rightpad=leftpad, ...) {
  if(!length(x)) return(x)
  lr <- strsplit(x, ".", fixed=TRUE)
  l <- sapply(lr, `[`, 1) ## l left
  r <- sapply(lr, `[`, 2) ## r right
  r[is.na(r)] <- ""
  nl <- nchar(l)
  nr <- nchar(r)
  ## substring() vectorizes (with 'nl'):
  l.blank <- substring(strrep(leftpad, max(nl)), 1L, max(nl) - nl)
  r.blank <- substring(strrep(rightpad,max(nr)), 1L, max(nr) - nr)
  paste0(l.blank, l, ".", r, r.blank)
}


##' currently still used in mpfr.Ncharacter(), but _not_ as a method
mpfr.Bcharacter <- function(x, precBits, scientific = NA, ...) {
    ## was scanBin()
    stopifnot(is.numeric(precBits))
    if (is.na(scientific)) ## we look for a "p" exponent..
        scientific <- any(grepl("p", x, fixed=TRUE))
    class(x) <- NULL
    if (!scientific) {
        x <- gsub("_", "0", x) ## TODO: chartr(.......)
    }
    mpfr(x, base = 2, precBits=precBits, ...)
}

## A mpfr() method for "Ncharacter"
mpfr.Ncharacter <- function(x, precBits = attr(x, "precBits"), ...) {
    class(x) <- NULL
    B <- attr(x, "base")
    if(B == 2) ## formatBin() gives very special format :
        mpfr.Bcharacter(x, precBits = precBits, ...)
    else
        mpfr(x, base = B, precBits = precBits, ...)
}
## was
## mpfr.Dcharacter <- function(x, precBits=attr(x, "bindigits")+1, ...) {
##     class(x) <- NULL
##     mpfr(gsub(" ", "", x), base = 10, precBits=precBits, ...)
## }

`[.Ncharacter` <- ## == base :: `[.listof`
    function (x, ...) structure(NextMethod("["), class = class(x))

## Don't seem to get these to work correctly (at least not easily):
## cbind.Bcharacter <- cbind.Hcharacter <-
##     function (...) structure(NextMethod("cbind"), class = class(..1))
## rbind.Bcharacter <- rbind.Hcharacter <-
##     function (...) structure(NextMethod("rbind"), class = class(..1))

## NB:  It *could* make sense to  set default     stringsAsFactors = FALSE   here..
##      but it would *not* be used when called from data.frame() which has its own default
as.data.frame.Ncharacter <- function (x, ...) {
    ## class(x) <- class(x)[class(x) != "Ncharacter"]
    ## as.data.frame(x, ...)
    NextMethod("as.data.frame")
}
