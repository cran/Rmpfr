#### Conversions   bigz  <-> mpfr   // also bigq <--> mpfr

if(packageVersion("gmp") < "0.5.8")## <-> ../NAMESPACE
    is.matrixZQ <- function(x) !is.null(attr(x, "nrow"))

## The following code is experimental, hence the "." :

### FIXME: we go via character.. which is not really efficient.
## ------  rather "should" use  MPFR Functions
## int mpfr_set_z (mpfr_t ROP, mpz_t OP, mpfr_rnd_t RND)
## int mpfr_set_q (mpfr_t ROP, mpq_t OP, mpfr_rnd_t RND)
##
## Set the value of ROP from OP, rounded toward the given direction RND.
##
## Directly in C, we'd need both Rmpfr and gmp's  C code (!)
## TODO(?:  gmp should "export" its C++ API ( -> inst/include/*.hh )
## and we should add  'LinkingTo: gmp' to DESCRIPTION and
##  then use C++ with "C" { ...} for those parts
.bigz2mpfr <- function(x, precB = NULL, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(inherits(x, "bigz"))
    ..bigz2mpfr(x, precB, rnd.mode)
}
## Fast, no-checking (and not exported) version:
..bigz2mpfr <- function(x, precB = NULL, rnd.mode = c('N','D','U','Z','A'))
    ## precB: 4 == log2(16) = log(base)
{
    stopifnot(is.character(rnd.mode <- toupper(rnd.mode)))
    rnd.mode <- match.arg(rnd.mode)
    b <- 16L
    cx <- .as.char.bigz(x, b)
    if(is.null(precB)) precB <- 4L*nchar(cx)
    if(is.matrixZQ(x))
	new("mpfrMatrix", .Call(str2mpfr1_list, cx, precB, b, rnd.mode),
	    Dim = as.integer(dim(x)))# "bigz" has no dimnames
    else
	new("mpfr", .Call(str2mpfr1_list, cx, precB, b, rnd.mode))
}
setAs("bigz", "mpfr", function(from) ..bigz2mpfr(from))


## FIXME: rather should use MPFR -- Function :
## ----   int mpfr_get_z (mpz_t ROP, mpfr_t OP, mpfr_rnd_t RND)
## Convert OP to a `mpz_t', after rounding it with respect to RND.  ....
## FIXME(2): should 'gmp' change as.bigz into an S3 generic, so this becomes S3 method?
as.bigz.mpfr <-
.mpfr2bigz <- function(x, mod=NA) {
    if(is.null(mod)) mod <- NA_integer_
    stopifnot(is.mpfr(x),
	      is.na(mod) || (length(mod) == 1L && is.numeric(mod)))
    dx <- dim(x)
### FIXME or rather  roundMpfr()  [or even round "RND" as in mpfr_get_z() above] ??
    cx <- format(trunc(x), scientific=FALSE, drop0trailing=TRUE)
    if(!is.null(dx <- dim(x))) dim(cx) <- dx ## needed?? {should *not* be, as in base R!}
    ..as.bigz(cx, mod)
}
setAs("mpfr", "bigz", function(from) .mpfr2bigz(from))


## Fast, no-checking (and not exported) version:
..bigq2mpfr <- function(x, precB = NULL, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(is.character(rnd.mode <- toupper(rnd.mode)))
    rnd.mode <- match.arg(rnd.mode)
    N <- numerator(x)
    D <- denominator(x)
    if(is.null(precB)) {
        eN <- frexpZ(N)$exp
        eD <- frexpZ(D)$exp
        precB <- pmax(128L, eN + eD + 1L) # precision of result
    }
    ..bigz2mpfr(N, precB, rnd.mode) / ..bigz2mpfr(D, precB, rnd.mode)
}

.bigq2mpfr <- function(x, precB = NULL, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(inherits(x, "bigq"))
    ..bigq2mpfr(x, precB, rnd.mode)
}
setAs("bigq", "mpfr", function(from) ..bigq2mpfr(from))



## not exported
##' @title get denominator 'd' of  m = n/d
##' @param m an mpfr number vector
##' @return the denominator 'd' (also mpfr vector)
##' @author Martin Maechler
getDenom <-  function(m) {
    ## stopifnot(is.mpfr(m))
    e <- pmax(0L, -.mpfr2exp(m)) # 2-exponents to multiply with; integer *iff* ....
    pre <- getPrec(m)
    mpfr(2, pre)^(e + pre) ## MM: it *seems* that (e + pre -1) works too ?
}

## relies on .mpfr2bigz()  above  {which has TODO s !}
.mpfr2bigq <- function(x) {
    stopifnot(is.mpfr(x))
    d <- getDenom(x)
    as.bigq(.mpfr2bigz(x*d),
            .mpfr2bigz( d ))
}

##---- Find "as small as possible" rational approximation to real number ---
## Adapted .rat() from  MASS/R/fractions.R
## Copyright (C) 1994-2005 W. N. Venables and B. D. Ripley
##
num2bigq <- function(x, cycles = 50L, max.denominator = 2^25, verbose = FALSE)
{
    n <- length(x <- as(x, "mpfr"))# precBits = 128  if(is.numeric(x))
    fin <- is.finite(x)
    a0 <- rep(0, n)
    Z1 <- as.bigz(1)
    b0 <- rep(Z1, n)
    A <- matrix(b0) # == b0
    bb <- .mpfr2bigz(x)
    r <- x - bb     # fractional part of x
    B <- matrix(bb) #  integer   part of x
    len <- 0L
    while(any(do <- fin & (r > 1/max.denominator)) && (len <- len + 1L) <= cycles) {
              a <- a0 # a[] will be in {0,1}
              b <- b0
              which <- which(do)
              a[which] <- 1
              r[which] <- 1/r[which]
              b[which] <- .mpfr2bigz(r[which]) # includes floor(.)
              r[which] <- r[which] - b[which]
### FIXME: bug in {gmp} ?  cbind(A, a, deparse.level = 0) ## adds a 0-column !!!
              A <- cbind(A, a) # is always in {0,1}
              B <- cbind(B, b) # always bigz
              if(verbose) { cat("it=", len,":  r="); print(r); cat("B =\n"); print(B) }
          }
    pq1 <- cbind(b0,     a0)
    pq  <- cbind(B[, 1], b0)
    len <- 1L
    while((len <- len + 1L) <= ncol(B)) {
        pq0 <- pq1
        pq1 <- pq
        pq <- B[, len] * pq1 + A[, len] * pq0
    }
    if(any(N <- !fin)) pq[N, 1L] <- .mpfr2bigz(x[N])
    as.bigq(pq[,1L], pq[,2L])
}

## The .rat() version  -- i.e. working with double()   ---- *not* exported; just for debugging ..
.num2bigq <- function(x, cycles = 50L, max.denominator = 2^25, verbose = FALSE)
{
    n <- length(x <- as.numeric(x))
    fin <- is.finite(x)
    a0 <- rep(0, n)
    Z1 <- 1 #N as.bigz(1)
    b0 <- rep(Z1, n)
    A <- matrix(b0) # == b0
    bb <- floor(x) #N .mpfr2bigz(x)
    r <- x - bb     # fractional part of x
    B <- matrix(bb) #  integer   part of x
    len <- 0L
    while(any(do <- fin & (r > 1/max.denominator)) && (len <- len + 1L) <= cycles) {
              a <- a0 # a[] will be in {0,1}
              b <- b0
              which <- which(do)
              a[which] <- 1
              r[which] <- 1/r[which]
              b[which] <- floor(r[which]) #N .mpfr2bigz(r[which]) # includes floor(.)
              r[which] <- r[which] - b[which]
              A <- cbind(A, a, deparse.level=0L) # is always in {0,1}
              B <- cbind(B, b, deparse.level=0L) # always bigz
              ## if(verbose) { cat("it=", len,":  r="); print(r); cat("A, B =\n"); print(A); print(B) }
              if(verbose) { cat("it=", len,":  r="); print(r); cat("B =\n"); print(B) }
          }
    pq1 <- cbind(b0,     a0, deparse.level=0L)
    pq  <- cbind(B[, 1], b0, deparse.level=0L)
    len <- 1L
    while((len <- len + 1L) <= ncol(B)) {
        pq0 <- pq1
        pq1 <- pq
        pq <- B[, len] * pq1 + A[, len] * pq0
    }
    pq[!fin, 1] <- x[!fin] #N .mpfr2bigz(x[!fin])
    pq ## list(rat = pq, x = x)
}
