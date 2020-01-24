#### Low level stuff - debugging etc
#### =========         =========

require("Rmpfr")
options(warn = 2)# warning -> error

identical3 <- function(x,y,z)	  identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d)   identical(a,b) && identical3(b,c,d)

## sane state [when re-source()ing this file]:
.mpfr_erange_set("Emin", -(2^30-1))
.mpfr_erange_set("Emax", +(2^30-1))

###----- _1_ mpfr1 , import, xport etc -----------------------------------------
i8 <- mpfr(-2:5, 32)
x4 <- mpfr(c(NA, NaN, -Inf, Inf), 32); x4 # NA -> NaN as well
stopifnot(identical3(is.na(x4), is.nan(x4), c(T,T,F,F)))

o1 <- as(x4[1], "mpfr1")
stopifnot(is(o1, "mpfr1")) # failed previously
validObject(o1)            # ditto (failed on 64-bit only)

stopifnot(
    getPrec("0xabc", base=16, doNumeric=FALSE) == 3*4,
    getPrec(  "abc", base=16, doNumeric=FALSE) == 3*4,
    getPrec("0b1001", base=2, doNumeric=FALSE) == 4,
    getPrec(  "1001", base=2, doNumeric=FALSE) == 4,
    identical3(mpfr("0b101", base= 2),
               mpfr(  "101", base= 2), mpfr(5, precBits = 3))
   ,
    identical3(mpfr("0xabc", base=16),
               mpfr(  "abc", base=16), mpfr(2748, base=16, precBits = 12))
)

## save initial (Emin, Emax) eranges  :
erangesOrig <- .mpfr_erange()

###----- _2_ Debugging, changing MPFR defaults, .. -----------------------------
##  NB: Currently mostly  *not* documented, not even .mpfr_erange()

stopifnot(Rmpfr:::.mpfr_debug() == 0 # the default level
	  ## Activate debugging level 1:
	  , Rmpfr:::.mpfr_debug(1) == 0 # the previous level
	  ## and check it :
	  , Rmpfr:::.mpfr_debug() == 1 # the current level
)

r <- mpfr(7, 100)^-1000
r
## (same as without debugging)

## where as this does print info: -- notably the very large values [3..6]:
.eranges <- function() sapply(.mpfr_erange_kinds, .mpfr_erange, USE.NAMES=FALSE)
## now, mpfr_erange() works with a *vector* of args:
.erange2 <- function() .mpfr_erange(.mpfr_erange_kinds)
## now returning *double* - which loses some precision [ending in '04' instead of '03']:
formatC(.eranges(), format="fg")
stopifnot(identical(.eranges(), .erange2()))

.mpfr_minPrec()
.mpfr_maxPrec()# debug printing shows the long integer (on 64 bit)

## Now, level 2 :
stopifnot(Rmpfr:::.mpfr_debug(2) == 1)
r
## with quite a bit of output

if(FALSE) # on Winbuilder [2019-08-08, both 32 and 64 bit]:
.mpfr_erange_set("Emax", 1073741823)

r2 <- r^100
r2
L <- r^-100000
L3 <- L^3
str(L3, internal=TRUE)
## Class 'mpfr' [package "Rmpfr"] of length 1 and precision 100
##  internally @.Data: List of 1
##  $ :Formal class 'mpfr1' [package "Rmpfr"] with 4 slots
##   .. ..@ prec: int 100
##   .. ..@ exp : int [1:2] 842206477 0
##   .. ..@ sign: int 1
##   .. ..@ d   : int [1:4] 268435456 761715680 1492345294 -1000766770
str(L3)
## lots of debugging output, then
## 1.00989692356e+253529412
##              ^^~~~~~~~~~ 10 ^ 253'529'412 that is humongous
if(!interactive()) # not seg.faulting,  but printing a *huge* line [no longer!]
  show(L3)
## segmentation fault -- randomly; 2017-06: no longer see any problem, not even with
if(FALSE) ## well, not really, definitely not interactively for now
if(interactive())
    for(i in 1:256) show(L3)
##

## quite platform dependent {valgrind ==> bug? even in mpfr/gmp/.. ?}
str(.mpfr2list(x4))
## slightly nicer ["uniformly not worse"] (still very similar) :
str(x4, internal=TRUE)
x4 ## "similar info" as .mpfr2list(.)

## Increase maximal exponent:

tools:::assertWarning(
    .mpfr_erange_set("Emax", 5e18)) # too large {FIXME why only warning and not error ??}
.mpfr_erange("Emax") # is unchanged
if(4e18 < .mpfr_erange("max.emax")) {
    .mpfr_erange_set("Emax", 4e18) # now ok:
    stopifnot(.mpfr_erange("Emax") == 4e18)
}


## revert to no debugging:
stopifnot(Rmpfr:::.mpfr_debug(0) == 2)
.mpfr_maxPrec()

L / (r2^-1000)# 1.00000....448  (could be more accurate?)

stopifnot(exprs = {
    all.equal(L, r2^-1000, tol= 1e-27) # why not more accurate?
    all.equal(log(L), -100000 * (-1000) * log(7), tol = 1e-15)
})

## Now, our experimental "transport vehicle":
stopifnot(length(rv <- c(r, r2, L)) == 3)

str(mpfrXport(rv))
str(mpfrXport(mpfr(2, 64)^(-3:3)))
str(mpfrXport(Const("pi")* 2^(-3:3)))

## and a very large one
mil <- mpfr(1025, 111)
str(mm <- mpfrXport(xx <- mil^(2^25)))
stopifnot(all.equal(log2(xx) * 2^-25, log2(mil), tol=1e-15))

## even larger -- strictly needs extended erange:
if(.mpfr_erange("min.emin") <= -2^40) {
    .mpfr_erange_set("Emin", - 2^40)
    show(xe <- 2^mpfr(-seq(1,70, by=3)*8e8, 64))
    ## used to print wrongly {because of integer overflow in .mpfr2str()$exp},
    ## with some exponents large positive
    stopifnot(exprs = {
        ! .mpfr_erange_is_int() # as 'exp's now are double
        (ee <- as.numeric(sub(".*e","", formatMpfr(xe)))) < -240e6
        (diff(ee) + 722471990) %in% 0:1
    })
} else {
    cat(sprintf(
     "Cannot set 'Emin' to -2^40 (= %g), as .mpfr_erange(\"min.emin\") is larger,
      namely %g.\n",
     - 2^40, .mpfr_erange("min.emin")))
}

## Bill Dunlap's example (with patch about convert S_alloc bug):
##               (precision increases, then decreases)
z <- c(mpfr(1,8)/19, mpfr(1,32)/19, mpfr(1,24)/19)
cbind(fz <- format(z))
stopifnot(identical(fz, rev(format(rev(z)))))
stopifnot(identical(fz, c("0.05273",
                          "0.052631578947",
                          "0.0526315793"))) # << smaller prec, again since 2019-08-09

e.xx. <- .mpfr2exp(xx)
e.z.  <- .mpfr2exp(z)

## revert to original 'erange' settings (which gives integer 'exp'):
.mpfr_erange_set("Emax", erangesOrig[["Emax"]]) # typically  2^30 - 1 = 1073741823
.mpfr_erange_set("Emin", erangesOrig[["Emin"]])

e.xx <- .mpfr2exp(xx)
e.z  <- .mpfr2exp(z)
stopifnot(exprs = {
    .mpfr_erange_is_int()
    e.xx == e.xx.
    e.xx == 335591572
    e.z  == e.z.
    e.z  == -4
    is.integer(e.xx) # but e.xx. is double
    is.integer(e.z)
})

k1 <- mpfr(  c(123, 1234, 12345, 123456), precBits=2)
(N1 <- asNumeric(k1))# 128  1024  12288  131072 -- correct
str(sk1    <- .mpfr2str(k1))
str(sk1.   <- .mpfr2str(k1, maybe.full=TRUE))
str(sk1.2  <- .mpfr2str(k1, digits=2,        base=2))
str(sk1.2F <- .mpfr2str(k1, maybe.full=TRUE, base=2))
stopifnot(exprs = {
    identical(sk1 [1:2], list(str = c("13", "10", "12", "13"), exp = 3:6))
    identical(sk1.[1:2], list(str = c("128", "1024", "12288", "131072"), exp = 3:6))
    identical(sk1.2, list(str = c("10", "10", "11", "10"),
                          exp = c( 8L,  11L,  14L,  18L),
                          finite = rep(TRUE, 4), is.0 = rep(FALSE, 4)))
    all.equal(sk1.2[2:4], .mpfr_formatinfo(k1), tol=0) # not identical(): int <-> double
    identical(formatMpfr(k1, base=2, digits=20, drop0trailing=TRUE),
              with(sk1.2, paste0(str, sapply(exp - nchar(str), strrep, x="0"))))
    identical(formatMpfr(k1, base=2, digits=2, exponent.plus=FALSE),
              c("1.0e7", "1.0e10", "1.1e13", "1.0e17"))
})
## MM: --> need_dig is fine  but is not used in the string that is returned !!

(fk1sF <- formatMpfr(k1, scientific=FALSE)) # "the bug" --- now fixed! ==> new "Bug" in new Rmpfr ????
## was "128."  "1024." "12288." "131072." , but now obeying internal precision gives
##     "1.e2"  "1.e3"  "1.e4"   "1.e5"
(fk1 <- formatMpfr(k1, digits=6))
stopifnot(exprs = {
    N1 == as.numeric(fk1)
    ## FIXME: This should change again        "1024"
    identical(format(k1, digits=3), c("128.", "1020.", "1.23e+4", "1.31e+5"))
})
##
digs <- setNames(1:6, 1:6)
## Each of these are  4 x 6  matrices
ffix <- sapply(digs, function(d) format(k1, digits = d, scientific = FALSE)) ## *not* good at all ..
## ==> need a maybe.full=TRUE   even here
ff   <- sapply(digs, function(d) format(k1, digits = d))# sci..fic = NA -- digits=1 failing for '128'
fsci <- sapply(digs, function(d) format(k1, digits = d, scientific = TRUE)) # perfect
stopifnot(exprs = {
    length(dd <- dim(ff)) == 2
    identical(dd, dim(ffix))
    identical(dd, dim(fsci))
    all.equal(asNumeric(fsci), asNumeric(ffix) -> dmat, tol=0)
    all.equal(asNumeric(ff),   asNumeric(ffix), tol=0)
})
rE <- 1 - dmat / asNumeric(k1)
i <- 1:5
summary(fm <- lm(log10(colMeans(abs(rE)))[i] ~ i))
stopifnot(exprs = {
    rE[ cbind(FALSE, upper.tri(rE)[,-6]) ] == 0
    abs(residuals(fm)) < 0.15
})

## formatting / printing :
tenth <- mpfr(-12:12, 52)/10
cents <- mpfr(-11:11, 64)/100
(kxi <- sort(c(k1, x4, i8, tenth, cents), na.last=FALSE))
mstr <- .mpfr2str       (kxi)
mfi  <- .mpfr_formatinfo(kxi)
es <- mstr$exp # base 10 ; with '0'    when  !is.finite or is0
ef <- mfi $exp # base  2 ; "undefined" when  !is.finite or is0
j2 <- c("finite", "is.0")
dxi <- cbind(x = asNumeric(kxi), prec = .getPrec(kxi),
             as.data.frame(mstr, stringsAsFactors = FALSE))
stopifnot(is.data.frame(dxi), identical(mstr$str, dxi[,"str"]),
          identical(mstr[j2], mfi[j2]),
          identical(ef, .mpfr2exp(kxi)))
dxi ## 2019-08-09: again *varying* size of 'str' rather than only growing !!
## Show that *order* no longer matters:
n <- length(ixk <- rev(kxi))
dix <- cbind(x = asNumeric(ixk), prec = .getPrec(ixk),
             as.data.frame(.mpfr2str(ixk), stringsAsFactors = FALSE))[n:1,]
attr(dix, "row.names") <- .set_row_names(n)
stopifnot(identical(dxi, dix))


## somewhat (but not so much) revealing :
cbind(prec = .getPrec(kxi), kxi = asNumeric(kxi), str = es,
      fi.10 = ceiling(ef/log2(10)), str.2 = as.integer(es*log2(10)), fi = ef)



## Bug example from RMH 2018-03-16 :
(x <- mpfr(c(65, 650, 6500, 65000, 650000), precBits=6))
data.frame(fDec = formatDec(x), f = formatMpfr(x))
x. <- as.numeric(xDec <- formatDec(x))
stopifnot(abs(x - x.) <= c(0, 0, 2, 12, 360))

cat("Checking compatibility  .mpfr_formatinfo()  <-->  .mpfr2str(*, base=2) :\n")
for(nm in ls())
    if(is(OO <- get(nm), "mpfr")) {
        cat(nm,": str(*) :\n"); str(OO); cat("compatibility: ")
        I <- .mpfr_formatinfo(OO)
        S <- .mpfr2str(OO, base = 2L)
        if(identical(I, S[-1]))
            cat("[Ok]\n")
        else {
            if(any(B <- !I$finite)) I$exp[B] <- S$exp[B]
            if(any(B <-  I $ is.0)) I$exp[B] <- S$exp[B]
            if(identical(I, S[-1]))
                cat(" after fixup [Ok]\n")
            else
                stop(".mpfr_formatinfo(*) and .mpfr2str(*, base=2) do not match")
        }
    }
