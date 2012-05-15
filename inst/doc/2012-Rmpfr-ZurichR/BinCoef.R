### R code from vignette source 'BinCoef.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=75)
try(Mlibrary(Rmpfr))
stopifnot(require("Rmpfr"))
## MM_FIXME_  add this function to  'sfsmisc' package!
capture.and.write <- function(EXPR, first, last = 2,
                              middle = NA, i.middle,
                              dotdots = " ....... ", n.dots = 2) {
    co <- capture.output(EXPR)
    writeLines(head(co, first))
    catDots <- function(M) cat(rep.int(paste(dotdots,"\n", sep=""), M), sep="")
    catDots(n.dots)
    if(is.numeric(middle)) {
        stopifnot(length(middle) == 1, middle >= 0, middle == round(middle))
        i0 <- first+2
        if(missing(i.middle)) {
            i.middle <- max(i0, length(co) %/% 2 - middle %/% 2)
        } else { ## !missing(i.middle)
            if(i.middle < i0)
                stop("'i.middle' is too small, should be at least ", i0)
        }
        writeLines(co[i.middle-1 + seq_len(middle)])
        catDots(n.dots)
    }
    writeLines(tail(co, last))
}


###################################################
### code chunk number 2: factorial-full
###################################################
noquote(sprintf("%-30.0f", factorial(24)))


###################################################
### code chunk number 3: factorial-ini
###################################################
ns <- mpfr(5:24, 120)


###################################################
### code chunk number 4: factorial-mpfr-show (eval = FALSE)
###################################################
## ns <- mpfr(5:24, 120) ; factorial(ns)


###################################################
### code chunk number 5: factorial-mpfr
###################################################
capture.and.write(factorial(ns), 5, 4)


###################################################
### code chunk number 6: chooseM-ex-fake (eval = FALSE)
###################################################
## chooseMpfr.all(n = 70)


###################################################
### code chunk number 7: chooseM-ex-do-hidden
###################################################
capture.and.write(chooseMpfr.all(n = 70), 5, middle = 4)


