require("Rmpfr")

## must take the *larger* of the two precisions:
stopifnot(format(mpfr(1, 60) / mpfr(7, 160)) ==
          "0.14285714285714285714285714285714285714285714285712")

(x <- mpfr(0:7, 100) / 7)
stopifnot( mpfr.is.0(x - x) ) # badly failed on 64-bit

## checking hexadecimal input :
stopifnot(mpfr("0xFFFFFFFFFFFFFFFFFFFF", base=16) + 1 == 2^80,
## sign(0) == 0:
          identical(sign(as(-1:1, "mpfr")), -1:1 + 0))

eps2 <- 2 * .Machine$double.eps
eps8 <- 8 * .Machine$double.eps
stopifnot(all.equal(as.numeric(x+ 1L),
                    as.numeric(x)+1L, tol = eps2),
          (3 * x)/3 <= x,
          all.equal(as.numeric(x * 2L),
                    as.numeric(x + x), tol = 0))

all.EQ <- function(x,y, tolerance = 2^-98, ...)
    all.equal(x, y, tolerance=tolerance, ...)

u <- mpfr(0:17, 128)/17
two <- mpfr(2,100)
stopifnot(all.EQ(u ^ two, u ^ 2),
          identical(u ^ 2, u ^ 2L),
          all.EQ(two ^ u, 2 ^ u),
          identical(2 ^ u, 2L ^ u),
          floor  (3*u) == floor  (3/17*(0:17)),
          ceiling(u*5) == ceiling(5/17*(0:17))
          )

i7  <- mpfr(0:7,  200)/ 7
i17 <- mpfr(0:17, 300)/17
stopifnot(all.equal(as.numeric(x+1),
		    as.numeric(x)+1),
	  all.equal(0:7,   7 * round ( i7, 25), tol = 2e-25),
	  all.equal(0:7,   7 * round ( i7, 50), tol = 2e-50),
	  all.equal(0:17, 17 * signif(i17,100), tol = 2e-100),
	  all.equal(0:17, 17 * signif(i17, 20), tol = 2e-20)
          )

## When we compute with 100 bits,
## we should compare relative errors with  2^-100 :
prettyNum(format(abs((x+pi)-pi - x) / 2^-100), drop0 = TRUE)

`%=N=%` <- function(x,y) x == y | (is.na(x) & is.na(y))
checkPmin <- function(x, nx = as(x, "numeric")) {
    stopifnot(all.equal(x, nx),
	      pmin(x, x, 1) %=N=% x, x %=N=% pmax(x, 0, x),
	      all.equal(x, pmin(x, nx, x, 1)),
	      all.equal(x, pmax(0, nx, x, round(x, 25), 0)),
	      all.equal(pmin(x, 0.75), pmin(nx, 0.75)),
	      all.equal(pmax(x, 0.25), pmax(nx, 0.25)))
}

checkPmin(x)

nx <- (0:7)/7
 x[c(2,5)] <- NA
nx[c(2,5)] <- NA
checkPmin(x, nx)

stopifnot(all.equal( round(x, 10),  round(nx, 10)),
          all.equal(signif(x, 10), signif(nx, 10)))
