stopifnot(require("Rmpfr"))

n <- 1000
head(x <- mpfr(0:n, 100) / n)

stopifnot(range(x) == 0:1
	  ,all.equal(as.numeric(j0(x)),
		     besselJ(as.numeric(x), 0), tol = 1e-14)
	  ,all.equal(as.numeric(j1(x)),
		     besselJ(as.numeric(x), 1), tol = 1e-14)
	  ,all.equal(as.numeric(y0(x)),
		     besselY(as.numeric(x), 0), tol = 1e-14)
	  ,all.equal(as.numeric(y1(x)),
		     besselY(as.numeric(x), 1), tol = 1e-14)
	  )

### pnorm() -> erf() :
u <- 7*x - 2
stopifnot(all.equal(pnorm(as.numeric(u)),
		    as.numeric(pnorm(u)), tol = 1e-14))

### Riemann's Zeta function:

## -- integer arguments --
stopifnot(all(mpfr.is.0(zeta(-2*(1:100)))))

k.neg <- 2*(-100:0) - 1
Z.neg <- zeta(k.neg)
plot(k.neg, abs(as.numeric(Z.neg)), type = "l", log="y")

Pi <- Const("pi", 128L)

## confirm published value of Euler's gamma to 100 digits
pub.g <-
    paste("0.5772156649", "0153286060", "6512090082", "4024310421", "5933593992",
          "3598805767", "2348848677", "2677766467", "0936947063", "2917467495",
          sep="")

## almost =
our.g <- Const("gamma", log2(10) * 100) # 100 digits
(ff.g <- .mpfr2str(our.g))


M <- function(x) mpfr(x, 128L)
stopifnot(all.equal(zeta( 0), -1/2,      tol = 2^-100)
          , all.equal(zeta(-1), -1/M(12),  tol = 2^-100)
          , all.equal(zeta(-3),  1/M(120), tol = 2^-100)
          ## positive ones :
          , all.equal(zeta(2),  Pi^2/6,   tol = 2^-100)
          , all.equal(zeta(4),  Pi^4/90,  tol = 2^-100)
          , all.equal(zeta(6),  Pi^6/945, tol = 2^-100)
          )

### Exponential Integral Ei(.)
curve(Ei, 0,5, n=5001)

### Utilities  hypot(), atan2() :

## TODO!


cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()
