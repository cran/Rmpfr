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
if(mpfrVersion() >= "3") { ## only available since MPFR 3.0.0
  ### Airy function Ai(.)
  curve(Ai, -10, 5, n=5001); abline(h=0,v=0, col="gray", lty=3)
}

### Utilities  hypot(), atan2() : --- TODO !

## beta(), lbeta()
## ---------------
## The simplistic "slow" versions:
B  <- function(a,b) { a <- as(a, "mpfr"); b <- as(b, "mpfr"); gamma(a)*gamma(b) / gamma(a+b) }
lB <- function(a,b) { a <- as(a, "mpfr"); b <- as(b, "mpfr"); lgamma(a)+lgamma(b) - lgamma(a+b) }

x <- 1:10 + 0
(b10 <- mpfr(x, 128L))

stopifnot(all.equal(	B(1,b10),  1/x),
	  all.equal(	B(2,b10),  1/(x*(x+1))),
	  all.equal( beta(1,b10),  1/x),
	  all.equal( beta(2,b10),  1/(x*(x+1))),
	  TRUE)

for(a in c(0.1, 1, 1.5, 2, 20)) {
    c10 <- b10 + 0.25
    stopifnot(all.equal( B(a,b10), (bb <- beta(a,b10))),
	      all.equal(lB(a,b10), (lb <- lbeta(a,b10))),
	      all.equal(lb, log(bb)),
	      all.equal( B(a,c10), (bb <- beta(a,c10))),
	      all.equal(lB(a,c10), (lb <- lbeta(a,c10))),
	      all.equal(lb, log(bb)),
	      TRUE)
}

## However, the speedup is *not* much (50%) when applied to vectors:
stopifnot(validObject(xx <- outer(b10, runif(20))),
          dim(xx) == c(length(b10), 20),
          validObject(vx <- as(xx, "mpfr")), class(vx) == "mpfr", is.null(dim(vx)))
C1 <- replicate(10, system.time(bb <<- beta(vx, vx+2)))
C2 <- replicate(10, system.time(b2 <<-    B(vx, vx+2)))
summary(1000*C1[1,]) ##  80.3 {cmath-5, 2009}
summary(1000*C2[1,]) ## 125.1 { " }
stopifnot(all.equal(bb, b2))
## and for a single number, the speedup is a factor 3:
x1 <- vx[1]; x2 <- x1+2
system.time(for(i in 1:100) bb <- beta(x1, x2))# .27
system.time(for(i in 1:100) b2 <-    B(x1, x2))# .83



cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()
