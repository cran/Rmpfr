\name{chooseMpfr}
\alias{chooseMpfr}
\alias{chooseMpfr.all}
\alias{pochMpfr}
\title{Binomial Coefficients and Pochhammer Symbol aka Rising Factorial}
\description{
  Compute binomial coefficients, \code{chooseMpfr(a,n)} being
  mathematically the same as \code{\link{choose}(a,n)}, but using high
  precision (MPFR) arithmetic.

  \code{chooseMpfr.all(n)} means the vector \code{\link{choose}(n, 1:n)},
  using enough bits for exact computation via MPFR.
  However, \code{chooseMpfr.all()} is now \bold{deprecated} in favor of
  \code{\link[gmp]{chooseZ}} from package \pkg{gmp}, as that is now
  vectorized.

  \code{pochMpfr()} computes the Pochhammer symbol or \dQuote{rising
  factorial}, also called the \dQuote{Pochhammer function},
  \dQuote{Pochhammer polynomial}, \dQuote{ascending factorial},
  \dQuote{rising sequential product} or \dQuote{upper factorial},

  \deqn{x^{(n)}=x(x+1)(x+2)\cdots(x+n-1)= \frac{(x+n-1)!}{(x-1)!} = \frac{\Gamma(x+n)}{\Gamma(x)}.
  }{x^(n) = x(x+1)(x+2)...(x+n-1) = (x+n-1)! / (x-1)! = Gamma(x+n) / Gamma(x).}
}
\usage{
chooseMpfr (a, n, rnd.mode = c("N","D","U","Z","A"))
chooseMpfr.all(n, precBits=NULL, k0=1, alternating=FALSE)
pochMpfr(a, n, rnd.mode = c("N","D","U","Z","A"))
}
\arguments{
  \item{a}{a numeric or \code{\linkS4class{mpfr}} vector.}
  \item{n}{an \code{\link{integer}} vector; if not of length one,
    \code{n} and \code{a} are recycled to the same length.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see \code{\link{mpfr}}.}
  \item{precBits}{integer or NULL for increasing the default precision
    of the result.}
  \item{k0}{integer scalar}
  \item{alternating}{logical, for \code{chooseMpfr.all()}, indicating if
    \emph{alternating sign} coefficients should be returned, see below.}
}
\note{
  Currently this works via a (C level) \code{for(i in 1:n)}-loop which
  really slow for large \code{n}, say \eqn{10^6}, with computational cost
  \eqn{O(n^2)}.  In such cases, if you need high precision
  \code{choose(a,n)} (or Pochhammer(a,n)) for large \code{n}, preferably
  work with the corresponding \code{\link{factorial}(mpfr(..))}, or
  \code{\link{gamma}(mpfr(..))} terms.
}
\value{
  For
  \describe{
    \item{\code{chooseMpfr()}, \code{pochMpfr()}:}{an
      \code{\linkS4class{mpfr}} vector of length \code{max(length(a),
	length(n))};}
    \item{\code{chooseMpfr.all(n, k0)}:}{a \code{mpfr} vector of length
      \code{n-k0+1}, of binomial coefficients \eqn{C_{n,m}}{C[n,m]}
      or, if \code{alternating} is true, \eqn{(-1)^m\cdot C_{n,m}}{(-1)^m * C[n,m]}
      for \eqn{m \in}{m in} \code{k0:n}.}
  }
}
\seealso{
  \code{\link{choose}(n,m)} (\pkg{base} \R) computes the binomial coefficient
  \eqn{C_{n,m}}{C[n,m]} which can also be expressed via Pochhammer
  symbol as
  \eqn{C_{n,m} = (n-m+1)^{(m)}/m!}{C[n,m] = (n-m+1)^(m) / m!}.

  \code{\link[gmp]{chooseZ}} from package \pkg{gmp};
  for now, %% as we "deprecate" it, in favor of  factorialZ() :
  \code{\link{factorialMpfr}}.

  For (alternating) binomial sums, directly use
  \code{\link{sumBinomMpfr}}, as that is potentially % once we do it in C
  more efficient.
}
\examples{
pochMpfr(100, 4) == 100*101*102*103 # TRUE
a <- 100:110
pochMpfr(a, 10) # exact (but too high precision)
x <- mpfr(a, 70)# should be enough
(px <- pochMpfr(x, 10)) # the same as above (needing only 70 bits)
stopifnot(pochMpfr(a, 10) == px,
          px[1] ==prod(mpfr(100:109, 100)))# used to fail

(c1 <- chooseMpfr(1000:997, 60)) # -> automatic "correct" precision
stopifnot(all.equal(c1, choose(1000:997, 60), tolerance=1e-12))

## --- Experimenting & Checking
n.set <- c(1:10, 20, 50:55, 100:105, 200:203, 300:303, 500:503,
           699:702, 999:1001)
if(!Rmpfr:::doExtras()) { ## speed up: smaller set
  n. <- n.set[-(1:10)]
  n.set <- c(1:10, n.[ c(TRUE, diff(n.) > 1)])
}
C1 <- C2 <- numeric(length(n.set))
for(i.n in seq_along(n.set)) {
  cat(n <- n.set[i.n],":")
  C1[i.n] <- system.time(c.c <- chooseMpfr.all(n) )[1]
  C2[i.n] <- system.time(c.2 <- chooseMpfr(n, 1:n))[1]
  stopifnot(is.whole(c.c), c.c == c.2,
            if(n > 60) TRUE else all.equal(c.c, choose(n, 1:n), tolerance = 1e-15))
  cat(" [Ok]\n")
}
matplot(n.set, cbind(C1,C2), type="b", log="xy",
        xlab = "n", ylab = "system.time(.)  [s]")
legend("topleft", c("chooseMpfr.all(n)", "chooseMpfr(n, 1:n)"),
       pch=as.character(1:2), col=1:2, lty=1:2, bty="n")

## Currently, chooseMpfr.all() is faster only for large n (>= 300)
## That would change if we used C-code for the *.all() version

## If you want to measure more:
measureMore <- TRUE
measureMore <- FALSE
if(measureMore) { ## takes ~ 2 minutes (on "lynne", Intel i7-7700T, ~2019)
  n.s <- 2^(5:20)
  r <- lapply(n.s, function(n) {
      N <- ceiling(10000/n)
      cat(sprintf("n=\%9g => N=\%d: ",n,N))
      ct <- system.time(C <- replicate(N, chooseMpfr(n, n/2)))
      cat("[Ok]\n")
      list(C=C, ct=ct/N)
  })
  print(ct.n <- t(sapply(r, `[[`, "ct")))
  hasSfS <- requireNamespace("sfsmisc")
  plot(ct.n[,"user.self"] ~ n.s, xlab=quote(n), ylab="system.time(.) [s]",
       main = "CPU Time for  chooseMpfr(n, n/2)",
       log ="xy", type = "b", axes = !hasSfS)
  if(hasSfS) for(side in 1:2) sfsmisc::eaxis(side)
  summary(fm <- lm(log(ct.n[,"user.self"]) ~ log(n.s), subset = n.s >= 10^4))
  ## --> slope ~= 2  ==>  It's  O(n^2)
  nn <- 2^seq(11,21, by=1/16) ; Lcol <- adjustcolor(2, 1/2)
  bet <- coef(fm)
  lines(nn, exp(predict(fm, list(n.s = nn))), col=Lcol, lwd=3)
  text(500000,1, substitute(AA \%*\% n^EE,
                            list(AA = signif(exp(bet[1]),3),
                                 EE = signif(    bet[2], 3))), col=2)
} # measure more
}% example

\keyword{arith}
