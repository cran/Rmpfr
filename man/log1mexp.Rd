\name{log1mexp}% originally ~/R/Pkgs/copula/man/log1mexp.Rd  __keep in sync !!__
\alias{log1pexp}
\alias{log1mexp}
\title{Compute  f(a) = \eqn{\mathrm{log}}{log}(1 +/- \eqn{\mathrm{exp}}{exp}(-a))
  Numerically Optimally}
\description{
  Compute f(a) = log(1 - exp(-a)), respectively
  g(x) = log(1 + exp(x))  quickly numerically accurately.
}
\usage{
log1mexp(a, cutoff = log(2))
log1pexp(x, c0 = -37, c1 = 18, c2 = 33.3)
}
\arguments{
  \item{a}{numeric (or \code{\link{mpfr}}) vector of positive values.}
  \item{x}{numeric vector, may also be an \code{"\link{mpfr}"} object.}
  \item{cutoff}{positive number; \code{log(2)} is \dQuote{optimal},
    %% see below, TODO
    but the exact value is unimportant, and anything in
    \eqn{[0.5, 1]} is fine.}
  \item{c0, c1, c2}{cutoffs for \code{log1pexp}; see below.}
}
\value{
\deqn{log1mexp(a) := f(a) = \log(1 - \exp(-a)) = \mathrm{log1p}(-\exp(-a)) = \log(-\mathrm{expm1}(-a))}{%
      log1mexp(a) := f(a) =  log(1 -  exp(-a)) = log1p(- exp(-a)) = log(- expm1(-a))}
or, respectively,

\deqn{log1pexp(x) := g(x) = \log(1 + \exp(x)) = \mathrm{log1p}(\exp(x))}{%
      log1pexp(x) := g(x) =  log(1 +  exp(x)) = log1p(exp(x))}
computed accurately and quickly.
}
\author{Martin Maechler, May 2002; \code{log1pexp()} in 2012}
\references{%  ~/R/Pkgs/Rmpfr/vignettes/log1mexp-note.Rnw
  Martin \enc{Mächler}{Maechler} (2012).
  Accurately Computing \eqn{\log(1-\exp(-|a|))};
  \url{https://CRAN.R-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}.
  %  see also <<log1mexp>> in ~/R/Pkgs/copula/inst/doc/Frank-Rmpfr.Rnw
}
\examples{
fExpr <- expression(
          DEF   = log(1 - exp(-a)),
          expm1 = log(-expm1(-a)),
          log1p = log1p(-exp(-a)),
          F     = log1mexp(a))
a. <- 2^seq(-58, 10, length = 256)
a <- a. ; str(fa <- do.call(cbind, as.list(fExpr)))
head(fa)# expm1() works here
tail(fa)# log1p() works here

## graphically:
lwd <- 1.5*(5:2); col <- adjustcolor(1:4, 0.4)
op <- par(mfcol=c(1,2), mgp = c(1.25, .6, 0), mar = .1+c(3,2,1,1))
  matplot(a, fa, type = "l", log = "x", col=col, lwd=lwd)
  legend("topleft", fExpr, col=col, lwd=lwd, lty=1:4, bty="n")
  # expm1() & log1mexp() work here

  matplot(a, -fa, type = "l", log = "xy", col=col, lwd=lwd)
  legend("left", paste("-",fExpr), col=col, lwd=lwd, lty=1:4, bty="n")
  # log1p() & log1mexp() work here
par(op)

aM <- 2^seqMpfr(-58, 10, length=length(a.)) # => default prec = 128
a <- aM; dim(faM <- do.call(cbind, as.list(fExpr))) # 256 x 4, "same" as 'fa'
## Here, for small 'a' log1p() and even 'DEF' is still good enough
l_f <- asNumeric(log(-faM))
all.equal(l_f[,"F"], l_f[,"log1p"], tol=0) # see TRUE (Lnx 64-bit)
io <- a. < 80 # for these, the differences are small
all.equal(l_f[io,"F"], l_f[io,"expm1"], tol=0) # see 6.662e-9
all.equal(l_f[io,"F"], l_f[io, "DEF" ], tol=0)
stopifnot(exprs = {
  all.equal(l_f[,"F"], l_f[,"log1p"],     tol= 1e-15)
  all.equal(l_f[io,"F"], l_f[io,"expm1"], tol= 1e-7)
  all.equal(l_f[io,"F"], l_f[io, "DEF" ], tol= 1e-7)
})
## For 128-bit prec, if we go down to 2^-130, "log1p" is no longer ok:
aM2 <- 2^seqMpfr(-130, 10, by = 1/2)
a <- aM2; fa2 <- do.call(cbind, as.list(fExpr))
head(asNumeric(fa2), 12)
tail(asNumeric(fa2), 12)

matplot(a, log(-fa2[,1:3]) -log(-fa2[,"F"]),  type="l", log="x",
        lty=1:3, lwd=2*(3:1)-1, col=adjustcolor(2:4, 1/3))
legend("top", colnames(fa2)[1:3], lty=1:3, lwd=2*(3:1)-1, col=adjustcolor(2:4, 1/3))

cols <- adjustcolor(2:4, 1/3); lwd <- 2*(3:1)-1
matplot(a, 1e-40+abs(log(-fa2[,1:3]) -log(-fa2[,"F"])),  type="o", log="xy",
        main = "log1mexp(a) -- approximation rel.errors, mpfr(*, prec=128)",
        pch=21:23, cex=.6, bg=5:7, lty=1:2, lwd=lwd, col=cols)
legend("top", colnames(fa2)[1:3], bty="n", lty=1:2, lwd=lwd, col=cols,
        pch=21:23, pt.cex=.6, pt.bg=5:7)


## -------------------------- log1pexp() [simpler] --------------------

curve(log1pexp, -10, 10, asp=1)
abline(0,1, h=0,v=0, lty=3, col="gray")

## Cutoff c1 for log1pexp() -- not often "needed":
curve(log1p(exp(x)) - log1pexp(x), 16, 20, n=2049)
## need for *some* cutoff:
x <- seq(700, 720, by=2)
cbind(x, log1p(exp(x)), log1pexp(x))

## Cutoff c2 for log1pexp():
curve((x+exp(-x)) - x, 20, 40, n=1025)
curve((x+exp(-x)) - x, 33.1, 33.5, n=1025)
}
\keyword{math}
