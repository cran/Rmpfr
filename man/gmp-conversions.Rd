\name{gmp-conversions}
\title{Conversion Utilities gmp <-> Rmpfr}
\alias{.bigq2mpfr}
\alias{.bigz2mpfr}
\alias{.mpfr2bigz}
\alias{.mpfr2bigq}
\alias{coerce,bigq,mpfr-method}
\alias{coerce,bigz,mpfr-method}
\description{
  Coerce from and to big integers (\code{\link[gmp]{bigz}}) and
  \code{\link{mpfr}} numbers.

  Further, coerce from big rationals (\code{\link[gmp]{bigq}}) to
  \code{\link{mpfr}} numbers.
}
\usage{
.bigz2mpfr(x, precB = NULL, rnd.mode = c('N','D','U','Z','A'))
.bigq2mpfr(x, precB = NULL, rnd.mode = c('N','D','U','Z','A'))
.mpfr2bigz(x, mod = NA)
.mpfr2bigq(x)

}% TODO as(x, "mpfr")

\arguments{
  \item{x}{an \R object of class \code{bigz}, \code{bigq} or \code{mpfr}
    respectively.}
  \item{precB}{precision in bits for the result.  The default,
    \code{NULL}, means to use the \emph{minimal} precision necessary for
    correct representation.}
  \item{rnd.mode}{a 1-letter string specifying how \emph{rounding}
    should happen at C-level conversion to MPFR, see details of
    \code{\link{mpfr}}.}
  \item{mod}{a possible modulus, see \code{\link[gmp]{as.bigz}} in
    package \pkg{gmp}.}
}
\details{
  Note that we also provide the natural (S4) coercions,
  \code{as(x, "mpfr")} for \code{x} inheriting from class \code{"bigz"}
  or \code{"bigq"}.
}
\value{
  a numeric vector of the same length as \code{x}, of the desired class.
}
\seealso{
  \code{\link{mpfr}()}, \code{\link[gmp]{as.bigz}} and
  \code{\link[gmp]{as.bigq}} in package \pkg{gmp}.
}
\examples{
 S <- gmp::Stirling2(50,10)
 show(S)
 SS <- S * as.bigz(1:3)^128
 stopifnot(all.equal(log2(SS[2]) - log2(S), 128, tolerance=1e-15),
           identical(SS, .mpfr2bigz(.bigz2mpfr(SS))))

 .bigz2mpfr(S)            # 148 bit precision
 .bigz2mpfr(S, precB=256) # 256 bit

 ## rational --> mpfr:
 sq <- SS / as.bigz(2)^100
 MP <- as(sq, "mpfr")
 stopifnot(identical(MP, .bigq2mpfr(sq)),
           SS == MP * as(2, "mpfr")^100)

 ## New since 2024-01-20:   mpfr --> big rational "bigq"
 Pi <- Const("pi", 128)
 m <- Pi* 2^(-5:5)
 (m <- c(m, mpfr(2, 128)^(-5:5)))

 getDenom <- Rmpfr:::getDenom
 stopifnot(is.whole(m * (d.m <- getDenom(m))))
 stopifnot(all.equal(m, mpfr(.mpfr2bigq(m), 130), tolerance = 2^-130)) # I see even
 all.equal(m, mpfr(.mpfr2bigq(m), 130), tolerance = 0) # TRUE

 m <- m * mpfr(2, 128)^200 # quite a bit larger
 stopifnot(is.whole(m * (d.m <- getDenom(m))))
 stopifnot(all.equal(m, mpfr(.mpfr2bigq(m), 130), tolerance = 2^-130)) # I see even
 all.equal(m, mpfr(.mpfr2bigq(m), 130), tolerance = 0) # TRUE

 m2 <- m * mpfr(2, 128)^20000 ## really huge
 stopifnot(is.whole(m2 * (d.m2 <- getDenom(m2))))
 stopifnot(all.equal(m2, mpfr(.mpfr2bigq(m2), 130), tolerance = 2^-130)) # I see even
           all.equal(m2, mpfr(.mpfr2bigq(m2), 130), tolerance = 0) # TRUE
}
\keyword{arith}% <- or better ?
