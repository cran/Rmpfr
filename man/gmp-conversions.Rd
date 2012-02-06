\name{gmp-conversions}
\alias{.bigz2mpfr}
\alias{.mpfr2bigz}
\title{Conversion Utilities  gmp <-> Rmpfr}
\description{
  Coerce from and to big integers (\code{\link[gmp]{bigz}}) and
  \code{\link{mpfr}} numbers.
}
\usage{
.bigz2mpfr(x)
.mpfr2bigz(x, mod = NA)
}
\arguments{
  \item{x}{an \R object of class \code{bigz} or \code{mpfr} respectively.}
  \item{mod}{a possible modulus, see \code{\link[gmp]{as.bigz}} in
    package \pkg{gmp}.}
}
% \details{
% }
\value{
  a numeric vector of the same length as \code{x}, of the desired class.
}
\seealso{
  \code{\link{mpfr}()} and \code{\link[gmp]{as.bigz}} in package
  \pkg{gmp}.
}
\examples{
## Only if you have a recent enough version of R package "gmp" :
if(!is.na(r <- suppressWarnings(packageDescription("gmp",
                                                   fields="Version")))
   && package_version(r) >= 0.5) {

 S <- gmp::Stirling2(50,10)
 show(S)
 SS <- S * gmp::as.bigz(1:3)^128
 stopifnot(all.equal(log2(SS[2]) - log2(S), 128, tol=1e-15),
           identical(SS, .mpfr2bigz(.bigz2mpfr(SS))))
}% end{if}
}
\keyword{arith}% <- or better ?
