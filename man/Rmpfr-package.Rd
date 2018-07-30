\name{Rmpfr-package}
\alias{Rmpfr-package}
\alias{Rmpfr}
\docType{package}
\title{R MPFR - Multiple Precision Floating-Point Reliable}
\description{
  Rmpfr provides S4 classes and methods for arithmetic
  including transcendental ("special") functions for arbitrary
  precision floating point numbers, here often called \dQuote{mpfr -
    numbers}.  To this end, it interfaces to the LGPL'ed MPFR (Multiple
  Precision Floating-Point Reliable) Library which itself is based on
  the GMP (GNU Multiple Precision) Library.
}
\details{
%  The DESCRIPTION file:
  \packageDESCRIPTION{Rmpfr}
  \packageIndices{Rmpfr}

The following (help pages) index does not really mention that we provide \emph{many}
methods for mathematical functions, including
\code{\link{gamma}}, \code{\link{digamma}}, etc, namely, all of \R's (S4)
\code{Math} group (with the only exception of \code{\link{trigamma}}),
see the list in the examples.
Additionally also \code{\link{pnorm}}, the \dQuote{error function},
and more, see the list in \code{\link{zeta}}, and
further note the first vignette (below).

%% MM: could try more systematically

\bold{\emph{Partial} index}:
\tabular{ll}{
  \code{\link{mpfr}}              \tab  Create "mpfr" Numbers (Objects) \cr
  \code{\link{mpfrArray}}         \tab  Construct "mpfrArray" almost as by \code{\link{array}()} \cr
  \code{\link{mpfr-class}}        \tab  Class "mpfr" of Multiple Precision Floating Point Numbers \cr
  \code{\link{mpfrMatrix-class}}  \tab  Classes "mpfrMatrix" and "mpfrArray" \cr
				  \tab \cr
  \code{\link{Bernoulli}}         \tab  Bernoulli Numbers in Arbitrary Precision \cr
  \code{\link{Bessel_mpfr}}       \tab  Bessel functions of Integer Order in multiple precisions \cr
  \code{\link{c.mpfr}}            \tab  MPFR Number Utilities \cr
  \code{\link{cbind}}             \tab  "mpfr" \code{...} - Methods for Functions cbind(), rbind() \cr
  \code{\link{chooseMpfr}}        \tab  Binomial Coefficients and Pochhammer Symbol aka \cr
				  \tab  Rising Factorial \cr
  \code{\link{factorialMpfr}}     \tab  Factorial 'n!'  in Arbitrary Precision \cr
  \code{\link{formatMpfr}}        \tab  Formatting MPFR (multiprecision) Numbers \cr
  \code{\link{getPrec}}           \tab  Rmpfr - Utilities for Precision Setting, Printing, etc \cr
  \code{\link{roundMpfr}}         \tab  Rounding to Binary bits, "mpfr-internally" \cr
  \code{\link{seqMpfr}}           \tab  "mpfr" Sequence Generation \cr
  \code{\link{sumBinomMpfr}}      \tab  (Alternating) Binomial Sums via Rmpfr \cr
  \code{\link{zeta}}              \tab  Special Mathematical Functions (MPFR) \cr
				  \tab \cr
  \code{\link{integrateR}}        \tab  One-Dimensional Numerical Integration - in pure R \cr
  \code{\link{unirootR}}          \tab  One Dimensional Root (Zero) Finding - in pure R \cr
  \code{\link{optimizeR}}         \tab  High Precisione One-Dimensional Optimization \cr
  \code{\link{hjkMpfr}}           \tab  Hooke-Jeeves Derivative-Free Minimization R (working for MPFR) \cr
}

Further information is available in the following vignettes:
\tabular{ll}{
  % An overview of how to use the package, including the most important functions:
  \code{Rmpfr-pkg} \tab Rmpfr (source, pdf)\cr
  \code{log1mexp-note} \tab Acccurately Computing log(1 - exp(.)) -- Assessed by Rmpfr (source, pdf)\cr
}
}%- end{details}
\author{Martin Maechler}
\references{
  MPFR (MP Floating-Point Reliable Library),
  \url{http://mpfr.org/}

  GMP (GNU Multiple Precision library),
  \url{http://gmplib.org/}

  and see the vignettes mentioned above.
}
\seealso{
  The \R package \code{\link[gmp:biginteger]{gmp}} for big integer and
  rational numbers (\code{\link[gmp]{bigrational}}) on which \pkg{Rmpfr}
  now depends.
}
\examples{
## Using  "mpfr" numbers instead of regular numbers...
n1.25 <- mpfr(5, precBits = 256)/4
n1.25

## and then "everything" just works with the desired chosen precision:hig
n1.25 ^ c(1:7, 20, 30) ## fully precise; compare with
print(1.25 ^ 30, digits=19)

exp(n1.25)

## Show all math functions which work with "MPFR" numbers (1 exception: trigamma)
getGroupMembers("Math")

## We provide *many* arithmetic, special function, and other methods:
showMethods(classes = "mpfr")
showMethods(classes = "mpfrArray")
}
\keyword{package}
