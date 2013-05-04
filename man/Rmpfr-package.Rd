\name{Rmpfr-package}
\alias{Rmpfr-package}
\alias{Rmpfr}
\docType{package}
\title{R MPFR - Multiple Precision Floating-Point Reliable}
\description{
  Rmpfr provides S4 classes and methods for arithmetic
  including transcendental ("special") functions for arbitrary
  precision floating point numbers.  To this end, it interfaces to
  the LGPL'ed MPFR (Multiple Precision Floating-Point Reliable)
  Library which itself is based on the GMP (GNU Multiple Precision) Library.
}
\details{
  \tabular{ll}{
    Package: \tab Rmpfr\cr
    %% Type: \tab Package\cr
    %% Version: \tab 0.5-0\cr
    %% Date: \tab 2012-07-07\cr
    SystemRequirements: \tab gmp (>= 4.2.3), mpfr (>= 3.0.0)\cr
                        \tab (\bold{C} (not R!) libraries; must be installed)\cr
    Depends: \tab methods, gmp (>= 0.5-2), R (>= 2.12.0)\cr
    Imports: \tab gmp, stats, utils\cr
    Suggests: \tab MASS, polynom, sfsmisc (>= 1.0-20)\cr
    SuggestNotes: \tab MASS, polynom, sfsmisc are only needed for vignette\cr
    URL: \tab http://rmpfr.r-forge.r-project.org/\cr
    License: \tab GPL (>= 2)\cr
    %% Packaged: \tab 2012-07-07 20:33:21 UTC; maechler\cr
    %% Repository: \tab CRAN\cr
    %% Date/Publication: \tab 2012-07-13 14:45:04\cr
    %% Built: \tab R 2.15.1; x86_64-unknown-linux-gnu; 2012-08-20 09:00:35 UTC; unix\cr
}

The following (help pages) index does not really mention that we provide \emph{many}
methods for mathematical functions, including
\code{\link{gamma}}, \code{\link{digamma}}, etc, namely, all of \R's (S4)
\code{Math} group (with the only exception of \code{\link{trigamma}}),
see the list in the examples.
Additionally also \code{\link{pnorm}} and the \dQuote{error function},
and more, see the list in \code{\link{zeta}}.

Further note the first vignette.
%% MM: could try more systematically


\bold{\emph{Partial} index}:
\preformatted{
mpfr                    Create "mpfr" Numbers (Objects)
mpfrArray               Construct "mpfrArray" almost as by \code{\link{array}()}
mpfr-class              Class "mpfr" of Multiple Precision Floating Point Numbers
mpfrMatrix-class        Classes "mpfrMatrix" and "mpfrArray"

Bernoulli               Bernoulli Numbers in Arbitrary Precision
Bessel_mpfr             Bessel functions of Integer Order in multiple precisions
c.mpfr                  MPFR Number Utilities
cbind                   "mpfr" \code{...} - Methods for Functions cbind(), rbind()
chooseMpfr              Binomial Coefficients and Pochhammer Symbol aka
                        Rising Factorial
factorialMpfr           Factorial 'n!'  in Arbitrary Precision
formatMpfr              Formatting MPFR (multiprecision) Numbers
getPrec                 Rmpfr - Utilities for Precision Setting, Printing, etc
hjkMpfr                 Hooke-Jeeves Derivative-Free Minimization R (working for MPFR)
roundMpfr               Rounding to Binary bits, "mpfr-internally"
seqMpfr                 "mpfr" Sequence Generation
sumBinomMpfr            (Alternating) Binomial Sums via Rmpfr
zeta                    Special Mathematical Functions (MPFR)

integrateR              One-Dimensional Numerical Integration - in pure R
optimizeR               High Precisione One-Dimensional Optimization
unirootR                One Dimensional Root (Zero) Finding - in pure R
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
