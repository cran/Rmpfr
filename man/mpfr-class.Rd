\name{mpfr-class}
\docType{class}
\title{Class "mpfr" of Multiple Precision Floating Point Numbers}
%
\alias{mpfr-class}
\alias{mpfr1-class}
%
\alias{[[,mpfr-method}
\alias{[,mpfr,ANY,missing,missing-method}
\alias{[<-,mpfr,missing,missing,ANY-method}
\alias{[<-,mpfr,ANY,missing,ANY-method}
\alias{[<-,mpfr,ANY,missing,mpfr-method}
% FIXME: many of these are just \alias{}ed, but *NOT* documented:
\alias{as.numeric,mpfr-method}
\alias{as.vector,mpfrArray-method}
\alias{as.integer,mpfr-method}
\alias{beta,ANY,mpfr-method}
\alias{beta,ANY,mpfrArray-method}
\alias{beta,mpfr,ANY-method}
\alias{beta,mpfrArray,ANY-method}
\alias{beta,mpfrArray,mpfrArray-method}
\alias{beta,mpfr,mpfr-method}
\alias{lbeta,ANY,mpfr-method}
\alias{lbeta,ANY,mpfrArray-method}
\alias{lbeta,mpfr,ANY-method}
\alias{lbeta,mpfrArray,ANY-method}
\alias{lbeta,mpfrArray,mpfrArray-method}
\alias{lbeta,mpfr,mpfr-method}
%
\alias{coerce,mpfr,character-method}
\alias{coerce,mpfr,numeric-method}
\alias{coerce,mpfr,integer-method}
\alias{coerce,mpfr1,numeric-method}
\alias{coerce,mpfr1,mpfr-method}
\alias{coerce,integer,mpfr-method}
\alias{coerce,logical,mpfr-method}
\alias{coerce,raw,mpfr-method}
\alias{coerce,numeric,mpfr-method}
\alias{coerce,numeric,mpfr1-method}
\alias{coerce,array,mpfr-method}
\alias{coerce,character,mpfr-method}
\alias{coerce,mpfr,mpfr1-method}
%
\alias{Ops,mpfr,ANY-method}
\alias{Ops,ANY,mpfr-method}
\alias{Ops,array,mpfr-method}
\alias{Ops,mpfr,array-method}
\alias{Ops,mpfr,vector-method}
\alias{Ops,vector,mpfr-method}
\alias{Arith,mpfr,array-method}
\alias{Arith,mpfr,missing-method}
\alias{Arith,mpfr,mpfr-method}
\alias{Arith,mpfr,integer-method}
\alias{Arith,mpfr,numeric-method}
\alias{Arith,integer,mpfr-method}
\alias{Arith,numeric,mpfr-method}
\alias{Arith,array,mpfr-method}
%
\alias{Compare,mpfr,mpfr-method}
\alias{Compare,mpfr,integer-method}
\alias{Compare,mpfr,numeric-method}
\alias{Compare,integer,mpfr-method}
\alias{Compare,numeric,mpfr-method}
\alias{Compare,mpfr,array-method}
\alias{Compare,array,mpfr-method}
%
\alias{Logic,mpfr,mpfr-method}
\alias{Summary,mpfr-method}
\alias{Math,mpfr-method}
\alias{Math2,mpfr-method}
\alias{abs,mpfr-method}
\alias{log,mpfr-method}
\alias{factorial,mpfr-method}
\alias{sign,mpfr-method}

\alias{format,mpfr-method}
\alias{is.finite,mpfr-method}
\alias{is.infinite,mpfr-method}
\alias{is.na,mpfr-method}
\alias{is.nan,mpfr-method}
\alias{unique,mpfr,missing-method}
%
\alias{all.equal,mpfr,mpfr-method}
\alias{all.equal,mpfr,ANY-method}
\alias{all.equal,ANY,mpfr-method}
\alias{mean,mpfr-method}

\alias{dim<-,mpfr-method}
\alias{t,mpfr-method}

\alias{\%*\%,array_or_vector,mpfr-method}
\alias{\%*\%,mpfr,array_or_vector-method}
\alias{crossprod,array_or_vector,mpfr-method}
\alias{crossprod,mpfr,array_or_vector-method}
\alias{tcrossprod,array_or_vector,mpfr-method}
\alias{tcrossprod,mpfr,array_or_vector-method}

\alias{show,mpfr-method}
\alias{show,mpfr1-method}
\alias{print.mpfr}
\alias{print.mpfr1}
%
\description{
  \code{"mpfr"} is the class of \bold{M}ultiple \bold{P}recision
  \bold{F}loatingpoint numbers with \bold{R}eliable arithmetic.

  For the high-level user, \code{"mpfr"} objects should behave
  as standard \R's \code{\link{numeric}} \emph{vectors}, just with
  prespecified (typically high) precision.

}
\section{Objects from the Class}{
  Objects are typically created by \code{\link{mpfr}(<number>, precBits)}.
}
\section{Slots}{
  Internally, \code{"mpfr"} objects just contain standard \R
  \code{\link{list}}s where each list element is of class
  \code{"mpfr1"}, representing \emph{one} MPFR number, in a structure
  with four slots, very much parallelizing the C \code{struc} in the
  \code{mpfr} C library to which the \pkg{Rmpfr} package interfaces.

  An object of class \code{"mpfr1"} has slots
  \describe{
    \item{\code{prec}:}{\code{"integer"} specifying the maxmimal
      precision in \bold{bits}.}
    \item{\code{exp}:}{\code{"integer"} specifying the base-\bold{2}
      exponent of the number.}
    \item{\code{sign}:}{\code{"integer"}, typically \code{-1} or
      \code{1}, specifying the sign (i.e. \code{\link{sign}(.)}) of the
      number.}
    \item{\code{d}:}{an \code{"integer"} vector (of 32-bit
      \dQuote{limbs}) which corresponds to the full mantissa of the
      number.}
  }% desc
}
% \section{Extends}{
% Class \code{"\linkS4class{list}"}, from data part.
% Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
% }
\section{Methods}{
  \describe{
    \item{abs}{\code{signature(x = "mpfr")}: ... }

    \item{beta}{\code{signature(a = "mpfr", b = "ANY")},}
    \item{lbeta}{\code{signature(a = "ANY", b = "mpfrArray")},}
    \item{beta}{\code{signature(a = "mpfr", b = "mpfr")}, \dots, etc:
      Compute the beta function \eqn{B(a,b)}, using high precision,
      building on internal \code{\link{gamma}} or \code{\link{lgamma}}.
      See the help for \R's base function \code{\link[base]{beta}} for
      more.}
    \item{dim<-}{\code{signature(x = "mpfr")}: Setting a dimension
      \code{\link{dim}} on an \code{"mpfr"} object makes it into an object
      of class \code{"\linkS4class{mpfrArray}"} or (more specifically)
      \code{"mpfrMatrix"} for a length-2 dimension, see their help page;
      note that \code{t(x)} (below) is a special case of this.}

    \item{Ops}{\code{signature(e1 = "mpfr", e2 = "ANY")}: ... }
    \item{Ops}{\code{signature(e1 = "ANY", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "missing")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "integer")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "numeric")}: ... }
    \item{Arith}{\code{signature(e1 = "integer", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "numeric", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "integer")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "numeric")}: ... }
    \item{Compare}{\code{signature(e1 = "integer", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "numeric", e2 = "mpfr")}: ... }

    \item{Logic}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }

    \item{Summary}{\code{signature(x = "mpfr")}: The S4
      \code{\link[methods:S4groupGeneric]{Summary}} group functions,
      \code{\link{max}}, \code{\link{min}}, \code{\link{range}},
      \code{\link{prod}}, \code{\link{sum}},
      \code{\link{any}}, and \code{\link{all}}
      are all defined for MPFR numbers.}

    \item{Math}{\code{signature(x = "mpfr")}: All the S4
      \code{\link[methods:S4groupGeneric]{Math}} group functions are
      defined, using multiple precision (MPFR) arithmetic, from
      \code{\link{getGroupMembers}("Math")}, these are (in alphabetical order):

      \code{\link{abs}}, \code{\link{sign}}, \code{\link{sqrt}},
      \code{\link{ceiling}}, \code{\link{floor}}, \code{\link{trunc}},
      \code{\link{cummax}}, \code{\link{cummin}}, \code{\link{cumprod}},
      \code{\link{cumsum}}, \code{\link{exp}}, \code{\link{expm1}},
      \code{\link{log}}, \code{\link{log10}}, \code{\link{log2}},
      \code{\link{log1p}}, \code{\link{cos}}, \code{\link{cosh}},
      \code{\link{sin}}, \code{\link{sinh}}, \code{\link{tan}},
      \code{\link{tanh}}, \code{\link{acos}}, \code{\link{acosh}},
      \code{\link{asin}}, \code{\link{asinh}}, \code{\link{atan}},
      \code{\link{atanh}}, \code{\link{gamma}}, \code{\link{lgamma}},
      \code{\link{digamma}}, and \code{\link{trigamma}}.

      Currently, \code{digamma}, and \code{trigamma} are not provided by
      the MPFR library, and hence the methods not implemented.
      \cr
      Further, the \code{cum*()} methods are \emph{not yet} implemented.}

    \item{factorial}{\code{signature(x = "mpfr")}: this will
      \code{\link{round}} the result when \code{x} is integer valued.
      Note however that \code{\link{factorialMpfr}(n)} for integer
      \code{n} is slightly more efficient, using the MPFR function
      \samp{mpfr_fac_ui}.}

    \item{Math2}{\code{signature(x = "mpfr")}: \code{\link{round}(x,
	digits)} and \code{\link{signif}(x, digits)} methods.}

    \item{as.numeric}{\code{signature(x = "mpfr")}: ... }

    \item{as.vector}{\code{signature(x = "mpfrArray")}: as for standard
      \code{\link{array}}s, this \dQuote{drops} the \code{dim} (and
      \code{dimnames}), i.e., transforms \code{x} into an \sQuote{MPFR}
      number vector, i.e., class \code{\linkS4class{mpfr}}.}

    \item{[[}{\code{signature(x = "mpfr", i = "ANY")}, and}
    \item{[}{\code{signature(x = "mpfr", i = "ANY", j = "missing", drop = "missing")}:
      subsetting aka \dQuote{indexing} happens as for numeric vectors.}

    \item{format}{\code{signature(x = "mpfr")}, further arguments
      \code{digits = NULL, scientific = NA}, etc:
      returns \code{\link{character}} vector of same length as \code{x};
      when \code{digits} is \code{NULL}, with \emph{enough} digits to
      recreate \code{x} accurately.  For details, see
      \code{\link{formatMpfr}}.}
    \item{is.finite}{\code{signature(x = "mpfr")}: ... }
    \item{is.infinite}{\code{signature(x = "mpfr")}: ... }
    \item{is.na}{\code{signature(x = "mpfr")}: ... }
    \item{is.nan}{\code{signature(x = "mpfr")}: ... }
    \item{log}{\code{signature(x = "mpfr")}: ... }
    \item{show}{\code{signature(object = "mpfr")}: ... }
    \item{sign}{\code{signature(x = "mpfr")}: ... }

    \item{all.equal}{\code{signature(target = "mpfr", current = "mpfr")},}
    \item{all.equal}{\code{signature(target = "mpfr", current = "ANY")}, and}
    \item{all.equal}{\code{signature(target = "ANY", current = "mpfr")}:
      methods for numerical (approximate) equality,
      \code{\link[base]{all.equal}} of multiple precision numbers.  Note
      that the default \code{tolerance} (argument) is taken to correspond
      to the (smaller of the two) precisions when both main arguments are
      of class \code{"mpfr"}, and hence can be considerably less than
      double precision machine epsilon \code{\link{.Machine}$double.eps}.}
    \item{coerce}{\code{signature(from = "numeric", to = "mpfr")}:
      \code{\link{as}(., "mpfr")} coercion methods are available for
      \code{\link{character}} strings, \code{\link{numeric}}, \code{\link{integer}},
      \code{\link{logical}}, and even \code{\link{raw}}.  Note however,
      that \code{\link{mpfr}(., precBits, base)} is more flexible.}
    \item{coerce}{\code{signature(from = "mpfr", to = "numeric")}: ... }
    \item{coerce}{\code{signature(from = "mpfr", to = "character")}: ... }

    \item{unique}{\code{signature(x = "mpfr")}: and}
    \item{duplicated}{\code{signature(x = "mpfr")}: just work as with numbers.}

    \item{t}{\code{signature(x = "mpfr")}: makes \code{x} into an \eqn{n
	\times 1}{n x 1} \code{\linkS4class{mpfrMatrix}}.}
  }
}
\author{Martin Maechler}

\note{Many more methods (\dQuote{functions}) automagically work for
  \code{"mpfr"} number vectors (and matrices, see the
  \code{\linkS4class{mpfrMatrix}} class doc),
  notably
  \code{\link{sort}}, \code{\link{order}}, \code{\link{quantile}},
  \code{\link{rank}}.
}

\seealso{
  The \code{"\linkS4class{mpfrMatrix}"} class, which extends the
  \code{"mpfr"} one.

  \code{\link{roundMpfr}} to \emph{change} precision of an \code{"mpfr"}
  object; \code{\link{is.whole}()} etc.

  Special mathematical functions such as some Bessel ones, e.g., \code{\link{jn}};
  further, \code{\link{zeta}(.)} \eqn{(= \zeta(.))}, \code{Ei()} etc.
  \code{\link{Bernoulli}} numbers and the Pochhammer function
  \code{\link{pochMpfr}}.
}
\examples{
## 30 digit precision
str(x <- mpfr(c(2:3, pi), prec = 30 * log2(10)))
x^2
x[1] / x[2] # 0.66666... ~ 30 digits

## indexing - as with numeric vectors
stopifnot(identical(x[2], x[[2]]),
	  ## indexing "outside" gives NA (well: "mpfr-NaN" for now):
	  is.na(x[5]),
	  ## whereas "[[" cannot index outside:
	  is(try(x[[5]]), "try-error"),
	  ## and only select *one* element:
	  is(try(x[[2:3]]), "try-error"))

## factorial() & lfactorial would work automagically via [l]gamma(),
## but factorial() additionally has an "mpfr" method which rounds
f200 <- factorial(mpfr(200, prec = 1500)) # need high prec.!
f200
as.numeric(log2(f200))# 1245.38 -- need precBits >~ 1246 for full precision

##--> see  factorialMpfr() for more such computations.
}
\keyword{classes}
