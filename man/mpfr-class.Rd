\name{mpfr-class}
\docType{class}
\title{Class "mpfr" of Multiple Precision Floating Point Numbers}
%
\alias{mpfr-class}
\alias{mpfr1-class}
%
\alias{[,mpfr,ANY,missing,missing-method}
\alias{[<-,mpfr,ANY,missing,ANY-method}
\alias{[<-,mpfr,ANY,missing,mpfr-method}
% FIXME: many of these are just \alias{}ed, but *NOT* documented:
\alias{as.numeric,mpfr-method}
\alias{as.integer,mpfr-method}
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

    \item{Summary}{\code{signature(x = "mpfr")}: ... }

    \item{Math}{\code{signature(x = "mpfr")}: ... }

    \item{Math2}{\code{signature(x = "mpfr")}: \code{\link{round}(x,
    digits)} and \code{\link{signif}(x, digits)} methods.}

    \item{as.numeric}{\code{signature(x = "mpfr")}: ... }
    \item{[}{\code{signature(x = "mpfr", i = "ANY", j = "missing", drop = "missing")}: ... }
    \item{format}{\code{signature(x = "mpfr")}: ... }
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
  \code{\link{sort}}, \code{\link{order}}, \code{\link{rank}}.
}

\seealso{
  The \code{"\linkS4class{mpfrMatrix}"} class, which extends the
  \code{"mpfr"} one.
}
\examples{
## 30 digit precision
str(x <- mpfr(c(2:3, pi), prec = 30 * log2(10)))
x^2
x[1] / x[2] # 0.66666... ~ 30 digits
}
\keyword{classes}
