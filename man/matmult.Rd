\name{matmult}
\title{(MPFR) Matrix (Vector) Multiplication}
%% Note : all the  %*%, crossprod(), and tcrossprod() methods are documented in
%% ----  ===> ./mpfrMatrix-class.Rd
%%              ^^^^^^^^^^^^^^^^^^^
\alias{matmult}
\alias{.matmult.R}% hidden, not exported
\description{
 Matrix / vector multiplication of \code{\link{mpfr}} (and \dQuote{simple}
 \code{\link{numeric}}) matrices and vectors.

 \code{ matmult (x,y, fPrec = 2)} or
 \code{crossprod(x,y, fPrec = 2)} use higher precision in underlying computations.
}
\usage{
matmult(x, y, ...)
%% .matmult.R(x,y, op = 0L, fPrec = 1, precBits = fPrec * max(getPrec(x), getPrec(y)))
}
\arguments{
  \item{x, y}{\code{\link{numeric}} or \code{\linkS4class{mpfrMatrix}}-classed \R
    objects, i.e. semantically numeric matrices or vectors.}
  \item{\dots}{arguments passed to the hidden underlying
    \code{.matmult.R()} work horse which is also underlying the
    \code{\link{\%*\%}}, \code{\link{crossprod}()}, and \code{\link{tcrossprod}()}
    methods, see the \code{\linkS4class{mpfrMatrix}} class documentation:
    \describe{
      \item{fPrec}{a multiplication factor, a positive number determining
	the number of bits \code{precBits} used for the underlying
	multiplication and summation arithmetic.  The default is
	\code{fPrec = 1}.  Setting \code{fPrec = 2} doubles the
	precision which has been recommended, e.g., by John Nash.}
      \item{precBits}{the number of bits used for the underlying
	multiplication and summation arithmetic; by default
	\code{precBits = fPrec * max(getPrec(x), getPrec(y))} which
	typically uses the same accuracy as regular
	\code{\link{mpfr}}-arithmetic would use.}
    }}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  a (base \R) \code{\link{matrix}} or \code{\linkS4class{mpfrMatrix}},
  depending on the classes of \code{x} and \code{y}.
}
%% \references{
%% }
\author{Martin Maechler}
\note{
  Using \code{matmult(x,y)} instead of \code{x \link{\%*\%} y}, makes sense
  mainly \emph{if} you use non-default \code{fPrec} or \code{precBits} arguments.

  The \code{\link{crossprod}()}, and \code{\link{tcrossprod}()} function
  have the \emph{identical} optional arguments \code{fPrec} or \code{precBits}.
}
\seealso{
 \code{\link{\%*\%}},
 \code{\link{crossprod}}, \code{\link{tcrossprod}}.
}
\examples{
%% ## The function is currently defined as
%% function (x, y, ...)
%% .matmult.R(x, y, op = 0L, ...)

## FIXME: add example

## 1)  matmult()  <-->  \%*\%

## 2)  crossprod() , tcrossprod()  %% <--> ./mpfrMatrix-class.Rd  examples  (!)
}
\keyword{arith}


