% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/Rmpfr/inst/NEWS.Rd"))
\name{NEWS}
\title{Rmpfr News}
\encoding{UTF-8}

\section{Changes in version 0.7-2 [2019-01-14, r298]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{str(<mpfr>)} no longer calls \code{formatMpfr(x, digits, *)}
      with a \code{digits} \emph{vector} of the same length as \code{x}
      (which never worked correctly).

      \item \code{seqMpfr(1, length.out=8)} now works correctly.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{unirootR()} gets an option to \emph{not} warn on
      non-convergence.

      \item Provide a \code{summary()} method for \code{"mpfr"} numbers
      closely modeled after \code{summary.default} for numeric.

      \item \code{mpfr(NULL)} now works, returning \code{mpfr(logical())}.

      \item a simple \code{sapplyMpfr()} function, showing how to work
      around the fact that \code{sapply()} does typically not work with
      \code{"mpfr"} numbers.
    }
  }
}

\section{Changes in version 0.7-1 [2018-07-24, r291]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{formatMpfr()} for large low-precision numbers now uses
      scientific representation, fixing the bug RMH was reporting March
      17 already.

      \item \code{outer()} is \dQuote{imported} from base, so it
      behaves as an \pkg{Rmpfr} function which dispatches e.g., when
      calling \code{\link{tcrossprod}()}.
    }
  }
}

\section{Changes in version 0.7-0 [2018-01-12, r284]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{.mpfr2list()} and \code{mpfrXport()} gain an option
      \code{names} (for nicer output).

      \item \code{formatMpfr()} and the \code{print()} method get a new
      option \code{max.digits} with default \code{9999} for the print
      methods, to limit the number of digits printed for high precision
      numbers.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item For non-\dQuote{regular} mpfr numbers, the \code{d} slot in
      the \code{"mpfr1"} representation is now empty instead of
      \dQuote{random}.   This also eliminates valgrind warnings about
      uninitialized values in C.
    }
  }
}
\section{Changes in version 0.6-2 [2017-05-29, r264], never on CRAN}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The S3 classes \code{"Hcharacter"} and \code{"Bcharacter"}
      resulting from \code{formatHex()} and \code{formatBin()} now
      \dQuote{inherit from} \code{"character"} formally.

      \item They also got a \code{`[`} method, so subsetting should
      work, including for \code{array}s of these.

      \item The \code{"mpfr"} method of \code{str()} gains option
      \code{internal}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item when \code{print()}ing mpfr numbers, the result no longer sometimes
      loses the last digit.

      \item \code{dnorm()} now works correctly with mpfr numbers;
      similarly \code{dbinom()} and \code{dpois()} should work in all
      cases, now.

      \item in \file{NAMESPACE}, also \code{exportMethods(apply)}, so
      \pkg{SNscan} works.

      \item \code{print(formatHex(..))}, \code{formatBin()} and
      \code{formatDec()} now look better and are more correct; the first
      two get a new option \code{expAlign} indicating to use the same
      number of digits for exponents (in \dQuote{scientific} cases).

      Notably, \code{mpfr(formatBin(mpx))} works for more \code{mpx}
      objects (of class \code{"mpfr"}).

      \item \code{format(mpfr(3,7), digits = 1, base = 2)} no longer
      crashes (from inside MPFR).

      \item \code{formatDec(mpfr(NA, 7))} now works.

      \item For non-\dQuote{regular} mpfr numbers, the \code{d} slot in
      the \code{"mpfr1"} representation is now empty instead of
      \dQuote{random}.   This also eliminates valgrind warnings about
      uninitialized values in C.
    }
  }
}

\section{Changes in version 0.6-1 [2016-11-15, r249]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{head()} and \code{tail()} methods for \code{"mpfrMatrix"}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item C-level \code{mpfr2str()} no longer calls S_realloc() with
      wrong "old size" (thanks to Bill Dunlap).
      \item \code{c()} now also works when its result is a length-0
      \code{"mpfr"} object.
    }
  }
}

\section{Changes in version 0.6-0 [2015-12-04, r237]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{mpfr()} now is S3 generic with several methods, notably a
      \code{"mpfr"} method allowing to change precision or rounding mode.

      \item \code{mpfr()}, \code{formatMpfr()}, etc, now work with bases from 2
      to 62 (using digits, upper and lower case ASCII letters, \code{62
      == 10 + 2*26} characters), as this has been in MPFR since version
      3.0.0 (see \code{\link{mpfrVersion}}), which is hence (implicitly)
      required for \code{base} greater than 36.

      \item \code{formatMpfr()} gets a new argument \code{base = 10} and can
      be used to produce in other bases, notably binary (\code{base = 2})
      or hexadecimal (\code{base = 16}).

      \item \code{str(<mpfr>, ....)} is now based on \code{formatMpfr()}
      and nicely shows numbers also out of the double precision range.

      Further, it now chooses a smart default for optional argument
      \code{vec.len}.

      \item \code{matrix(mp, ..)} now also works when \code{mp} is of
      class \code{"mpfr"}.
      \item new matrix \code{norm()} for several \code{kind}s.

      \item new functions \code{formatHex()} and \code{formatBin()}
      thanks to Rich Heiberger.

      \item \code{mpfr(x)} also works as \emph{inverse} of
      \code{formatBin} and \code{formatHex}.

      \item \code{roundMpfr()} and mathematical functions such as
      \code{jn}, or \code{chooseMpfr()} get new optional argument
      \code{rnd.mode} passed to the corresponding MPFR function.

      \item \code{median(x)}, \code{mean(x, trim)} for \code{trim > 0}
      now work fine for \code{"mpfr"} x, and \code{quantile(x, *)} no
      longer needs \code{names=FALSE} to avoid a warning.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{pnorm(.)}, \code{j0()} and similar special functions
      now preserve \code{mpfrMatrix} and \code{mpfrArray} classes.

      \item similarly, \code{is.finite()} etc keep the
      \code{dim()}ensionality for \code{"mpfrArray"} arguments.

      \item \code{mpfr("0xabc", base=16)} and \code{mpfr("0b101",
	base=2)} and corresponding \code{getPrec()} now give the correct
      precBits instead of too many.

      \item \code{str(<0-length mpfr>)} now works correctly.
    }
  }
}



\section{Changes in version 0.5-7 [2014-11-27, r205]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{as.integer()} now rounds \dQuote{to zero} as for
      regular \R numbers (it accidentally did round \dQuote{to nearest}
      previously).
    }
  }
}

\section{Changes in version 0.5-6 [2014-09-05, r203]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item experimental \code{mpfrImport()}, \code{mpfrXport()}
      utilities -- as we found cases, where save() \code{"mpfr"} objects
      were \emph{not} portable between different platforms.

      \item \code{as(*,"mpfr")} now also supports rounding mode
      \code{"A"} (\dQuote{\bold{A}way from zero}).

      \item Several hidden low level utilities also get a
      \code{rnd.mode} option.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 0.5-5 [2014-06-19, r190]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The result of \code{integrateR()} now prints even if a
      warning happened.

      \item \code{pbetaI(x, a,b)}, the arbitrarily accurate
      \code{pbeta()} computation for \emph{integer} \eqn{a} and \eqn{b},
      now works for larger \eqn{(a,b)}.

      \item Newly providing \code{mpfr}-versions of \code{dbinom()},
      \code{dpois()}, and \code{dnorm()}.

      \item New utility functions \code{mpfr_default_prec()},
      \code{.mpfr.minPrec()}, etc, to get, check, set default exponent
      ranges and precision.

      \item New \code{sinpi()}, \code{cospi()} etc, notably for \R >= 3.0.1.
    }
  }
  %% \subsection{BUG FIXES}{
  %%   \itemize{
  %%     \item .
  %%   }
  %% }
}

\section{Changes in version 0.5-4 [2013-10-22, r173]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}
