## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_Rmpfr_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(sprintf("C code of R package 'Rmpfr': GMP using %d bits per limb\n",
				  .mpfr_gmp_numbbits()))
}

.onLoad <- function(libname, pkgname) {
    if(mpfrVersion() < "3.0.0")
	warning("MPFR C library version ", format(mpfrVersion()),
		" is outdated, and minor functionality will be missing.\n",
		"  Consider installing a newer version of MPFR (e.g., from mpfr.org),\n",
		"  and re-install the R package Rmpfr after that.", call.=FALSE)
}


if(packageVersion("gmp") < "0.6-1") local({ ## need c_bigz() and c_bigq() already now
    env <- asNamespace("gmp")
    getGmp <- function(x) get(x, envir=env, inherits=FALSE)
    biginteger_c  <- getGmp("biginteger_c")
    bigrational_c <- getGmp("bigrational_c")
    rm(env, getGmp)
    c_bigz <<- function(L) .Call(biginteger_c, L)
    c_bigq <<- function(L) .Call(bigrational_c, L)
})

if(getRversion() < "4.0") {
    ## deparse(.) returning *one* string
    deparse1 <- function(expr, collapse = " ", width.cutoff = 500L, ...)
        paste(deparse(expr, width.cutoff, ...), collapse=collapse)
}
