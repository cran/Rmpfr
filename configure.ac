# autoconf this file --> configure
AC_INIT([Rmpfr])
AC_CONFIG_SRCDIR([DESCRIPTION])

dnl Select an optional include path, from a configure option
dnl or from an environment variable  >> mpfr_INCLUDE <<
AC_ARG_WITH([mpfr-include],
            AC_HELP_STRING([--with-mpfr-include=INCLUDE_PATH],
                           [the location of MPFR header files]),
            [mpfr_include_path=$withval])

if test -n "$mpfr_include_path" ; then
   MPFR_CPPFLAGS="-I${mpfr_include_path}"
elif test -n "${mpfr_INCLUDE}" ; then
   MPFR_CPPFLAGS="-I${mpfr_INCLUDE}"
fi

dnl Ditto for a library path  environment variable  >> mpfr_LDFLAGS <<
AC_ARG_WITH([mpfr-lib],
            AC_HELP_STRING([--with-mpfr-lib=LIB_PATH],
                           [the location of MPFR libraries]),
            [mpfr_lib_path=$withval])

if test -n "$mpfr_lib_path" ; then
   MPFR_LDFLAGS="-L$mpfr_lib_path ${LDFLAGS}"
elif test -n "${mpfr_LDFLAGS}" ; then
   MPFR_LDFLAGS="-L${mpfr_LDFLAGS} ${LDFLAGS}"
fi

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

dnl Setup the compilers
AC_PROG_CC
AC_PROG_CPP

dnl use the include and lib
CPPFLAGS="${CPPFLAGS} ${MPFR_CPPFLAGS}"
LDFLAGS="${LDFLAGS} ${MPFR_LDFLAGS}"

PKG_CPPFLAGS=$CPPFLAGS
PKG_CFLAGS=$CFLAGS
PKG_LDFLAGS=$LDFLAGS

AC_CHECK_HEADER([mpfr.h], ,
    [AC_MSG_ERROR([Header file mpfr.h not found; maybe use --with-mpfr-include=INCLUDE_PATH])])

AC_CHECK_HEADER([gmp.h], ,
    [AC_MSG_ERROR([Header file gmp.h not found; maybe use --with-mpfr-include=INCLUDE_PATH])])

AC_CHECK_LIB(gmp, __gmpz_init, ,
    [AC_MSG_ERROR([GNU MP not found, see README])])

AC_CHECK_LIB(mpfr, mpfr_init, ,
    [AC_MSG_ERROR([MPFR Library not found, see README])])

dnl Our Version check <= 3 -- etc seems to fail nowadays
dnl	src/Rmpfr_utils.h
dnl	~~~~~~~~~~~~~~~~~
AC_CHECK_LIB(mpfr, mpfr_digamma, ,
    [AC_MSG_ERROR([MPFR Library must be at least version 3.0.0, see README])])

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LDFLAGS)

AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT
