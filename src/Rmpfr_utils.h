#ifndef R_MPFR_MUTILS_H
#define R_MPFR_MUTILS_H

/* #ifdef __cplusplus */
/* extern "C" { */
/* #endif */

#include <ctype.h>

#include <stdarg.h>
/* for va_list ..*/

#include <R.h>  /* includes Rconfig.h */
#include <Rversion.h>
/* for NEW_OBJECT(), GET_SLOT + Rinternals.h : */
#include <Rdefines.h>

/* must come *after* the above, e.g., for
   mpfr_out_str()  (which needs stdio): */
#include <gmp.h>
#include <mpfr.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Rmpfr", String)
#else
#define _(String) (String)
#endif

/*----------------------------------------*/

#ifdef _in_Rmpfr_init_
/* global */ int R_mpfr_debug_ = 0;
#else
extern       int R_mpfr_debug_;
#endif

/* A version of Rprintf() .. but only printing when .. is 'TRUE' :*/
static R_INLINE void R_mpfr_dbg_printf(const char *format, ...)
{
    va_list(ap);
    if(R_mpfr_debug_) {
	Rprintf("mpfr.debug[%d]: ", R_mpfr_debug_);
	va_start(ap, format);
	REvprintf(format, ap);
	va_end(ap);
    }
}



/* This is from Matrix/src/Mutils.h : */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
    SEXP val = allocVector(type, length);

    SET_SLOT(obj, nm, val);
    return val;
}

static R_INLINE int R_mpfr_nr_limbs(mpfr_t r)
{
    int d = (int)mpfr_get_prec(r),
	nr = d/mp_bits_per_limb;
    if (d % mp_bits_per_limb) nr++;
    return nr;
}

#define R_mpfr_prec(x) INTEGER(GET_SLOT(x, Rmpfr_precSym))[0]

#define N_LIMBS(_PREC_) ceil(((double)_PREC_)/mp_bits_per_limb)


#define MISMATCH_WARN							\
    if (mismatch)							\
	warning(_("longer object length is not a multiple of shorter object length"))

#define SET_MISMATCH					\
    if (nx == ny || nx == 1 || ny == 1) mismatch = 0;	\
    else if (nx > 0 && ny > 0) {			\
	if (nx > ny) mismatch = nx % ny;		\
	else mismatch = ny % nx;			\
    }



/* ./convert.c : */
SEXP d2mpfr  (SEXP x, SEXP prec);
SEXP d2mpfr1 (SEXP x, SEXP prec);
SEXP d2mpfr1_(double x, int i_prec);
SEXP d2mpfr1_list(SEXP x, SEXP prec);
SEXP mpfr2d(SEXP x);
SEXP mpfr2i(SEXP x);
SEXP mpfr2str(SEXP x, SEXP digits);
SEXP str2mpfr1_list(SEXP x, SEXP prec, SEXP base);

SEXP print_mpfr (SEXP x, SEXP digits);
SEXP print_mpfr1(SEXP x, SEXP digits);

SEXP Math_mpfr(SEXP x, SEXP op);
SEXP Arith_mpfr(SEXP x, SEXP y, SEXP op);
SEXP Arith_mpfr_i(SEXP x, SEXP y, SEXP op);
SEXP Arith_i_mpfr(SEXP x, SEXP y, SEXP op);
SEXP Arith_mpfr_d(SEXP x, SEXP y, SEXP op);
SEXP Arith_d_mpfr(SEXP x, SEXP y, SEXP op);

SEXP Compare_mpfr(SEXP x, SEXP y, SEXP op);
SEXP Compare_mpfr_i(SEXP x, SEXP y, SEXP op);
SEXP Compare_mpfr_d(SEXP x, SEXP y, SEXP op);

SEXP Summary_mpfr(SEXP x, SEXP na_rm, SEXP op);

#ifdef __NOT_ANY_MORE__
/* deprecated: */
SEXP exp_mpfr1(SEXP x);
SEXP log_mpfr1(SEXP x);
#endif

void R_asMPFR(SEXP x, mpfr_ptr r);
SEXP MPFR_as_R(mpfr_t r);

/* ./utils.c */
SEXP R_mpfr_set_debug(SEXP I);
SEXP R_mpfr_set_default_prec(SEXP prec);
SEXP R_mpfr_get_default_prec(void);


SEXP const_asMpfr(SEXP I, SEXP prec);

SEXP R_mpfr_is_finite(SEXP x);
SEXP R_mpfr_is_infinite(SEXP x);
SEXP R_mpfr_is_integer(SEXP x);
SEXP R_mpfr_is_na(SEXP x);
SEXP R_mpfr_is_zero(SEXP x);

SEXP R_mpfr_jn(SEXP x, SEXP n);
SEXP R_mpfr_yn(SEXP x, SEXP n);
SEXP R_mpfr_atan2(SEXP x, SEXP y);
SEXP R_mpfr_hypot(SEXP x, SEXP y);


/* #ifdef __cplusplus */
/* } */
/* #endif */

#endif /* R_MPFR_MUTILS_H_ */
