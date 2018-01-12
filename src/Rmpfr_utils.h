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
#include <R_ext/Print.h>

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

#if (MPFR_VERSION < MPFR_VERSION_NUM(3,0,0))
/* define back-compatibility types:*/
# define MPFR_RNDD GMP_RNDD
# define MPFR_RNDN GMP_RNDN
# define MPFR_RNDU GMP_RNDU
# define MPFR_RNDZ GMP_RNDZ
// # define MPFR_RNDA GMP_RNDA

# define mpfr_exp_t mp_exp_t

#endif

/*----------------------------------------*/

#ifdef _in_Rmpfr_init_
/* global */ int R_mpfr_debug_ = 0;
#else
extern       int R_mpfr_debug_;
#endif

/* A version of Rprintf() .. but only printing when .. is 'TRUE' :*/
static R_INLINE void R_mpfr_dbg_printf(int dbg_level, const char *format, ...)
{
    va_list(ap);
    if(R_mpfr_debug_ && R_mpfr_debug_ >= dbg_level) {
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

#define N_LIMBS(_PREC_) ceil(((double)_PREC_)/mp_bits_per_limb)

static R_INLINE int R_mpfr_nr_limbs(mpfr_t r)
{
    return N_LIMBS(mpfr_get_prec(r));
}

// Note: "in theory" we could set precBits > INT_MAX, but currently not in Rmpfr:
static R_INLINE void R_mpfr_check_prec(int prec)
{
    if(prec == NA_INTEGER)
	error("Precision(bit) is NA (probably from coercion)");
    if(prec < MPFR_PREC_MIN)
	error("Precision(bit) = %d < %ld (= MPFR_PREC_MIN)", prec, (long)MPFR_PREC_MIN);
/* 2018-01-01 gives a WARNING with clang:
 Found the following significant warnings:
   ./Rmpfr_utils.h:89:13: warning: comparison of constant 9223372036854775807 with expression of type 'int' is always false [-Wtautological-constant-out-of-range-compare]

... of course, I don't want a  WARN  in the CRAN checks, hence disable (for now):
    if(prec > MPFR_PREC_MAX)
	error("Precision(bit) = %d > %ld (= MPFR_PREC_MAX)", prec, (long)MPFR_PREC_MAX);
*/
    return;
}

#define R_mpfr_prec(x) INTEGER(GET_SLOT(x, Rmpfr_precSym))[0]


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
mpfr_rnd_t R_rnd2MP(SEXP rnd_mode);
SEXP d2mpfr1 (SEXP x, SEXP prec, SEXP rnd_mode);
SEXP d2mpfr1_(double x, int i_prec, mpfr_rnd_t rnd);
SEXP d2mpfr1_list(SEXP x, SEXP prec, SEXP rnd_mode);
SEXP mpfr2d(SEXP x, SEXP rnd_mode);
SEXP mpfr2i(SEXP x, SEXP rnd_mode);
SEXP mpfr2str(SEXP x, SEXP digits, SEXP maybe_full, SEXP base);
SEXP str2mpfr1_list(SEXP x, SEXP prec, SEXP base, SEXP rnd_mode);

#ifdef R_had_R_Outputfile_in_API
# ifndef WIN32
SEXP print_mpfr (SEXP x, SEXP digits);
SEXP print_mpfr1(SEXP x, SEXP digits);
# endif
#endif

SEXP Rmpfr_minus(SEXP x);
SEXP Rmpfr_abs(SEXP x);
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
SEXP R_mpfr_sumprod(SEXP x, SEXP y, SEXP minPrec, SEXP alternating);

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
SEXP R_mpfr_get_erange(SEXP kind);
SEXP R_mpfr_set_erange(SEXP kind, SEXP val);
SEXP R_mpfr_prec_range(SEXP ind);
SEXP R_mpfr_get_version(void);
SEXP R_mpfr_get_GMP_numb_bits(void);

SEXP const_asMpfr(SEXP I, SEXP prec, SEXP rnd_mode);

SEXP R_mpfr_is_finite(SEXP x);	SEXP R_mpfr_is_finite_A(SEXP x);
SEXP R_mpfr_is_infinite(SEXP x);SEXP R_mpfr_is_infinite_A(SEXP x);
SEXP R_mpfr_is_integer(SEXP x);	SEXP R_mpfr_is_integer_A(SEXP x);
SEXP R_mpfr_is_na(SEXP x);	SEXP R_mpfr_is_na_A(SEXP x);
SEXP R_mpfr_is_zero(SEXP x);    SEXP R_mpfr_is_zero_A(SEXP x);

SEXP R_mpfr_atan2(SEXP x, SEXP y, SEXP rnd_mode);
SEXP R_mpfr_hypot(SEXP x, SEXP y, SEXP rnd_mode);
SEXP R_mpfr_beta (SEXP x, SEXP y, SEXP rnd_mode);
SEXP R_mpfr_lbeta(SEXP x, SEXP y, SEXP rnd_mode);

SEXP R_mpfr_jn(SEXP x, SEXP n, SEXP rnd_mode);
SEXP R_mpfr_yn(SEXP x, SEXP n, SEXP rnd_mode);
SEXP R_mpfr_fac   (SEXP n, SEXP prec, SEXP rnd_mode);
SEXP R_mpfr_choose(SEXP a, SEXP n, SEXP rnd_mode);
SEXP R_mpfr_poch  (SEXP a, SEXP n, SEXP rnd_mode);
SEXP R_mpfr_round (SEXP x, SEXP prec, SEXP rnd_mode);

/* #ifdef __cplusplus */
/* } */
/* #endif */

#endif /* R_MPFR_MUTILS_H_ */
