/*
 * MPFR - Multiple Precision Floating-Point Reliable Library
 * ----   -        -         -              -
 */
#include <Rmath.h>
/* imax2() */

#include "Rmpfr_utils.h"
#include "Syms.h"

/*------------------------------------------------------------------------*/

SEXP R_mpfr_get_version(void) {
    return mkString(mpfr_get_version());
}

SEXP R_mpfr_set_debug(SEXP I) {
    /* Set or get the C-global debugging level : */
    if(LENGTH(I) < 1 || INTEGER(I)[0] == NA_INTEGER)
	return ScalarInteger(R_mpfr_debug_);
    /* else : */
    R_mpfr_debug_ = asInteger(I);
    return I;
}

SEXP R_mpfr_get_default_prec(void) {
    return ScalarInteger((int) mpfr_get_default_prec());
}

SEXP R_mpfr_set_default_prec(SEXP prec) {
    SEXP ans = ScalarInteger((int) mpfr_get_default_prec());
    mpfr_set_default_prec((mp_prec_t) asInteger(prec));
    return ans;
}


#define INIT_1_SETUP(_X_, _R_)			\
    mpfr_t _R_;					\
						\
    mpfr_init2(_R_, R_mpfr_prec(_X_));		\
    R_asMPFR(_X_, _R_)

#define FINISH_1_RETURN(_R_, val)		\
    val = PROTECT(MPFR_as_R(_R_));		\
    mpfr_clear (_R_);				\
    mpfr_free_cache();				\
    UNPROTECT(1);				\
    return val

SEXP const_asMpfr(SEXP I, SEXP prec)
{
    SEXP val;
    mpfr_t r;
    mpfr_init2(r, asInteger(prec));

    switch(asInteger(I)) {
    case 1: mpfr_const_pi     (r, GMP_RNDN); break;
    case 2: mpfr_const_euler  (r, GMP_RNDN); break;
    case 3: mpfr_const_catalan(r, GMP_RNDN); break;
    default:
	error("invalid integer code {const_asMpfr()}"); /* -Wall */
    }

    FINISH_1_RETURN(r, val);
}

#define R_MPFR_Logic_Function(_FNAME, _MPFR_NAME)			\
SEXP _FNAME(SEXP x) {							\
    SEXP D = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));/* R list() */	\
    int n = length(D), i;						\
    SEXP val = PROTECT(allocVector(LGLSXP, n));				\
    int *r = LOGICAL(val);						\
    mpfr_t r_i;								\
    mpfr_init(r_i);							\
									\
    for(i=0; i < n; i++) {						\
	R_asMPFR(VECTOR_ELT(D, i), r_i);				\
	r[i] = _MPFR_NAME (r_i);					\
    }									\
									\
    mpfr_clear (r_i);							\
    mpfr_free_cache();							\
    UNPROTECT(2);							\
    return val;								\
}

R_MPFR_Logic_Function(R_mpfr_is_finite,   mpfr_number_p)
R_MPFR_Logic_Function(R_mpfr_is_infinite, mpfr_inf_p)
R_MPFR_Logic_Function(R_mpfr_is_integer,  mpfr_integer_p)
R_MPFR_Logic_Function(R_mpfr_is_na,       mpfr_nan_p)
R_MPFR_Logic_Function(R_mpfr_is_zero,     mpfr_zero_p)

#ifdef __NOT_ANY_MORE__

#define R_MPFR_1_Numeric_Function(_FNAME, _MPFR_NAME)			\
SEXP _FNAME(SEXP x) {							\
    SEXP D = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));/* an R list() */	\
    int n = length(D), i;						\
    SEXP val = PROTECT(allocVector(VECSXP, n));				\
    mpfr_t r_i;								\
    mpfr_init(r_i);							\
									\
    for(i=0; i < n; i++) {						\
	R_asMPFR(VECTOR_ELT(D, i), r_i);				\
	_MPFR_NAME(r_i, r_i, GMP_RNDN);					\
	SET_VECTOR_ELT(val, i, MPFR_as_R(r_i));				\
    }									\
									\
    mpfr_clear (r_i);							\
    mpfr_free_cache();							\
    UNPROTECT(2);							\
    return val;								\
}

R_MPFR_1_Numeric_Function(R_mpfr_erf, mpfr_erf)
R_MPFR_1_Numeric_Function(R_mpfr_erfc, mpfr_erfc)
R_MPFR_1_Numeric_Function(R_mpfr_zeta, mpfr_zeta)

R_MPFR_1_Numeric_Function(R_mpfr_eint, mpfr_eint)
R_MPFR_1_Numeric_Function(R_mpfr_j0, mpfr_j0)
R_MPFR_1_Numeric_Function(R_mpfr_j1, mpfr_j1)
R_MPFR_1_Numeric_Function(R_mpfr_y0, mpfr_y0)
R_MPFR_1_Numeric_Function(R_mpfr_y1, mpfr_y1)

#endif
/* __NOT_ANY_MORE__ */


#define R_MPFR_2_Numeric_Function(_FNAME, _MPFR_NAME)		\
SEXP _FNAME(SEXP x, SEXP y) {					\
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));		\
    SEXP yD = PROTECT(GET_SLOT(y, Rmpfr_Data_Sym));		\
    int nx = length(xD), ny = length(yD), i,			\
	n = (nx * ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;	\
    SEXP val = PROTECT(allocVector(VECSXP, n));			\
    mpfr_t x_i, y_i;						\
    mpfr_init(x_i); /* with default precision */		\
    mpfr_init(y_i); /* with default precision */		\
								\
    if (nx == ny || nx == 1 || ny == 1) mismatch = 0;		\
    else if (nx > 0 && ny > 0) {				\
	if (nx > ny) mismatch = nx % ny;			\
	else mismatch = ny % nx;				\
    }								\
								\
    for(i=0; i < n; i++) {					\
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);			\
	R_asMPFR(VECTOR_ELT(yD, i % ny), y_i);			\
	_MPFR_NAME(x_i, x_i, y_i, GMP_RNDN);			\
	SET_VECTOR_ELT(val, i, MPFR_as_R(x_i));			\
    }								\
								\
    MISMATCH_WARN;						\
    mpfr_clear (x_i); mpfr_clear (y_i);				\
    mpfr_free_cache();						\
    UNPROTECT(3);						\
    return val;							\
}

R_MPFR_2_Numeric_Function(R_mpfr_atan2, mpfr_atan2)
R_MPFR_2_Numeric_Function(R_mpfr_hypot, mpfr_hypot)


#define R_MPFR_2_Num_Long_Function(_FNAME, _MPFR_NAME)			\
SEXP _FNAME(SEXP x, SEXP y) {						\
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));			\
    int *yy = INTEGER(y);						\
    int nx = length(xD), ny = length(y), i,				\
	n = (nx * ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;		\
    SEXP val = PROTECT(allocVector(VECSXP, n));				\
    mpfr_t x_i;								\
									\
    if(TYPEOF(y) != INTSXP)						\
	error("MPFR_2_Num_Long...(d,mpfr): 'y' is not \"integer\"");	\
    mpfr_init(x_i); /* with default precision */			\
									\
    if (nx == ny || nx == 1 || ny == 1) mismatch = 0;			\
    else if (nx > 0 && ny > 0) {					\
	if (nx > ny) mismatch = nx % ny;				\
	else mismatch = ny % nx;					\
    }									\
									\
    for(i=0; i < n; i++) {						\
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);				\
	_MPFR_NAME(x_i, (long) yy[i % ny], x_i, GMP_RNDN);		\
	SET_VECTOR_ELT(val, i, MPFR_as_R(x_i));				\
    }									\
									\
    MISMATCH_WARN;							\
    mpfr_clear (x_i);							\
    mpfr_free_cache();							\
    UNPROTECT(2);							\
    return val;								\
}

R_MPFR_2_Num_Long_Function(R_mpfr_jn, mpfr_jn)
R_MPFR_2_Num_Long_Function(R_mpfr_yn, mpfr_yn)

