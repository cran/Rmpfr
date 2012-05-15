/*
 * MPFR - Multiple Precision Floating-Point Reliable Library
 * ----   -        -         -              -
 */
#include <Rmath.h>
/* imax2() */

#include "Rmpfr_utils.h"
#include "Syms.h"

#ifdef DEBUG_Rmpfr
/* ONLY for debugging  !! */
# ifndef WIN32
#  include <Rinterface.h>
# endif
#define R_PRT(_X_) mpfr_out_str (R_Outputfile, 10, 0, _X_, GMP_RNDD)
#endif

int my_mpfr_beta (mpfr_t ROP, mpfr_t X, mpfr_t Y, mp_rnd_t RND);
int my_mpfr_lbeta(mpfr_t ROP, mpfr_t X, mpfr_t Y, mp_rnd_t RND);

int my_mpfr_choose(mpfr_t ROP, long n,    mpfr_t X, mp_rnd_t RND);
int my_mpfr_poch  (mpfr_t ROP, long n,    mpfr_t X, mp_rnd_t RND);
int my_mpfr_round (mpfr_t ROP, long prec, mpfr_t X, mp_rnd_t RND);
/* argument order above must match the one of mpfr_jn() etc .. */

/*------------------------------------------------------------------------*/
int my_mpfr_beta (mpfr_t R, mpfr_t X, mpfr_t Y, mp_rnd_t RND)
{
    /* NOTA BENE: When called, R is typically *identical* to X
     *            ==> can use  R  only at the very end! */
    int ans;
    mpfr_t s,t;
    mp_prec_t p_X = mpfr_get_prec(X), p_Y = mpfr_get_prec(Y);
    if(p_X < p_Y) p_X = p_Y;
    mpfr_init2(s, p_X);
    mpfr_init2(t, p_X);
    /* FIXME: check each 'ans' below, and return when not ok ... */
    ans = mpfr_gamma(s, X, RND);
    ans = mpfr_gamma(t, Y, RND);
    ans = mpfr_mul(s, s, t, RND); /* s = gamma(X) * gamma(Y) */

#ifdef DEBUG_Rmpfr
    Rprintf("my_mpfr_beta(): t = gamma(Y)= "); R_PRT(t);
    Rprintf("\n   s = G(X) * G(Y) = ");        R_PRT(s); Rprintf("\n");
    Rprintf("\n   X = "); R_PRT(X);
    Rprintf("\n   Y = "); R_PRT(Y);
#endif
    ans = mpfr_add(X, X, Y, RND);
    ans = mpfr_gamma(Y, X, RND);  /* Y. = gamma(X + Y) */
#ifdef DEBUG_Rmpfr
    Rprintf("\n  X. = X + Y =  "); R_PRT(X);
    Rprintf("\n  Y. = gamma(X.)= "); R_PRT(Y);
    Rprintf("\n");
#endif
    ans = mpfr_div(R, s, Y, RND);
    mpfr_clear (s);
    mpfr_clear (t);
    /* mpfr_free_cache() must be called in the caller !*/
    return ans;
}

int my_mpfr_lbeta(mpfr_t R, mpfr_t X, mpfr_t Y, mp_rnd_t RND)
{
    int ans;
    mpfr_t s,t;
    mp_prec_t p_X = mpfr_get_prec(X), p_Y = mpfr_get_prec(Y);
    if(p_X < p_Y) p_X = p_Y;
    mpfr_init2(s, p_X);
    mpfr_init2(t, p_X);
    /* FIXME: check each 'ans' below, and return when not ok ... */
    ans = mpfr_lngamma(s, X, RND);
    ans = mpfr_lngamma(t, Y, RND);
    ans = mpfr_add(s, s, t, RND); /* s = lgamma(X) + lgamma(Y) */
    ans = mpfr_add(X, X, Y, RND);
    ans = mpfr_lngamma(Y, X, RND);/* Y = lgamma(X + Y) */
    ans = mpfr_sub(R, s, Y, RND);
    mpfr_clear (s);
    mpfr_clear (t);
    /* mpfr_free_cache() must be called in the caller !*/
    return ans;
}

/** Binomial Coefficient --
 * all initialization and cleanup is called in the caller
 * @result R = choose(X, n)
 */
int my_mpfr_choose (mpfr_t R, long n, mpfr_t X, mp_rnd_t RND)
{
    int ans;
    long i;
    mpfr_t r, x;
    mp_prec_t p_X = mpfr_get_prec(X);

    mpfr_init2(x, p_X); mpfr_set(x, X, RND);
    mpfr_init2(r, p_X);
    if(n > 0) {
	mpfr_set(r, X, RND);
	for(i=1; i < n; ) {
	    mpfr_sub_si(x, x, 1L, RND); // x = X - i
	    mpfr_mul   (r, r, x, RND); // r := r * x = X(X-1)..(X-i)
	    mpfr_div_si(r, r, ++i, RND);
	    // r := r / (i+1) =  X(X-1)..(X-i) / (1*2..*(i+1))
#ifdef DEBUG_Rmpfr
	    Rprintf("my_mpfr_choose(): X (= X_0 - %d)= ", i); R_PRT(x);
	    Rprintf("\n --> r ="); R_PRT(r); Rprintf("\n");
#endif
	}
    }
    else // n = 0
	mpfr_set_si(r, (long) 1, RND);
    ans = mpfr_set(R, r, RND);
    mpfr_clear (x);
    mpfr_clear (r);
    return ans;
}

/** Pochhammer Symbol -- *rising* factorial   x * (x+1) * ... (x+n-1)
 * all initialization and cleanup is called in the caller
 */
int my_mpfr_poch (mpfr_t R, long n, mpfr_t X, mp_rnd_t RND)
{
    int ans;
    long i;
    mpfr_t r, x;
    mp_prec_t p_X = mpfr_get_prec(X);

    mpfr_init2(x, p_X); mpfr_set(x, X, RND);
    mpfr_init2(r, p_X);
    if(n > 0) {
	mpfr_set(r, X, RND);
	for(i=1; i < n; i++) {
	    mpfr_add_si(x, x, 1L, RND); // x = X + i
	    mpfr_mul(r, r, x, RND); // r := r * x = X(X+1)..(X+i)
#ifdef DEBUG_Rmpfr
	    Rprintf("my_mpfr_poch(): X (= X_0 + %d)= ", i); R_PRT(x);
	    Rprintf("\n --> r ="); R_PRT(r); Rprintf("\n");
#endif
	}
    }
    else // n = 0
	mpfr_set_si(r, (long) 1, RND);
    ans = mpfr_set(R, r, RND);
    mpfr_clear (x);
    mpfr_clear (r);
    return ans;
}

/** round to (binary) bits, not (decimal) digits
 */
int my_mpfr_round (mpfr_t R, long prec, mpfr_t X, mp_rnd_t RND)
{
    int ans;
    if(prec < MPFR_PREC_MIN)
	error("prec = %d < %d  is too small", prec, MPFR_PREC_MIN);
    if(prec > MPFR_PREC_MAX)
	error("prec = %d > %d  is too large", prec, MPFR_PREC_MAX);
    mpfr_set(R, X, RND);
    ans = mpfr_prec_round(R, (mp_prec_t) prec, RND);
    return ans;
}

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
    case 4: mpfr_const_log2   (r, GMP_RNDN); break;
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

SEXP R_mpfr_fac (SEXP n_, SEXP prec)
{
    int n = length(n_), i, *nn;
    SEXP n_t, val = PROTECT(allocVector(VECSXP, n)); int nprot = 1;
    mpfr_t r_i;
    if(TYPEOF(n_) != INTSXP) {
	PROTECT(n_t = coerceVector(n_, INTSXP)); nprot++;/* or bail out*/
	nn = INTEGER(n_t);
    } else {
	nn = INTEGER(n_);
    }
    mpfr_init2(r_i, (mp_prec_t) asInteger(prec));
    for(i=0; i < n; i++) {
	// never happens when called from R:
	if(nn[i] < 0) error("R_mpfr_fac(%d): negative n.", nn[i]);
	mpfr_fac_ui(r_i, nn[i], GMP_RNDN);
	SET_VECTOR_ELT(val, i, MPFR_as_R(r_i));
    }

    mpfr_clear(r_i);
    mpfr_free_cache();
    UNPROTECT(nprot);
    return val;
}

#ifdef __NOT_ANY_MORE__
//       ------------ as we deal with these "as if Math() group"
// via Math_mpfr() in ./Ops.c

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

R_MPFR_1_Numeric_Function(R_mpfr_ai, mpfr_ai)

#endif
/* __NOT_ANY_MORE__ */


#define R_MPFR_2_Numeric_Function(_FNAME, _MPFR_NAME)	\
SEXP _FNAME(SEXP x, SEXP y) {				\
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));	\
    SEXP yD = PROTECT(GET_SLOT(y, Rmpfr_Data_Sym));	\
    int nx = length(xD), ny = length(yD), i,		\
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny);	\
    SEXP val = PROTECT(allocVector(VECSXP, n));		\
    mpfr_t x_i, y_i;					\
    mpfr_init(x_i); /* with default precision */	\
    mpfr_init(y_i); /* with default precision */	\
							\
    for(i=0; i < n; i++) {				\
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);		\
	R_asMPFR(VECTOR_ELT(yD, i % ny), y_i);		\
	_MPFR_NAME(x_i, x_i, y_i, GMP_RNDN);		\
	SET_VECTOR_ELT(val, i, MPFR_as_R(x_i));		\
    }							\
							\
    mpfr_clear (x_i); mpfr_clear (y_i);			\
    mpfr_free_cache();					\
    UNPROTECT(3);					\
    return val;						\
}

R_MPFR_2_Numeric_Function(R_mpfr_atan2, mpfr_atan2)
R_MPFR_2_Numeric_Function(R_mpfr_hypot, mpfr_hypot)

R_MPFR_2_Numeric_Function(R_mpfr_beta,  my_mpfr_beta)
R_MPFR_2_Numeric_Function(R_mpfr_lbeta, my_mpfr_lbeta)


#define R_MPFR_2_Num_Long_Function(_FNAME, _MPFR_NAME)			\
SEXP _FNAME(SEXP x, SEXP y) {						\
    SEXP xD, yt, val;							\
    int *yy, n, nx, ny = length(y), i, nprot = 0;			\
    mpfr_t x_i;								\
									\
    if(TYPEOF(y) != INTSXP) {						\
	PROTECT(yt = coerceVector(y, INTSXP)); nprot++;/* or bail out*/ \
	yy = INTEGER(yt);						\
    } else {								\
	yy = INTEGER(y);						\
    }									\
    PROTECT(xD = GET_SLOT(x, Rmpfr_Data_Sym));	nprot++;		\
    nx = length(xD);							\
    n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny);			\
    PROTECT(val = allocVector(VECSXP, n)); 	nprot++;		\
    mpfr_init(x_i); /* with default precision */			\
									\
    for(i=0; i < n; i++) {						\
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);				\
	_MPFR_NAME(x_i, (long) yy[i % ny], x_i, GMP_RNDN);		\
	SET_VECTOR_ELT(val, i, MPFR_as_R(x_i));				\
    }									\
									\
    mpfr_clear (x_i);							\
    mpfr_free_cache();							\
    UNPROTECT(nprot);							\
    return val;								\
}

R_MPFR_2_Num_Long_Function(R_mpfr_jn, mpfr_jn)
R_MPFR_2_Num_Long_Function(R_mpfr_yn, mpfr_yn)
R_MPFR_2_Num_Long_Function(R_mpfr_choose, my_mpfr_choose)
R_MPFR_2_Num_Long_Function(R_mpfr_poch, my_mpfr_poch)
R_MPFR_2_Num_Long_Function(R_mpfr_round, my_mpfr_round)

