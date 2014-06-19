/*
 * MPFR - Multiple Precision Floating-Point Reliable Library
 * ----   -        -         -              -
 *
 * Arithmetic, Math, etc
 */

#include <limits.h>
#include <Rmath.h>
/* imax2() */

#include "Rmpfr_utils.h"
#include "Syms.h"


SEXP Rmpfr_minus(SEXP x)
{
    int n = length(x);
    SEXP val = PROTECT(duplicate(x));
    for(int i=0; i < n; i++) {
	int sign = asInteger(GET_SLOT(VECTOR_ELT(x,i), Rmpfr_signSym));
	SEXP r_i = VECTOR_ELT(val, i);
	SET_SLOT(r_i, Rmpfr_signSym, ScalarInteger(-sign));
	SET_VECTOR_ELT(val, i, r_i);
    }

    UNPROTECT(1);
    return val;
} /* Rmpfr_minus() */

SEXP Rmpfr_abs(SEXP x)
{
    int n = length(x);
    SEXP val = PROTECT(duplicate(x));
    for(int i=0; i < n; i++) {
	SEXP r_i = VECTOR_ELT(val, i);
	SET_SLOT(r_i, Rmpfr_signSym, ScalarInteger(1));
	SET_VECTOR_ELT(val, i, r_i);
    }
    UNPROTECT(1);
    return val;
} /* Rmpfr_abs() */


/*------------------------------------------------------------------------*/

SEXP Math_mpfr(SEXP x, SEXP op)
{
#ifdef using_Data_slot
    SEXP D = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));
#else
# define D x
#endif
    mpfr_prec_t current_prec = mpfr_get_default_prec();
    int n = length(D), i_op = asInteger(op), i;

    SEXP val = PROTECT(allocVector(VECSXP, n));
    mpfr_t R_i, cum;
    Rboolean is_cum = (71 <= i_op && i_op <= 74);

    mpfr_init(R_i); /* with default precision */
    if(is_cum) { // cummax, cumsum, ...
	mpfr_init(cum);
	switch(i_op) {
	case 71: /* cummax */  mpfr_set_inf(cum, -1);/* := -Inf */; break;
	case 72: /* cummin */  mpfr_set_inf(cum, +1);/* := +Inf */; break;
	case 73: /* cumprod */ mpfr_set_d(cum, 1., MPFR_RNDZ);/* := 1 */; break;
	case 74: /* cumsum */  mpfr_set_d(cum, 0., MPFR_RNDZ);/* := 0 */; break;
	}
    }
    for(i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(D, i), R_i);
	if(is_cum) { /* hence using cum */
	    mpfr_prec_t i_prec = mpfr_get_prec(R_i);
	    if(current_prec < i_prec) /* increase precision */ {
		current_prec = i_prec;
		mpfr_prec_round(cum, i_prec, MPFR_RNDN);
	    }
	}


#define NOT_YET error("Math op. %d not yet implemented", i_op)
	switch(i_op) {
	    /* Note we assign use R_i as "input and output" ==> *same*
	       precision, even though in some cases the result may
	       need higher precision */

	case  0: /* trunc */ mpfr_trunc(R_i, R_i); break;
	case  1: /* floor */ mpfr_floor(R_i, R_i); break;
	case  2: /* ceiling*/ mpfr_ceil(R_i, R_i); break;
	case  3: /* sqrt */  mpfr_sqrt(R_i, R_i, MPFR_RNDN); break;
	case  4: /* sign */
	    error("'sign' is dealt with in R.  Should not happen, please report");
	    break;
	case 10: /* exp */   mpfr_exp(R_i, R_i, MPFR_RNDN); break;
	case 11: /* expm1 */ mpfr_expm1(R_i, R_i, MPFR_RNDN); break;
	case 12: /* log1p */ mpfr_log1p(R_i, R_i, MPFR_RNDN); break;
	case 13: /* log */   mpfr_log  (R_i, R_i, MPFR_RNDN); break;
	case 14: /* log2 */  mpfr_log2 (R_i, R_i, MPFR_RNDN); break;
	case 15: /* log10 */ mpfr_log10(R_i, R_i, MPFR_RNDN); break;

	case 20: /* cos */   mpfr_cos  (R_i, R_i, MPFR_RNDN); break;
	case 21: /* sin */   mpfr_sin  (R_i, R_i, MPFR_RNDN); break;
	case 22: /* tan */   mpfr_tan  (R_i, R_i, MPFR_RNDN); break;
	case 23: /* acos */  mpfr_acos (R_i, R_i, MPFR_RNDN); break;
	case 24: /* asin */  mpfr_asin (R_i, R_i, MPFR_RNDN); break;
	case 25: /* atan */  mpfr_atan (R_i, R_i, MPFR_RNDN); break;
	case 30: /* cosh */  mpfr_cosh (R_i, R_i, MPFR_RNDN); break;
	case 31: /* sinh */  mpfr_sinh (R_i, R_i, MPFR_RNDN); break;
	case 32: /* tanh */  mpfr_tanh (R_i, R_i, MPFR_RNDN); break;
	case 33: /* acosh */ mpfr_acosh(R_i, R_i, MPFR_RNDN); break;
	case 34: /* asinh */ mpfr_asinh(R_i, R_i, MPFR_RNDN); break;
	case 35: /* atanh */ mpfr_atanh(R_i, R_i, MPFR_RNDN); break;

	case 40: /* lgamma */ { int sgn[1];
		mpfr_lgamma(R_i, sgn, R_i, MPFR_RNDN); break; }
	case 41: /* gamma */ mpfr_gamma(R_i, R_i, MPFR_RNDN); break;
	case 42: /* digamma */
#if (MPFR_VERSION < MPFR_VERSION_NUM(3,0,0))
	    error("digamma() not implemented in oldish MPFR library version '%s'",
			MPFR_VERSION_STRING);
#else
	    mpfr_digamma(R_i, R_i, MPFR_RNDN); break;
#endif
	case 43: /* trigamma */ NOT_YET; break;

	case 47: /* cospi */ {
	    mpfr_prec_t i_prec = mpfr_get_prec(R_i);
	    mpfr_t tmp;
	    mpfr_init2(tmp, i_prec);
	    mpfr_abs(R_i, R_i, MPFR_RNDN); // R_i :=  | R_i |
	    mpfr_set_si(tmp, (long) 2, MPFR_RNDN); // tmp := 2
	    // R_i := R_i mod 2 :
	    mpfr_fmod(R_i, R_i, tmp, MPFR_RNDN);

	    if(mpfr_cmp_d(R_i, 0.5) == 0 || mpfr_cmp_d(R_i, 1.5) == 0)
		mpfr_set_zero(R_i, +1);
	    else if(mpfr_cmp_si(R_i, (long) 1) == 0)
		mpfr_set_si(R_i, (long) -1, MPFR_RNDN);
	    else if(mpfr_cmp_si(R_i, (long) 0) == 0)
		mpfr_set_si(R_i, (long)  1, MPFR_RNDN);
	    else { // otherwise return  cos(pi * x):
		mpfr_const_pi (tmp, MPFR_RNDN);
		mpfr_mul(R_i, R_i, tmp, MPFR_RNDN);
		mpfr_cos(R_i, R_i, MPFR_RNDN);
	    }
	    break;
	}

	case 48: /* sinpi */  {
	    mpfr_prec_t i_prec = mpfr_get_prec(R_i);
	    mpfr_t tmp;
	    mpfr_init2(tmp, i_prec);
	    mpfr_set_si(tmp, (long) 2, MPFR_RNDN); // tmp := 2
	    // R_i := R_i mod 2 :
	    mpfr_fmod(R_i, R_i, tmp, MPFR_RNDN);
	    // map (-2,2) --> (-1,1] :
	    if(mpfr_cmp_si(R_i, (long) -1) <= 0)
		mpfr_add(R_i, R_i, tmp, MPFR_RNDN);
	    else if(mpfr_cmp_si(R_i, (long) 1) > 0)
		mpfr_sub(R_i, R_i, tmp, MPFR_RNDN);

	    if(mpfr_integer_p(R_i)) // x = 0 or 1 : ==> sin(pi*x) = 0
		mpfr_set_zero(R_i, +1);
	    else if(mpfr_cmp_d(R_i, 0.5) == 0)
		mpfr_set_si(R_i, (long) 1, MPFR_RNDN);
	    else if(mpfr_cmp_d(R_i, -0.5) == 0)
		mpfr_set_si(R_i, (long) -1, MPFR_RNDN);
	    else { // otherwise return  sin(pi * x):
		mpfr_const_pi (tmp, MPFR_RNDN);
		mpfr_mul(R_i, R_i, tmp, MPFR_RNDN);
		mpfr_sin(R_i, R_i, MPFR_RNDN);
	    }
	    break;
	}

	case 49: /* tanpi */  {
	    mpfr_prec_t i_prec = mpfr_get_prec(R_i);
	    mpfr_t tmp;
	    mpfr_init2(tmp, i_prec);
	    mpfr_set_si(tmp, (long) 1, MPFR_RNDN); // tmp := 1
	    // R_i := R_i mod 1 :
	    mpfr_fmod(R_i, R_i, tmp, MPFR_RNDN);
	    // map (-1,1) --> (-1/2, 1/2] :
	    if(mpfr_cmp_d(R_i, (double) -0.5) <= 0)
		mpfr_add(R_i, R_i, tmp, MPFR_RNDN);
	    else if(mpfr_cmp_d(R_i, (double) 0.5) > 0)
		mpfr_sub(R_i, R_i, tmp, MPFR_RNDN);

	    if(mpfr_zero_p(R_i)) // x = 0 : ==> tan(pi*x) = 0
		mpfr_set_zero(R_i, +1);
	    else if(mpfr_cmp_d(R_i, 0.5) == 0)
		mpfr_set_si(R_i, (long) 1, MPFR_RNDN);
	    else if(mpfr_cmp_d(R_i, -0.5) == 0)
		mpfr_set_si(R_i, (long) -1, MPFR_RNDN);
	    else {
		// otherwise return  tan(pi * x):
		mpfr_const_pi (tmp, MPFR_RNDN);
		mpfr_mul(R_i, R_i, tmp, MPFR_RNDN);
		mpfr_tan(R_i, R_i, MPFR_RNDN);
	    }
	    break;
	}


	case 71: /* cummax */ mpfr_max(cum, cum, R_i, MPFR_RNDN); break;
	case 72: /* cummin */ mpfr_min(cum, cum, R_i, MPFR_RNDN); break;
	case 73: /* cumprod*/ mpfr_mul(cum, cum, R_i, MPFR_RNDN); break;
	case 74: /* cumsum */ mpfr_add(cum, cum, R_i, MPFR_RNDN); break;

/*---   more functions from the  mpfr - library but not in R "Math" : ---*/
	case 101: mpfr_erf (R_i, R_i, MPFR_RNDN); break;
	case 102: mpfr_erfc(R_i, R_i, MPFR_RNDN); break;

	case 104: mpfr_zeta(R_i, R_i, MPFR_RNDN); break;

	case 106: mpfr_eint(R_i, R_i, MPFR_RNDN); break;
	case 107:
#if (MPFR_VERSION < MPFR_VERSION_NUM(2,4,0))
	    error("Li2() not implemented in oldish MPFR library version '%s'",
		  MPFR_VERSION_STRING);
#else
	    mpfr_li2 (R_i, R_i, MPFR_RNDN); break;
#endif
	case 111: mpfr_j0(R_i, R_i, MPFR_RNDN); break;
	case 112: mpfr_j1(R_i, R_i, MPFR_RNDN); break;
	case 113: mpfr_y0(R_i, R_i, MPFR_RNDN); break;
	case 114: mpfr_y1(R_i, R_i, MPFR_RNDN); break;
	case 120:
#if (MPFR_VERSION < MPFR_VERSION_NUM(3,0,0))
	    error("Ai() not implemented in oldish MPFR library version '%s'",
			MPFR_VERSION_STRING);
#else
	    mpfr_ai(R_i, R_i, MPFR_RNDN); break;
#endif

	default:
	    error("invalid op code (%d) in Math_mpfr", i_op);
	} // end{switch()}

	if(is_cum)
	    SET_VECTOR_ELT(val, i, MPFR_as_R(cum));
	else
	    SET_VECTOR_ELT(val, i, MPFR_as_R(R_i));
    }

    mpfr_clear (R_i);
    if(is_cum) mpfr_clear(cum);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(2);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Math_mpfr() */
#undef NOT_YET


//  %%   operator -- do what R does: ~/R/D/r-devel/R/src/main/arithmetic.c
// -----    --> it uses  %% (only sometimes!)  and myfmod();
//  .... ok, now checked in ../tests/arith-ex.R
//                          ~~~~~~~~~~~~~~~~~~~

// CARE: caller can use  R_mpfr_mod(x, x, y) -- i.e., r == x as pointers!
static int R_mpfr_mod(mpfr_t r, mpfr_t x, mpfr_t y, mpfr_rnd_t RND)
{
    if(mpfr_nan_p(y) || mpfr_nan_p(x)) {
	mpfr_set_nan(r); return 0;
    }
    int s_y = mpfr_sgn(y);// --> {-1, 0, 1}
    if(s_y == 0) { // y = 0  |-> NaN :
	mpfr_set_nan(r); return 0;
    }

    mpfr_t rr; mpfr_init_set (rr, r, RND);// a copy

    int s = mpfr_fmod(r, x, y, RND);// CARE: if(r is x) will thrash x
    if((s_y > 0 && mpfr_sgn(r) < 0) || // as R :  (-5) %%   3   |-->  1
       (s_y < 0 && mpfr_sgn(r) > 0))   // as R :    5  %% (-3)  |--> -1
	s += mpfr_add(r, r, y, RND);

    return s;
}


SEXP Arith_mpfr(SEXP x, SEXP y, SEXP op)
{
#ifdef using_Data_slot
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym)),
	 yD = PROTECT(GET_SLOT(y, Rmpfr_Data_Sym));
#else
# define xD x
# define yD y
#endif
    int nx = length(xD), ny = length(yD), i_op = asInteger(op), i,
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;

    SEXP val = PROTECT(allocVector(VECSXP, n));
    mpfr_t x_i, y_i;

    mpfr_init(x_i); /* with default precision */
    mpfr_init(y_i);

    SET_MISMATCH;
    for(i=0; i < n; i++) {
	mpfr_prec_t x_prec, y_prec;

	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i); x_prec = mpfr_get_prec(x_i);
	R_asMPFR(VECTOR_ELT(yD, i % ny), y_i); y_prec = mpfr_get_prec(y_i);

	if(x_prec < y_prec) {/* increase it, since it will store the result */
	    mpfr_prec_round (x_i, y_prec, MPFR_RNDN);
	    x_prec = y_prec;
	}
	switch(i_op) {
	    /* Note we assign use x_i as "input and output" ==> *same*
	       precision, even though in some cases the result may
	       need higher precision */
	case  1: /*  +  */  mpfr_add (x_i, x_i, y_i, MPFR_RNDN); break;
	case  2: /*  -  */  mpfr_sub (x_i, x_i, y_i, MPFR_RNDN); break;
	case  3: /*  *  */  mpfr_mul (x_i, x_i, y_i, MPFR_RNDN); break;
	case  4: /*  ^  */  mpfr_pow (x_i, x_i, y_i, MPFR_RNDN); break;
	case  5: /* %%  */ R_mpfr_mod(x_i, x_i, y_i, MPFR_RNDN); break;
	case  6: /* %/% */ {
	    mpfr_t r;
	    mpfr_init(r);
	    if(mpfr_get_prec(r) < x_prec)
		mpfr_set_prec (r, x_prec);

	    // want to ensure   x  ==  (x %% y) +  y * ( x %/% y )
	    //        <==>      x - (x %% y)  ==   y * ( x %/% y )
	    //        <==>    [ x - (x %% y) ] / y  == ( x %/% y )
	    R_mpfr_mod(r, x_i, y_i, MPFR_RNDN);// r := x %% y,
	    mpfr_sub (x_i, x_i, r, MPFR_RNDN); // x~ = x - r =  x - (x %% y)
	    mpfr_div (x_i, x_i,y_i,MPFR_RNDN); // x = x~ / y = (x - (x %% y))/y
	    mpfr_clear(r); break;
	}
	case  7: /*  /  */ mpfr_div(x_i, x_i, y_i, MPFR_RNDN); break;

	default:
	    error("invalid op code (%d) in Arith_mpfr", i_op);
	}

	SET_VECTOR_ELT(val, i, MPFR_as_R(x_i));
    }
    MISMATCH_WARN;

    mpfr_clear (x_i); mpfr_clear (y_i);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(3);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Arith_mpfr */


SEXP Arith_mpfr_i(SEXP x, SEXP y, SEXP op)
{
#ifdef using_Data_slot
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));
#else
# define xD x
#endif
    int *yy = INTEGER(y);
    int nx = length(xD), ny = length(y), i_op = asInteger(op), i,
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;

    if(TYPEOF(y) != INTSXP)
	error("Arith[%d](mpfr,i): 'y' is not a \"integer\"", i_op);

    SEXP val = PROTECT(allocVector(VECSXP, n));
    mpfr_t x_i;
    mpfr_init(x_i); /* with default precision */

    SET_MISMATCH;
    for(i=0; i < n; i++) {
	int i_ = i % ny;
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);
	switch(i_op) {
	    /* Note we assign use x_i as "input and output" ==> *same*
	       precision, even though in some cases the result may
	       need higher precision */
	case  1: /*  +  */ mpfr_add_si(x_i, x_i, (long) yy[i_], MPFR_RNDN); break;
	case  2: /*  -  */ mpfr_sub_si(x_i, x_i, (long) yy[i_], MPFR_RNDN); break;
	case  3: /*  *  */ mpfr_mul_si(x_i, x_i, (long) yy[i_], MPFR_RNDN); break;
	case  4: /*  ^  */ mpfr_pow_si(x_i, x_i, (long) yy[i_], MPFR_RNDN); break;
	case  5: /* %%  */ {
	    mpfr_t yy_i;
	    mpfr_init_set_si(yy_i, (long) yy[i_], MPFR_RNDN);
	    R_mpfr_mod(x_i, x_i, yy_i, MPFR_RNDN);
	    mpfr_clear(yy_i); break;
	}
	case  6: /* %/% */ {
	    mpfr_t r, yy_i;
	    mpfr_init(r);
	    mpfr_prec_t x_prec = mpfr_get_prec(x_i);
	    if(mpfr_get_prec(r) < x_prec)
		mpfr_set_prec (r, x_prec);

	    mpfr_init_set_si(yy_i, (long) yy[i_], MPFR_RNDN);
	    R_mpfr_mod(r, x_i, yy_i, MPFR_RNDN);
	    mpfr_sub (x_i, x_i,  r, MPFR_RNDN); // x~ = x - r =  x - (x %% y)
	    mpfr_div (x_i, x_i,yy_i,MPFR_RNDN); // x = x~ / y = (x - (x %% y))/y
	    mpfr_clear(r); mpfr_clear(yy_i); break;
	}

	case  7: /*  /  */ mpfr_div_si(x_i, x_i, (long) yy[i_], MPFR_RNDN); break;

	default:
	    error("invalid op code (%d) in Arith_mpfr", i_op);
	}

	SET_VECTOR_ELT(val, i, MPFR_as_R(x_i));
    }
    MISMATCH_WARN;

    mpfr_clear (x_i);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(2);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Arith_mpfr_i */

SEXP Arith_i_mpfr(SEXP x, SEXP y, SEXP op)
{
#ifdef using_Data_slot
    SEXP yD = PROTECT(GET_SLOT(y, Rmpfr_Data_Sym));
#else
# define yD y
#endif
    int *xx = INTEGER(x);
    int nx = length(x), ny = length(yD), i_op = asInteger(op), i,
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;

    if(TYPEOF(x) != INTSXP)
	error("Arith[%d](i,mpfr): 'x' is not a \"integer\"", i_op);

    SEXP val = PROTECT(allocVector(VECSXP, n));
    mpfr_t y_i;
    mpfr_init(y_i); /* with default precision */

    SET_MISMATCH;
    for(i=0; i < n; i++) {
	int i_ = i % nx;
	R_asMPFR(VECTOR_ELT(yD, i % ny), y_i);
	switch(i_op) {
	    /* Note we assign use y_i as "input and output" ==> *same*
	       precision, even though in some cases the result may
	       need higher precision */
	case  1: /*  +  */ mpfr_add_si(y_i, y_i, (long) xx[i_], MPFR_RNDN); break;
	case  2: /*  -  */ mpfr_si_sub(y_i, (long) xx[i_], y_i, MPFR_RNDN); break;
	case  3: /*  *  */ mpfr_mul_si(y_i, y_i, (long) xx[i_], MPFR_RNDN); break;
	case  4: /*  ^  */ {
#define R_MPFR_SI_POW(_XXI, _YI)					\
	    long _x = (long) _XXI;					\
	    if(_x >= 0)							\
		mpfr_ui_pow(_YI, (unsigned long) _x, _YI, MPFR_RNDN);	\
	    else if(mpfr_integer_p(_YI)) { /* <neg. x> ^ <integer> */	\
		mpfr_ui_pow(_YI, (unsigned long) -_x, _YI, MPFR_RNDN);	\
		mpfr_neg(_YI, _YI, MPFR_RNDN);				\
	    }								\
	    else /* <neg. x>  ^  <non-integer>  |-> NaN : */		\
		mpfr_set_nan (_YI);					\
	    break

	    R_MPFR_SI_POW(xx[i_], y_i);
	}
	case  5: /* %%  */ {
	    mpfr_t xx_i, r;
	    mpfr_init_set_si(xx_i, (long) xx[i_], MPFR_RNDN);
	    mpfr_init(r);
	    R_mpfr_mod(r, xx_i, y_i, MPFR_RNDN);
	    mpfr_set(y_i, r, MPFR_RNDN);
	    mpfr_clear(r); mpfr_clear(xx_i); break;
	}
	case  6: /* %/% */ {
	    mpfr_t r, xx_i; mpfr_init(r);
	    mpfr_prec_t y_prec = mpfr_get_prec(y_i);
	    if(mpfr_get_prec(r) < y_prec)
		mpfr_set_prec (r, y_prec);
	    mpfr_init_set_si(xx_i, (long) xx[i_], MPFR_RNDN);
	    R_mpfr_mod(r, xx_i, y_i, MPFR_RNDN);
	    mpfr_sub (xx_i,xx_i, r, MPFR_RNDN); // x~ = x - r =  x - (x %% y)
	    mpfr_div (y_i, xx_i,y_i,MPFR_RNDN); // y = x~ / y = (x - (x %% y))/y
	    mpfr_clear(r); mpfr_clear(xx_i); break;
	}
	case  7: /*  /  */ mpfr_si_div(y_i, (long) xx[i_], y_i, MPFR_RNDN); break;

	default:
	    error("invalid op code (%d) in Arith_mpfr", i_op);
	}

	SET_VECTOR_ELT(val, i, MPFR_as_R(y_i));
    }
    MISMATCH_WARN;

    mpfr_clear (y_i);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(2);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Arith_i_mpfr */

SEXP Arith_mpfr_d(SEXP x, SEXP y, SEXP op)
{
#ifdef using_Data_slot
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));
#else
# define xD x
#endif
    double *yy = REAL(y);
    int nx = length(xD), ny = length(y), i_op = asInteger(op), i,
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;

    if(TYPEOF(y) != REALSXP)
	error("Arith[%d](mpfr,d): 'y' is not a \"double\"", i_op);

    SEXP val = PROTECT(allocVector(VECSXP, n));
    mpfr_t x_i, yy_i;

    mpfr_init(x_i);
    mpfr_init(yy_i); /* with default precision */

    SET_MISMATCH;
    for(i=0; i < n; i++) {
	double yi = yy[i % ny];
	int y_is_int = (yi == trunc(yi) && LONG_MIN <= yi && yi <= LONG_MAX);

	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);
	if(y_is_int) { /* can use <mpfr> o <integer>  routines */
	    switch(i_op) {
	    case  1: /*  +  */ mpfr_add_si(x_i, x_i, (long)yi, MPFR_RNDN); break;
	    case  2: /*  -  */ mpfr_sub_si(x_i, x_i, (long)yi, MPFR_RNDN); break;
	    case  3: /*  *  */ mpfr_mul_si(x_i, x_i, (long)yi, MPFR_RNDN); break;
	    case  4: /*  ^  */ mpfr_pow_si(x_i, x_i, (long)yi, MPFR_RNDN); break;
	    case  5: /* %%  */ {
		mpfr_set_si(yy_i, (long)yi, MPFR_RNDN);
		R_mpfr_mod(x_i, x_i, yy_i, MPFR_RNDN);
		break;
	    }
	    case  6: /* %/% */ {
		mpfr_t r;
		mpfr_init(r);
		mpfr_prec_t x_prec = mpfr_get_prec(x_i);
		if(mpfr_get_prec(r) < x_prec)
		    mpfr_set_prec (r, x_prec);

		mpfr_set_si(yy_i, (long) yi, MPFR_RNDN);
		R_mpfr_mod(r,  x_i, yy_i, MPFR_RNDN);
		mpfr_sub (x_i, x_i,  r,   MPFR_RNDN); // x~ = x - r =  x - (x %% y)
		mpfr_div (x_i, x_i, yy_i, MPFR_RNDN); // x = x~ / y = (x - (x %% y))/y
		mpfr_clear(r); break;
	    }

	    case  7: /*  /  */ mpfr_div_si(x_i, x_i, (long)yi, MPFR_RNDN); break;
	    default:
		error("invalid op code (%d) in Arith_mpfr_d", i_op);
	    }
	}
	else {
	    mpfr_set_d (yy_i, yi, MPFR_RNDD);
	    switch(i_op) {
		/* Note we assign use x_i as "input and output" ==> *same*
		   precision, even though in some cases the result may
		   need higher precision */
	    case  1: /*  +  */ mpfr_add(x_i, x_i, yy_i, MPFR_RNDN); break;
	    case  2: /*  -  */ mpfr_sub(x_i, x_i, yy_i, MPFR_RNDN); break;
	    case  3: /*  *  */ mpfr_mul(x_i, x_i, yy_i, MPFR_RNDN); break;
	    case  4: /*  ^  */ mpfr_pow(x_i, x_i, yy_i, MPFR_RNDN); break;
	    case  5: /* %% */ R_mpfr_mod(x_i, x_i, yy_i, MPFR_RNDN); break;
	    case  6: /* %/% */ {
		mpfr_t r; mpfr_init(r);
		mpfr_prec_t x_prec = mpfr_get_prec(x_i);
		if(mpfr_get_prec(r) < x_prec)
		    mpfr_set_prec (r, x_prec);
		R_mpfr_mod(r,  x_i, yy_i, MPFR_RNDN);
		mpfr_sub (x_i, x_i,  r,   MPFR_RNDN); // x~ = x - r =  x - (x %% y)
		mpfr_div (x_i, x_i, yy_i, MPFR_RNDN); // x = x~ / y = (x - (x %% y))/y
		mpfr_clear(r); break;
	    }
	    case  7: /*  /  */ mpfr_div(x_i, x_i, yy_i, MPFR_RNDN); break;
	    default:
		error("invalid op code (%d) in Arith_mpfr_d", i_op);
	    }
	}
	SET_VECTOR_ELT(val, i, MPFR_as_R(x_i));
    }
    MISMATCH_WARN;

    mpfr_clear (x_i); mpfr_clear (yy_i);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(2);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Arith_mpfr_d */

SEXP Arith_d_mpfr(SEXP x, SEXP y, SEXP op)
{
#ifdef using_Data_slot
    SEXP yD = PROTECT(GET_SLOT(y, Rmpfr_Data_Sym));
#else
# define yD y
#endif
    double *xx = REAL(x);
    int nx = length(x), ny = length(yD), i_op = asInteger(op), i,
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;

    if(TYPEOF(x) != REALSXP)
	error("Arith[%d](d,mpfr): 'x' is not a \"double\"", i_op);

    SEXP val = PROTECT(allocVector(VECSXP, n));
    mpfr_t y_i;
    mpfr_init(y_i);

    SET_MISMATCH;
    for(i=0; i < n; i++) {
	double xi = xx[i % nx];
	int x_is_int = (xi == trunc(xi) && LONG_MIN <= xi && xi <= LONG_MAX);

	R_asMPFR(VECTOR_ELT(yD, i % ny), y_i);
	if(x_is_int) { /* can use  <integer> o <mpfr>  routines */
/* 	    REprintf("x[i] (= %g) is int: (long)* = %ld\n", xi, (long)xi); */
	    switch(i_op) {
	    case  1: /*  +  */ mpfr_add_si(y_i, y_i, (long)xi, MPFR_RNDN); break;
	    case  2: /*  -  */ mpfr_si_sub(y_i, (long)xi, y_i, MPFR_RNDN); break;
	    case  3: /*  *  */ mpfr_mul_si(y_i, y_i, (long)xi, MPFR_RNDN); break;
	    case  4: /*  ^  */ {
		R_MPFR_SI_POW((long)xi, y_i);
	    }
	    case  5: /* %%  */ {
		mpfr_t xx_i, r;
		mpfr_init_set_si(xx_i, (long)xi, MPFR_RNDN);
		mpfr_init(r);
		R_mpfr_mod(r, xx_i, y_i, MPFR_RNDN);
		mpfr_set(y_i, r, MPFR_RNDN);
		mpfr_clear(r); mpfr_clear(xx_i); break;
	    }
	    case  6: /* %/% */ {
		mpfr_t r, xx_i; mpfr_init(r);
		mpfr_prec_t y_prec = mpfr_get_prec(y_i);
		if(mpfr_get_prec(r) < y_prec)
		    mpfr_set_prec (r, y_prec);
		mpfr_init_set_si(xx_i, (long) xi, MPFR_RNDN);
		R_mpfr_mod(r,  xx_i, y_i, MPFR_RNDN);
		mpfr_sub (xx_i,xx_i,  r,  MPFR_RNDN); // x~ = x - r =  x - (x %% y)
		mpfr_div (y_i, xx_i, y_i, MPFR_RNDN); // y = x~ / y = (x - (x %% y))/y
		mpfr_clear(r); mpfr_clear(xx_i); break;
	    }
	    case  7: /*  /  */ mpfr_si_div(y_i, (long)xi, y_i, MPFR_RNDN); break;
	    default:
		error("invalid op code (%d) in Arith_d_mpfr", i_op);
	    }
	}
	else {
	    mpfr_t xx_i;
	    mpfr_init_set_d (xx_i, xi, MPFR_RNDD);
	    switch(i_op) {
		/* Note we assign use y_i as "input and output" ==> *same*
		   precision, even though in some cases the result may
		   need higher precision */
	    case  1: /*  +  */ mpfr_add(y_i, xx_i, y_i, MPFR_RNDN); break;
	    case  2: /*  -  */ mpfr_sub(y_i, xx_i, y_i, MPFR_RNDN); break;
	    case  3: /*  *  */ mpfr_mul(y_i, xx_i, y_i, MPFR_RNDN); break;
	    case  4: /*  ^  */ mpfr_pow(y_i, xx_i, y_i, MPFR_RNDN); break;
	    case  5: /* %%  */ {
		mpfr_t r; mpfr_init(r);
		R_mpfr_mod(r, xx_i, y_i, MPFR_RNDN);
		mpfr_set(y_i, r, MPFR_RNDN);
		mpfr_clear(r); break;
	    }
	    case  6: /* %/% */ {
		mpfr_t r; mpfr_init(r);
		mpfr_prec_t y_prec = mpfr_get_prec(y_i);
		if(mpfr_get_prec(r) < y_prec)
		    mpfr_set_prec (r, y_prec);
		R_mpfr_mod(r,  xx_i, y_i, MPFR_RNDN);
		mpfr_sub (xx_i,xx_i,  r,  MPFR_RNDN); // x~ = x - r =  x - (x %% y)
		mpfr_div (y_i, xx_i, y_i, MPFR_RNDN); // y = x~ / y = (x - (x %% y))/y
		mpfr_clear(r); break;
	    }
	    case  7: /*  /  */ mpfr_div(y_i, xx_i, y_i, MPFR_RNDN); break;
	    default:
		error("invalid op code (%d) in Arith_d_mpfr", i_op);
	    }
	    mpfr_clear(xx_i);
	}
	SET_VECTOR_ELT(val, i, MPFR_as_R(y_i));
    }
    MISMATCH_WARN;

    mpfr_clear (y_i);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(2);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Arith_d_mpfr */





SEXP Compare_mpfr(SEXP x, SEXP y, SEXP op)
{
#ifdef using_Data_slot
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym)),
	yD = PROTECT(GET_SLOT(y, Rmpfr_Data_Sym));
#else
# define xD x
# define yD y
#endif
    int nx = length(xD), ny = length(yD), i_op = asInteger(op), i,
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;

    SEXP val = PROTECT(allocVector(LGLSXP, n));
    int *r = LOGICAL(val);
    mpfr_t x_i, y_i;
    mpfr_init(x_i); /* with default precision */
    mpfr_init(y_i); /* with default precision */

    SET_MISMATCH;
    for(i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);
	R_asMPFR(VECTOR_ELT(yD, i % ny), y_i);

	switch(i_op) {
	case 1: /* == */ r[i] = mpfr_equal_p(x_i, y_i); break;
	case 2: /* >  */ r[i] = mpfr_greater_p(x_i, y_i); break;
	case 3: /* <  */ r[i] = mpfr_less_p(x_i, y_i); break;
	case 4: /* != */ r[i] = mpfr_lessgreater_p(x_i, y_i); break;
	case 5: /* <= */ r[i] = mpfr_lessequal_p(x_i, y_i); break;
	case 6: /* >= */ r[i] = mpfr_greaterequal_p(x_i, y_i); break;
	default:
	    error("invalid op code (%d) in Compare_mpfr", i_op);
	}
    }
    MISMATCH_WARN;

    mpfr_clear (x_i); mpfr_clear (y_i);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(3);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Compare_mpfr */

SEXP Compare_mpfr_i(SEXP x, SEXP y, SEXP op)
{
#ifdef using_Data_slot
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));
#else
# define xD x
#endif
    int *yy = INTEGER(y);
    int nx = length(xD), ny = length(y), i_op = asInteger(op), i,
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;

    SEXP val = PROTECT(allocVector(LGLSXP, n));
    int *r = LOGICAL(val);
    mpfr_t x_i;
    mpfr_init(x_i);

    SET_MISMATCH;
    for(i=0; i < n; i++) {
	int yi = yy[i % ny], c;
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);
	c = mpfr_cmp_si(x_i, (long) yi);/* gives  c > or == or < 0 */
	if(c == 0 && /* also includes case where an operand is NaN */
	   (yi == NA_INTEGER || mpfr_nan_p(x_i))) {
	    r[i] = NA_LOGICAL;
	} else {
	    switch(i_op) {
	    case 1: /* == */ r[i] = (c == 0); break;
	    case 2: /* >  */ r[i] = (c >  0); break;
	    case 3: /* <  */ r[i] = (c <  0); break;
	    case 4: /* != */ r[i] = (c != 0); break;
	    case 5: /* <= */ r[i] = (c <= 0); break;
	    case 6: /* >= */ r[i] = (c >= 0); break;
	    default:
		error("invalid op code (%d) in Compare_mpfr_i", i_op);
	    }
	}
    }
    MISMATCH_WARN;

    mpfr_clear (x_i);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(2);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Compare_mpfr_i */

SEXP Compare_mpfr_d(SEXP x, SEXP y, SEXP op)
{
#ifdef using_Data_slot
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));
#else
# define xD x
#endif
    double *yy = REAL(y);
    int nx = length(xD), ny = length(y), i_op = asInteger(op), i,
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny), mismatch = 0;

    SEXP val = PROTECT(allocVector(LGLSXP, n));
    int *r = LOGICAL(val);
    mpfr_t x_i;
    mpfr_init(x_i);

    SET_MISMATCH;
    for(i=0; i < n; i++) {
	double yi = yy[i % ny];
	int c;
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);
	c = mpfr_cmp_d(x_i, yi);/* gives  c > or == or < 0 */
	if(c == 0 && /* also includes case where an operand is NaN */
	   (ISNAN(yi) || mpfr_nan_p(x_i))) {
	    r[i] = NA_LOGICAL;
	} else {
	    switch(i_op) {
	    case 1: /* == */ r[i] = (c == 0); break;
	    case 2: /* >  */ r[i] = (c >  0); break;
	    case 3: /* <  */ r[i] = (c <  0); break;
	    case 4: /* != */ r[i] = (c != 0); break;
	    case 5: /* <= */ r[i] = (c <= 0); break;
	    case 6: /* >= */ r[i] = (c >= 0); break;
	    default:
		error("invalid op code (%d) in Compare_mpfr_d", i_op);
	    }
	}
    }
    MISMATCH_WARN;

    mpfr_clear (x_i);
    mpfr_free_cache();
#ifdef using_Data_slot
    UNPROTECT(2);
#else
    UNPROTECT(1);
#endif
    return val;
} /* Compare_mpfr_d */





#ifdef __NOT_ANY_MORE__

/* Not really used anymore : */

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


SEXP exp_mpfr1(SEXP x)
{
    SEXP val;
    INIT_1_SETUP(x, r);
    mpfr_exp(r, r, MPFR_RNDN);
    /*       -  -  ((result may need higher precision)) .. */
    FINISH_1_RETURN(r, val);
}

SEXP log_mpfr1(SEXP x)
{
    SEXP val;
    INIT_1_SETUP(x, r); mpfr_log(r, r, MPFR_RNDN);
    FINISH_1_RETURN(r, val);
}


/* Unused */

#define INIT_2_SETUP(_X_, _R_, _S_)		\
    mpfr_t _R_, _S_;				\
						\
    mpfr_init2(_R_, R_mpfr_prec(_X_));		\
    /* _S_ should get same precision as _R_ :*/	\
    mpfr_init2(_S_, mpfr_get_prec(_R_));	\
    R_asMPFR(_X_, _R_)

#define FINISH_2_RETURN(_R_, _S_, val)		\
    val = PROTECT(MPFR_as_R(_R_));		\
    mpfr_clear(_R_); mpfr_clear(_S_);		\
    mpfr_free_cache();				\
    UNPROTECT(1);				\
    return val


#endif
/* __NOT_ANY_MORE__ */
