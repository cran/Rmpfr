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

/*------------------------------------------------------------------------*/

SEXP Summary_mpfr(SEXP x, SEXP na_rm, SEXP op)
{
    enum { MAX = 1, MIN, RANGE, PROD, SUM, ANY = 10, ALL } i_op = asInteger(op);
/* MUST be sync'ed with  ../R/Summary.R
 *                       ~~~~~~~~~~~~~~ where  .Summary.codes <-
 * c("max" = 1,  "min" = 2, "range" = 3, "prod" = 4, "sum" = 5,
 *   "any" = 10, "all" = 11)
*/
    mpfr_prec_t current_prec = mpfr_get_default_prec();
    int n = length(x),
	return_list = (i_op < ANY),
	remove_na = asLogical(na_rm), n_valid = 0, i;

    SEXP val = R_NilValue;
    int *ans = NULL; /*"Wall", will be := LOGICAL(val)   for any() / all() only */
    mpfr_t R_i,
	Summ, Sum2; /* for range(), max(), min(), sum(), prod() */

    mpfr_init(R_i); /* with default precision */
    if (return_list)
	mpfr_init(Summ);

#define Rmpfr_set(_N_) val = PROTECT(allocVector(VECSXP, _N_)); break

    switch(i_op) {

    case MAX: mpfr_set_inf(Summ, -1);/* := -Inf */; Rmpfr_set(1);
    case MIN: mpfr_set_inf(Summ, +1);/* := +Inf */; Rmpfr_set(1);
    case RANGE:
	mpfr_init(Sum2);
	mpfr_set_inf(Summ, +1);/* := +Inf for min() */
	mpfr_set_inf(Sum2, -1);/* := -Inf for max() */
	Rmpfr_set(2);
    case PROD: mpfr_set_d (Summ, 1., MPFR_RNDZ); Rmpfr_set(1);
    case SUM: mpfr_set_d (Summ, 0., MPFR_RNDZ); Rmpfr_set(1);

    case ANY: val = ScalarLogical(FALSE); break;
    case ALL: val = ScalarLogical(TRUE);  break;

    default:
	error("invalid op code (%d) in Summary_mpfr", i_op);
    }

    if (!return_list) /* will return logical */
	ans = LOGICAL(val);

    for(i=0; i < n; i++) {
	SEXP xi = VECTOR_ELT(x, i);
	R_asMPFR(xi, R_i);

	if(mpfr_nan_p(R_i)) { /* handling does not depend on i_op */
/* 	    REprintf("Summary_mpfr(), i=%d :  R_i is NaN\n", i); */
	    if(remove_na) /* skip this NA / NAN entry: */
		continue;
	    else { /* result is NA */
		/* should not be needed: R_i *is* NaN already :
		   mpfr_set_nan(R_i); */
		switch(i_op) {
		case MAX:
		case MIN:
		case PROD:
		case SUM:
		    SET_VECTOR_ELT(val, 0, xi);
		    break;
		case RANGE:
		    SET_VECTOR_ELT(val, 0, xi);
		    SET_VECTOR_ELT(val, 1, xi);
		    break;
		    /*---------------------------------------------*/
		case ANY:
		    if(*ans == FALSE) *ans = NA_LOGICAL;
		    break;
		case ALL:
		    if(*ans == TRUE) *ans = NA_LOGICAL;
		    break;
		}
	    }
	    if(return_list) { /* return()  *unless* for  any()/all() : */
		mpfr_free_cache();
		UNPROTECT(1);
		return val;
	    }
	}
	else n_valid++;

	if(return_list) { /* hence using  Summ */
	    mpfr_prec_t i_prec = mpfr_get_prec(R_i);
	    if(current_prec < i_prec) /* increase precision */ {
		current_prec = i_prec;
		mpfr_prec_round(Summ, i_prec, MPFR_RNDN);
		if(i_op == RANGE)
		    mpfr_prec_round(Sum2, i_prec, MPFR_RNDN);
	    }
	}

	switch(i_op) {
	    /* Note we assign use R_i as "input and output" ==> *same*
	       precision, even though in some cases the result may
	       need higher precision */

	case MAX: mpfr_max(Summ, Summ, R_i, MPFR_RNDN); break;
	case MIN: mpfr_min(Summ, Summ, R_i, MPFR_RNDN); break;
	case RANGE:
	    mpfr_min(Summ, Summ, R_i, MPFR_RNDN);
	    mpfr_max(Sum2, Sum2, R_i, MPFR_RNDN);
	    break;

	case PROD: mpfr_mul(Summ, Summ, R_i, MPFR_RNDN); break;
	case SUM:  mpfr_add(Summ, Summ, R_i, MPFR_RNDN); break;

	case ANY: if(!mpfr_zero_p(R_i)) *ans = TRUE; break;
	case ALL: if( mpfr_zero_p(R_i)) *ans = FALSE; break;

	}
    } /* for(i .. n) */

    mpfr_clear (R_i);
    switch(i_op) {
    case MAX:
    case MIN:
    case PROD:
    case SUM:
	SET_VECTOR_ELT(val, 0, MPFR_as_R(Summ));
	mpfr_clear (Summ);
	break;
    case RANGE:
	SET_VECTOR_ELT(val, 0, MPFR_as_R(Summ));
	SET_VECTOR_ELT(val, 1, MPFR_as_R(Sum2));
	mpfr_clear (Summ);
	mpfr_clear (Sum2);
	break;
    case ANY:
    case ALL:
	/* nothing to be done */
	break;
    }

    mpfr_free_cache();
    if(return_list) UNPROTECT(1);
    return val;
} /* Summary_mpfr() */

/**
   Compute sum(x * y)  ==  x %*% y   for two [mpfr-]vectors of the same length
   Both x and y can be in  {mpfr, double, integer} !
 */
SEXP R_mpfr_sumprod(SEXP x, SEXP y, SEXP minPrec, SEXP alternating_)
{
    int n = length(x);
    if(length(y) != n)
	error("%d == length(x) != length(y) == %d", n, length(y));
    int i_min_prec = asInteger(minPrec), nprot = 1;
    Rboolean alternating = asLogical(alternating_);
    // Simplification (FIXME: more efficient -> use 6 cases; s/ M_n / M_d and M_i /)
    if(isInteger(x)) { PROTECT(x = coerceVector(x, REALSXP)); nprot++; }
    if(isInteger(y)) { PROTECT(y = coerceVector(y, REALSXP)); nprot++; }
    if(isReal(x) && isReal(y))
	error("R_mpfr_sumprod(x,y, .): either x or y must be \"mpfr\", but both are numeric");
    // --> three cases:
    //  M_M: both mpfr,
    //  n_M: x numeric, y mpfr
    //  M_n: x mpfr   , y numeric
    enum { M_M, n_M, M_n } R_case = isReal(x) ? n_M : isReal(y) ? M_n : M_M;
    Rboolean use_r = alternating && R_case == M_M;

    mpfr_t Summ, x_i, y_i, r;
    mpfr_inits(Summ, x_i, y_i, (mpfr_ptr) 0); /* with default precision */
    mpfr_set_d(Summ, 0., MPFR_RNDZ);
    double *xx = NULL, *yy = NULL;
    if(R_case == n_M)
	xx = REAL(x);
    else if(R_case == M_n)
	yy = REAL(y);

    mpfr_prec_t min_prec = MPFR_PREC_MIN, xy_prec, S_prec = mpfr_get_prec(Summ);
    if(i_min_prec != NA_INTEGER && min_prec < i_min_prec)
	min_prec = (mpfr_prec_t) i_min_prec;
    if(S_prec < min_prec) {
	mpfr_prec_round (Summ, min_prec, MPFR_RNDN);
	S_prec = min_prec;
    }
    if(use_r)
	mpfr_init2(r, S_prec);

    for(int i=0; i < n; i++)
    {
	Rboolean NA_res;
	double xi = 0., yi = 0.; // Wall
	switch(R_case) {
	case M_M :
	    R_asMPFR(VECTOR_ELT(x, i), x_i);
	    R_asMPFR(VECTOR_ELT(y, i), y_i);
	    NA_res = (mpfr_nan_p(x_i) || mpfr_nan_p(y_i));
	    xy_prec = imax2(mpfr_get_prec(x_i),
			    mpfr_get_prec(y_i));
	    break;
	case M_n :
	    R_asMPFR(VECTOR_ELT(x, i), x_i);
	    yi = yy[i];
	    NA_res = (mpfr_nan_p(x_i) || ISNA(yi));
	    xy_prec = imax2(mpfr_get_prec(x_i), 53);
	    break;
	case n_M :
	    xi = xx[i];
	    R_asMPFR(VECTOR_ELT(y, i), y_i);
	    NA_res = (ISNA(xi) || mpfr_nan_p(y_i));
	    xy_prec = imax2(53, mpfr_get_prec(y_i));
	    break;
	} // switch()

	if(NA_res) {
	    mpfr_set_nan(Summ);
	    continue; // no need to continue the loop
	}

	if(S_prec < xy_prec) {/* increase it, since it will store the result */
	    mpfr_prec_round (Summ, xy_prec, MPFR_RNDN);
	    S_prec = xy_prec;
	    if(use_r) mpfr_set_prec(r, S_prec);
	}

	if(alternating && (i % 2)) { // Summ := Summ - (x_i * y_i)
	    switch(R_case) {
	    case M_M :
		/* mpfr_fms (ROP, OP1, OP2, OP3, RND)
		 * Set ROP to (OP1 times OP2) - OP3 rounded in the direction RND. */
		mpfr_fms (r, x_i, y_i, Summ, MPFR_RNDN); // r = x_i * y_i - Summ
		mpfr_neg (Summ, r, MPFR_RNDN);
		break;
	    case M_n :
		mpfr_mul_d(x_i, x_i, yi, MPFR_RNDN);
		mpfr_sub(Summ, Summ, x_i,MPFR_RNDN);
		break;
	    case n_M :
		mpfr_mul_d(y_i, y_i, xi, MPFR_RNDN);
		mpfr_sub(Summ, Summ, y_i,MPFR_RNDN);
		break;
	    }
	}
	else { // Summ := Summ + (x_i * y_i)
	    switch(R_case) {
	    case M_M :
		/* mpfr_fma (ROP, OP1, OP2, OP3, RND)
		 * Set ROP to (OP1 times OP2) + OP3  rounded in the direction RND. */
		mpfr_fma (Summ, x_i, y_i, Summ, MPFR_RNDN);
		break;
	    case M_n :
		mpfr_mul_d(x_i, x_i, yi, MPFR_RNDN);
		mpfr_add(Summ, Summ, x_i,MPFR_RNDN);
		break;
	    case n_M :
		mpfr_mul_d(y_i, y_i, xi, MPFR_RNDN);
		mpfr_add(Summ, Summ, y_i,MPFR_RNDN);
		break;
	    }
	}
    } // for( i )

    // val <- list( Summ ) :
    SEXP val = PROTECT(allocVector(VECSXP, 1));
    SET_VECTOR_ELT(val, 0, MPFR_as_R(Summ));

    mpfr_clears(Summ, x_i, y_i, (mpfr_ptr) 0);
    if(use_r) mpfr_clear(r);
    mpfr_free_cache();
    UNPROTECT(nprot);
    return val;
} // R_mpfr_sumprod
