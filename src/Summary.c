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

#define R_mpfr_prec(x) INTEGER(GET_SLOT(x, Rmpfr_precSym))[0]

SEXP Summary_mpfr(SEXP x, SEXP na_rm, SEXP op)
{
/* FIXME: use  enum type for the op codes */

    SEXP D = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));/* an R list() of length */
    mp_prec_t current_prec = mpfr_get_default_prec();
    int n = length(D), i_op = asInteger(op),
	return_list = (i_op <= 5),
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

    case 1: /* max */ mpfr_set_inf(Summ, -1);/* := -Inf */; Rmpfr_set(1);
    case 2: /* min */ mpfr_set_inf(Summ, +1);/* := +Inf */; Rmpfr_set(1);
    case 3: /* range */
	mpfr_init(Sum2);
	mpfr_set_inf(Summ, +1);/* := +Inf for min() */
	mpfr_set_inf(Sum2, -1);/* := -Inf for max() */
	Rmpfr_set(2);
    case 4: /* prod */mpfr_set_d (Summ, 1., GMP_RNDZ); Rmpfr_set(1);
    case 5: /* sum */ mpfr_set_d (Summ, 0., GMP_RNDZ); Rmpfr_set(1);

    case 6: /* any */ val = ScalarLogical(FALSE); break;
    case 7: /* all */ val = ScalarLogical(TRUE);  break;

    default:
	UNPROTECT(1);
	error("invalid op code (%d) in Summary_mpfr", i_op);
    }

    if (!return_list) /* will return logical */
	ans = LOGICAL(val);

    for(i=0; i < n; i++) {
	SEXP Di = VECTOR_ELT(D, i);
	R_asMPFR(Di, R_i);

	if(mpfr_nan_p(R_i)) { /* handling does not depend on i_op */
/* 	    REprintf("Summary_mpfr(), i=%d :  R_i is NaN\n", i); */
	    if(remove_na) /* skip this NA / NAN entry: */
		continue;
	    else { /* result is NA */
		/* should not be needed: R_i *is* NaN already :
		   mpfr_set_nan(R_i); */
		switch(i_op) {
		case 1: /* max */
		case 2: /* min */
		case 4: /* prod */
		case 5: /* sum */
		    SET_VECTOR_ELT(val, 0, Di);
		    break;
		case 3: /* range */
		    SET_VECTOR_ELT(val, 0, Di);
		    SET_VECTOR_ELT(val, 1, Di);
		    break;
		    /*---------------------------------------------*/
		case 6: /* any */
		    if(*ans == FALSE) *ans = NA_LOGICAL;
		    break;
		case 7: /* all */
		    if(*ans == TRUE) *ans = NA_LOGICAL;
		    break;
		}
	    }
	    if(i_op <= 5) { /* return()  *unless* for  any()/all() : */
		mpfr_free_cache();
		UNPROTECT(2);
		return val;
	    }
	}
	else n_valid++;

	if(return_list) { /* hence using  Summ */
	    mp_prec_t i_prec = mpfr_get_prec(R_i);
	    if(current_prec < i_prec) /* increase precision */ {
		current_prec = i_prec;
		mpfr_prec_round(Summ, i_prec, GMP_RNDN);
		if(i_op == 3)
		    mpfr_prec_round(Sum2, i_prec, GMP_RNDN);
	    }
	}

	switch(i_op) {
	    /* Note we assign use R_i as "input and output" ==> *same*
	       precision, even though in some cases the result may
	       need higher precision */

	case 1: /* max */ mpfr_max(Summ, Summ, R_i, GMP_RNDN); break;
	case 2: /* min */ mpfr_min(Summ, Summ, R_i, GMP_RNDN); break;
	case 3: /* range */
	    mpfr_min(Summ, Summ, R_i, GMP_RNDN);
	    mpfr_max(Sum2, Sum2, R_i, GMP_RNDN);
	    break;

	case 4: /* prod */ mpfr_mul(Summ, Summ, R_i, GMP_RNDN); break;
	case 5: /* sum */  mpfr_add(Summ, Summ, R_i, GMP_RNDN); break;

	case 6: /* any */ if(!mpfr_zero_p(R_i)) *ans = TRUE; break;
	case 7: /* all */ if( mpfr_zero_p(R_i)) *ans = FALSE; break;

	}
    } /* for(i .. n) */

    mpfr_clear (R_i);
    switch(i_op) {
    case 1: /* max */
    case 2: /* min */
    case 4: /* prod */
    case 5: /* sum */
	SET_VECTOR_ELT(val, 0, MPFR_as_R(Summ));
	mpfr_clear (Summ);
	break;
    case 3: /* range */
	SET_VECTOR_ELT(val, 0, MPFR_as_R(Summ));
	SET_VECTOR_ELT(val, 1, MPFR_as_R(Sum2));
	mpfr_clear (Summ);
	mpfr_clear (Sum2);
	break;
    case 6: /* any */
    case 7: /* all */
	/* nothing to be done */
	break;
    }

    mpfr_free_cache();
    UNPROTECT(return_list ? 2 : 1);
    return val;
} /* Summary_mpfr() */
