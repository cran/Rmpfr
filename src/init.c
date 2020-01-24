/* Setup here "copied" from Matrix package */

#include <R_ext/Rdynload.h>

#define _in_Rmpfr_init_

#include "Rmpfr_utils.h"
#include "Syms.h"

#undef _in_Rmpfr_init_

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

    CALLDEF(d2mpfr1, 3),
    CALLDEF(d2mpfr1_list, 3),
#ifdef Have_interface_Rmpfr_gmp
    CALLDEF(mpz2mpfr1, 3),
    CALLDEF(mpz2mpfr1_list, 3),
#endif
#ifdef R_had_R_Outputfile_in_API
#ifndef WIN32
    /* only works on "unix-alikes" */
    CALLDEF(print_mpfr, 2),
    CALLDEF(print_mpfr1, 2),
#endif
#endif
    CALLDEF(mpfr2d, 2),
    CALLDEF(mpfr2i, 2),
    CALLDEF(mpfr2str, 4),
    CALLDEF(R_mpfr_formatinfo, 1),
    CALLDEF(R_mpfr_2exp, 1),
    CALLDEF(str2mpfr1_list, 4),

    CALLDEF(Rmpfr_minus, 1),
    CALLDEF(Rmpfr_abs, 1),

    CALLDEF(Math_mpfr, 2),
    CALLDEF(Arith_mpfr, 3),
    CALLDEF(Arith_mpfr_i, 3),
    CALLDEF(Arith_i_mpfr, 3),
    CALLDEF(Arith_mpfr_d, 3),
    CALLDEF(Arith_d_mpfr, 3),

    CALLDEF(Compare_mpfr, 3),
    CALLDEF(Compare_mpfr_i, 3),
    CALLDEF(Compare_mpfr_d, 3),

    CALLDEF(Summary_mpfr, 3),
    CALLDEF(R_mpfr_sumprod, 4),

    CALLDEF(R_mpfr_set_debug, 1),
    CALLDEF(R_mpfr_set_default_prec, 1),
    CALLDEF(R_mpfr_get_default_prec, 0),
    CALLDEF(R_mpfr_prec_range, 1),
    CALLDEF(R_mpfr_get_erange, 1),
    CALLDEF(R_mpfr_set_erange, 2),
    CALLDEF(R_mpfr_erange_int_p, 0),

    CALLDEF(R_mpfr_get_version, 0),
    CALLDEF(R_mpfr_get_GMP_numb_bits, 0),

    CALLDEF(const_asMpfr, 3),

    CALLDEF(R_mpfr_is_finite, 1),	CALLDEF(R_mpfr_is_finite_A, 1),
    CALLDEF(R_mpfr_is_infinite, 1),	CALLDEF(R_mpfr_is_infinite_A, 1),
    CALLDEF(R_mpfr_is_integer, 1),	CALLDEF(R_mpfr_is_integer_A, 1),
    CALLDEF(R_mpfr_is_na, 1),		CALLDEF(R_mpfr_is_na_A, 1),
    CALLDEF(R_mpfr_is_zero, 1),      	CALLDEF(R_mpfr_is_zero_A, 1),

    CALLDEF(R_mpfr_atan2, 3),
    CALLDEF(R_mpfr_igamma, 3),
    CALLDEF(R_mpfr_hypot, 3),
    CALLDEF(R_mpfr_beta, 3),
    CALLDEF(R_mpfr_lbeta, 3),

    CALLDEF(R_mpfr_jn, 3),
    CALLDEF(R_mpfr_yn, 3),
    CALLDEF(R_mpfr_fac, 3),
    CALLDEF(R_mpfr_choose, 3),
    CALLDEF(R_mpfr_poch, 3),
    CALLDEF(R_mpfr_round, 3),

    {NULL, NULL, 0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_Rmpfr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

#define RREGDEF(name)  R_RegisterCCallable("Rmpfr", #name, (DL_FUNC) name)

    RREGDEF(d2mpfr1);
    RREGDEF(d2mpfr1_list);
#ifdef Have_interface_Rmpfr_gmp
    RREGDEF(mpz2mpfr1);
    RREGDEF(mpz2mpfr1_list);
#endif
#ifdef R_had_R_Outputfile_in_API
#ifndef WIN32
    RREGDEF(print_mpfr);
    RREGDEF(print_mpfr1);
#endif
#endif
    RREGDEF(mpfr2d);
    RREGDEF(mpfr2i);
    RREGDEF(mpfr2str);
    RREGDEF(str2mpfr1_list);

    RREGDEF(Rmpfr_minus);
    RREGDEF(Rmpfr_abs);
    RREGDEF(Math_mpfr);
    RREGDEF(Arith_mpfr);
    RREGDEF(Arith_mpfr_i);
    RREGDEF(Arith_i_mpfr);
    RREGDEF(Arith_mpfr_d);
    RREGDEF(Arith_d_mpfr);
    RREGDEF(Compare_mpfr);
    RREGDEF(Compare_mpfr_i);
    RREGDEF(Compare_mpfr_d);
    RREGDEF(Summary_mpfr);
    RREGDEF(R_mpfr_sumprod);

    RREGDEF(R_mpfr_set_debug);
    RREGDEF(R_mpfr_set_default_prec);
    RREGDEF(R_mpfr_get_default_prec);
    RREGDEF(R_mpfr_get_version);
    RREGDEF(R_mpfr_get_GMP_numb_bits);

    RREGDEF(const_asMpfr);

    RREGDEF(R_mpfr_is_finite);	 RREGDEF(R_mpfr_is_finite_A);
    RREGDEF(R_mpfr_is_infinite); RREGDEF(R_mpfr_is_infinite_A);
    RREGDEF(R_mpfr_is_integer);	 RREGDEF(R_mpfr_is_integer_A);
    RREGDEF(R_mpfr_is_na);	 RREGDEF(R_mpfr_is_na_A);
    RREGDEF(R_mpfr_is_zero);	 RREGDEF(R_mpfr_is_zero_A);

    RREGDEF(R_mpfr_jn);
    RREGDEF(R_mpfr_yn);
    RREGDEF(R_mpfr_atan2);
    RREGDEF(R_mpfr_hypot);
    RREGDEF(R_mpfr_igamma);
    RREGDEF(R_mpfr_beta);
    RREGDEF(R_mpfr_lbeta);
    RREGDEF(R_mpfr_fac);
    RREGDEF(R_mpfr_choose);
    RREGDEF(R_mpfr_poch);
    RREGDEF(R_mpfr_round);


/* Sync this with declarations in ./Syms.h : */
    Rmpfr_precSym = install("prec");
    Rmpfr_signSym = install("sign");
    Rmpfr_expSym = install("exp");
    Rmpfr_d_Sym = install("d");
    Rmpfr_Data_Sym = install(".Data");
    Rmpfr_Dim_Sym = install("Dim");
    Rmpfr_Dimnames_Sym = install("Dimnames");

/* not suppressable, hence moved to suppressable R startup code:
    Rprintf("Loading C code of R package 'Rmpfr': GMP using %d bits per limb\n",
	    GMP_NUMB_BITS);
*/
/*     if(GMP_NUMB_BITS != 32) */
/* 	error("The Rmpfr package currently needs 32-bit limbs"); */
}

/* void R_unload_Rmpfr(DllInfo *dll) */
/* { */

/* } */
