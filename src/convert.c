/*
 * MPFR - Multiple Precision Floating-Point Reliable Library
 * ----   -        -         -              -
 */

#include "Rmpfr_utils.h"
#include "Syms.h"

/*------------------------------------------------------------------------*/

#if GMP_NUMB_BITS == 32
# define R_mpfr_nr_ints nr_limbs
# define R_mpfr_exp_size 1
#elif GMP_NUMB_BITS == 64
# define R_mpfr_nr_ints (2*nr_limbs)
# define R_mpfr_exp_size 2
#else
# error "R <-> C Interface *not* implemented for GMP_NUMB_BITS=" ## GMP_NUMB_BITS
#endif

#define R_mpfr_MPFR_2R_init(_V_)					\
    SEXP _V_ = PROTECT(NEW_OBJECT(MAKE_CLASS("mpfr1")));		\
    SEXP prec_R = ALLOC_SLOT(_V_, Rmpfr_precSym, INTSXP, 1);		\
    SEXP sign_R = ALLOC_SLOT(_V_, Rmpfr_signSym, INTSXP, 1);		\
    SEXP exp_R  = ALLOC_SLOT(_V_, Rmpfr_expSym,  INTSXP, R_mpfr_exp_size); \
    SEXP d_R    = ALLOC_SLOT(_V_, Rmpfr_d_Sym,   INTSXP, R_mpfr_nr_ints); \
    /* the integer vector which makes up the mantissa: */		\
    int *dd = INTEGER(d_R),						\
        *ex = INTEGER(exp_R) /* the one for the exponent */

/*------------------------*/
#if GMP_NUMB_BITS == 32
/*                  ---- easy : a gmp-limb is an int <--> R */
# define R_mpfr_FILL_DVEC(i)					\
    R_mpfr_dbg_printf("r..d[i=%d] = 0x%lx\n",i,r->_mpfr_d[i]);	\
    dd[i] = (int) r->_mpfr_d[i]

# define R_mpfr_GET_DVEC(i)					\
    r->_mpfr_d[i] = (mp_limb_t) dd[i];				\
    R_mpfr_dbg_printf("dd[%d] = %10lu -> r..d[i=%d]= 0x%lx\n",	\
		      i, dd[i], i,r->_mpfr_d[i])

# define R_mpfr_FILL_EXP ex[0] = (int)r->_mpfr_exp
# define R_mpfr_GET_EXP  r->_mpfr_exp = (mp_exp_t) ex[0]



/*------------------------*/
#elif GMP_NUMB_BITS == 64
/*                    ---- here a gmp-limb is 64-bit (long long):
 * ==> one limb  <---> 2 R int.s : */

# define RIGHT_HALF(_LONG_) ((long long)(_LONG_) & 0x00000000FFFFFFFF)
/*					  1  4	 8|  4	 8 */
# define LEFT_SHIFT(_LONG_) (((long long)(_LONG_)) << 32)

# define R_mpfr_FILL_DVEC(i)					\
    R_mpfr_dbg_printf("r..d[i=%d] = 0x%lx\n",i,r->_mpfr_d[i]);	\
    dd[2*i]  = (int) RIGHT_HALF(r->_mpfr_d[i]);			\
    dd[2*i+1]= (int) (r->_mpfr_d[i] >> 32)

# define R_mpfr_GET_DVEC(i)						\
    r->_mpfr_d[i] = (mp_limb_t)(RIGHT_HALF(dd[2*i]) | LEFT_SHIFT(dd[2*i+1]));	\
    R_mpfr_dbg_printf("dd[%d:%d]= (%10lu,%10lu) -> r..d[i=%d]= 0x%lx\n", \
             2*i,2*i+1, dd[2*i],dd[2*i+1], i,r->_mpfr_d[i])

# define R_mpfr_FILL_EXP				\
    R_mpfr_dbg_printf("_exp = 0x%lx\n",r->_mpfr_exp);	\
    ex[0] = (int) RIGHT_HALF(r->_mpfr_exp);		\
    ex[1] = (int) (r->_mpfr_exp >> 32)

# define R_mpfr_GET_EXP							\
    r->_mpfr_exp = (mp_exp_t) (RIGHT_HALF(ex[0]) | LEFT_SHIFT(ex1));	\
    R_mpfr_dbg_printf("ex[0:1]= (%10lu,%10lu) -> _exp = 0x%lx\n",	\
                      ex[0],ex1, r->_mpfr_exp)

/*------------------------*/
#else
# error "will not happen"
#endif




#define R_mpfr_MPFR_2R_fill			\
    /* now fill the slots of val */		\
    INTEGER(prec_R)[0] = (int)r->_mpfr_prec;	\
    INTEGER(sign_R)[0] = (int)r->_mpfr_sign;	\
    R_mpfr_FILL_EXP; 				\
    /* the full *vector* of limbs : */		\
    for(i=0; i < nr_limbs; i++) {		\
	R_mpfr_FILL_DVEC(i);			\
    }

/* Return an R "mpfr1" object corresponding to mpfr input: */
SEXP MPFR_as_R(mpfr_t r) {

    int nr_limbs = R_mpfr_nr_limbs(r), i;

    R_mpfr_MPFR_2R_init(val);

    R_mpfr_MPFR_2R_fill;

    UNPROTECT(1);
    return val;
}

SEXP d2mpfr1_(double x, int i_prec)
{
    mpfr_t r;
    int nr_limbs = N_LIMBS(i_prec), i;

    R_mpfr_MPFR_2R_init(val);

    mpfr_init2 (r, (mpfr_prec_t)i_prec);
    mpfr_set_d (r, x, GMP_RNDD);

    R_mpfr_MPFR_2R_fill;

    /* free space used by the MPFR variables */
    mpfr_clear (r);
    mpfr_free_cache(); /* <- Manual 4.8 "Memory Handling" strongly advises ...*/

    UNPROTECT(1);
    return val;
}/* d2mpfr1_ */

SEXP d2mpfr1(SEXP x, SEXP prec) {
    if(LENGTH(x) != 1)
	error("length(x) (=%d) is not 1", LENGTH(x));
    return d2mpfr1_(asReal(x), asInteger(prec));
}

SEXP d2mpfr1_list(SEXP x, SEXP prec)
{
    int *iprec, n = LENGTH(x), np = LENGTH(prec), i, nprot = 1;
    SEXP val = PROTECT(allocVector(VECSXP, n));
    double *dx;

    if(!isReal(x))       { PROTECT(x    = coerceVector(x,   REALSXP)); nprot++; }
    if(!isInteger(prec)) { PROTECT(prec = coerceVector(prec, INTSXP)); nprot++; }
    dx = REAL(x);
    iprec = INTEGER(prec);
    for(i = 0; i < n; i++) {
	/* FIXME: become more efficient by doing R_..._2R_init() only once*/
	SET_VECTOR_ELT(val, i, d2mpfr1_(dx[i], iprec[i % np]));
    }

    UNPROTECT(nprot);
    return val;
}

/* From the MPFR (2.3.2, 2008) doc :
 -- Function:

 int mpfr_set_str (mpfr_t ROP, const char *S, int BASE, mp_rnd_t RND)

     Set ROP to the value of the whole string S in base BASE, rounded
     in the direction RND.  See the documentation of `mpfr_strtofr' for
     a detailed description of the valid string formats.  This function
     returns 0 if the entire string up to the final null character is a
     valid number in base BASE; otherwise it returns -1, and ROP may
     have changed.
*/

SEXP str2mpfr1_list(SEXP x, SEXP prec, SEXP base)
{
/* NB:  prec is "recycled"  within 'x' */
    int ibase = asInteger(base), *iprec,
	n = LENGTH(x), np = LENGTH(prec), i, nprot = 1;
    SEXP val = PROTECT(allocVector(VECSXP, n));
    mpfr_t r_i;
    mpfr_init(r_i);

    if(!isString(x))     { PROTECT(x    = coerceVector(x,    STRSXP)); nprot++; }
    if(!isInteger(prec)) { PROTECT(prec = coerceVector(prec, INTSXP)); nprot++; }

    iprec = INTEGER(prec);

    for(i = 0; i < n; i++) {
	int ierr;
	mpfr_set_prec(r_i, (mpfr_prec_t) iprec[i % np]);
	ierr = mpfr_set_str(r_i, CHAR(STRING_ELT(x, i)), ibase, GMP_RNDD);
	if(ierr)
	    error("str2mpfr1_list(x, *): x[%d] cannot be made into MPFR",
		  i+1);
	/* FIXME: become more efficient by doing R_..._2R_init() only once*/
	SET_VECTOR_ELT(val, i, MPFR_as_R(r_i));
    }
    mpfr_clear (r_i);
    mpfr_free_cache();
    UNPROTECT(nprot);
    return val;
}

#undef R_mpfr_MPFR_2R_init
#undef R_mpfr_MPFR_2R_fill


/* This does not work: gives *empty* .Data slot [bug in NEW_OBJECT()? ] */
SEXP d2mpfr(SEXP x, SEXP prec)
{
    int i_prec = asInteger(prec), n = LENGTH(x), i, nprot = 1;
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("mpfr"))),
	lis = ALLOC_SLOT(val, Rmpfr_Data_Sym, VECSXP, n);
    double *dx;

    if(!isReal(x)) { PROTECT(x = coerceVector(x, REALSXP)); nprot++; }
    REprintf("d2mpfr(x, prec): length(x) = %d, prec = %d -> length(lis) = %d\n",
	     n, i_prec, LENGTH(lis));
    dx = REAL(x);
    for(i = 0; i < n; i++) {
	SET_VECTOR_ELT(lis, i, duplicate(d2mpfr1_(dx[i], i_prec)));
    }

    UNPROTECT(nprot);
    return val;
}


/* The inverse of  MPFR_as_R() :
 * From an R  "mpfr1" object, (re)build an mpfr one: */
void R_asMPFR(SEXP x, mpfr_ptr r)
{
    SEXP prec_R = GET_SLOT(x, Rmpfr_precSym);
    SEXP sign_R = GET_SLOT(x, Rmpfr_signSym);
    SEXP exp_R  = GET_SLOT(x, Rmpfr_expSym);
    SEXP d_R    = GET_SLOT(x, Rmpfr_d_Sym);

    int x_prec = INTEGER(prec_R)[0],
	nr_limbs = N_LIMBS(x_prec), i;
    int *dd = INTEGER(d_R),/* the vector which makes up the mantissa */
        *ex = INTEGER(exp_R), ex1; /* the one for the exponent */

    if(length(d_R) != R_mpfr_nr_ints)
	error("nr_limbs(x_prec)= nr_limbs(%d)= %d : length(<d>) == %d != R_mpfr_nr_ints == %d",
	      x_prec, nr_limbs, length(d_R), R_mpfr_nr_ints);
    if(length(exp_R) < R_mpfr_exp_size) {
	if(length(exp_R) == 0)
	    error("'exp' slot has length 0");
	/* else: we got a 32-bit one in a 64-bit system */
	ex1 = 0;
    } else ex1 = ex[1];

    mpfr_set_prec(r, (mpfr_prec_t) x_prec);
    r->_mpfr_sign = (mpfr_sign_t) INTEGER(sign_R)[0];
    R_mpfr_GET_EXP;
    /* the full *vector* of limbs : */
    for(i=0; i < nr_limbs; i++) {
	R_mpfr_GET_DVEC(i);
    }
    return;
}

#ifndef WIN32
/* This only works on  "unix-alikes" ... but we don't really need it */
/* for R_Outputfile : */
#include <Rinterface.h>


SEXP print_mpfr1(SEXP x, SEXP digits)
{
    mpfr_t r;
    Rboolean use_x_digits = INTEGER(digits)[0] == NA_INTEGER;

    mpfr_init2(r, R_mpfr_prec(x));
    R_asMPFR(x, r);
/*     Rprintf(" * [dbg] after R_asMPFR() ..\n"); */
    mpfr_out_str (R_Outputfile, 10,
		  use_x_digits ? 0 : asInteger(digits),
		  r, GMP_RNDD);
    /* prints the value of s in base 10, rounded towards -Inf, where the third
       argument 0 means that the number of printed digits is automatically
       chosen from the precision of s; */
    Rprintf("\n");

    mpfr_clear (r);
    mpfr_free_cache(); /* <- Manual 4.8 "Memory Handling" strongly advises ...*/
    return x;
}

SEXP print_mpfr(SEXP x, SEXP digits)
{
    SEXP D = GET_SLOT(x, Rmpfr_Data_Sym);/* an R list() of length n */
    int n = length(D), i;
    mpfr_t r;
    Rboolean use_x_digits = INTEGER(digits)[0] == NA_INTEGER;
/* #if MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0) */
/*     char buf[R_BUFSIZE], *p = buf; */
/* #endif */

    mpfr_init(r); /* with default precision */
    for(i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(D, i), r);
/* #if MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0) */

/* 	Rprintf */
/* #else  /\* requires R_Outputfile  from  R's Interfaces.h  ___Unix-alike only__ *\/ */
	mpfr_out_str (R_Outputfile, 10, use_x_digits ? 0 : asInteger(digits),
		      r, GMP_RNDD);
/* #endif */
	Rprintf("\n");
    }

    mpfr_clear (r);
    mpfr_free_cache(); /* <- Manual 4.8 "Memory Handling" strongly advises ...*/
    return x;
}

#endif
/* ^^^ Unix-alike only */

/* Convert R "mpfr" object (list of "mpfr1")  to R "double" vector : */
SEXP mpfr2d(SEXP x) {
    SEXP D = GET_SLOT(x, Rmpfr_Data_Sym);/* an R list() of length */
    int n = length(D), i;
    SEXP val = PROTECT(allocVector(REALSXP, n));
    double *r = REAL(val);
    mpfr_t R_i;

    mpfr_init(R_i); /* with default precision */

    for(i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(D, i), R_i);
	r[i] = mpfr_get_d(R_i, GMP_RNDD);
    }

    mpfr_clear (R_i);
    mpfr_free_cache();
    UNPROTECT(1);
    return val;
}

/* Convert R "mpfr" object (list of "mpfr1")  to R "integer" vector : */
SEXP mpfr2i(SEXP x) {
    SEXP D = GET_SLOT(x, Rmpfr_Data_Sym);/* an R list() of length */
    int n = length(D), i;
    SEXP val = PROTECT(allocVector(INTSXP, n));
    int *r = INTEGER(val);
    mpfr_t R_i;
    mpfr_init(R_i); /* with default precision */

    for(i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(D, i), R_i);
	if(!mpfr_fits_sint_p(R_i, GMP_RNDD)) {
	    warning("NAs introduced by coercion from \"mpfr\" [%d]", i+1);
	    r[i] = NA_INTEGER;
	}
	else {
	    long lr = mpfr_get_si(R_i, GMP_RNDD);
	    r[i] = (int) lr;
	}
    }
    mpfr_clear (R_i);
    mpfr_free_cache();
    UNPROTECT(1);
    return val;
}

/* Convert R "mpfr" object (list of "mpfr1")  to R "character" vector,
 * using precision 'prec' which can be NA/NULL in which case
 * "full precision" (as long as necessary) is used : */
SEXP mpfr2str(SEXP x, SEXP digits) {
    SEXP D = GET_SLOT(x, Rmpfr_Data_Sym);/* an R list() of length */
    int n = length(D), i;
    int n_dig = isNull(digits) ? 0 : asInteger(digits);
    SEXP val = PROTECT(allocVector(VECSXP, 4)),
	nms, str, exp, fini, zero;
    int *i_exp, *is_fin, *is_0;

    mpfr_t R_i;

    if(n_dig < 0)
	error("'digits' must be NULL or integer >= 0");

    SET_VECTOR_ELT(val, 0, str = allocVector(STRSXP, n));
    SET_VECTOR_ELT(val, 1, exp = allocVector(INTSXP, n));
    SET_VECTOR_ELT(val, 2, fini = allocVector(LGLSXP, n));
    SET_VECTOR_ELT(val, 3, zero = allocVector(LGLSXP, n));
    setAttrib(val, R_NamesSymbol, nms = allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, mkChar("str"));
    SET_STRING_ELT(nms, 1, mkChar("exp"));
    SET_STRING_ELT(nms, 2, mkChar("finite"));
    SET_STRING_ELT(nms, 3, mkChar("is.0"));
    i_exp = INTEGER(exp);
    is_fin= INTEGER(fini);
    is_0  = INTEGER(zero);

    mpfr_init(R_i); /* with default precision */

    for(i=0; i < n; i++) {
	mp_exp_t exp = (mp_exp_t) 0;
	mp_exp_t *exp_ptr = &exp;
	char *ch;

	R_asMPFR(VECTOR_ELT(D, i), R_i);

	/* char * mpfr_get_str (char *STR, mp_exp_t *EXPPTR, int B,
	 *			size_t N, mpfr_t OP, mp_rnd_t RND) */
	ch = mpfr_get_str(NULL, exp_ptr, /* B = */ 10,
			  (size_t) n_dig, R_i, GMP_RNDN);
	SET_STRING_ELT(str, i, mkChar(ch));
	i_exp[i] = exp_ptr[0]; /* "FIXME": coerce to 'int' here; ok for non-0 */
	is_fin[i]= mpfr_number_p(R_i);
	is_0 [i] = mpfr_zero_p(R_i);
	mpfr_free_str(ch);
    }

    mpfr_clear (R_i);
    mpfr_free_cache();
    UNPROTECT(1);
    return val;
}
