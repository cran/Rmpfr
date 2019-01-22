/*
 * MPFR - Multiple Precision Floating-Point Reliable Library
 * ----   -        -         -              -
 */
#include <Rmath.h>
/* for imax2() */

#include "Rmpfr_utils.h"
#include "Syms.h"


/*------------------------------------------------------------------------*/

/* NB:  int nr_limbs = R_mpfr_nr_limbs(r)  [ in MPFR_as_R() ]  or
                     = N_LIMBS(i_prec)     [ in d2mpfr1_()  ]
*/
#if GMP_NUMB_BITS == 32
# define R_mpfr_nr_ints nr_limbs
# define R_mpfr_exp_size 1
#elif GMP_NUMB_BITS == 64
# define R_mpfr_nr_ints (2*nr_limbs)
# define R_mpfr_exp_size 2
#else
# error "R <-> C Interface *not* implemented for GMP_NUMB_BITS=" ## GMP_NUMB_BITS
#endif

// Initialize contents (4 slots) of a "mpfr1" R object
#define R_mpfr_MPFR_2R_init(_V_, _d_length_)				\
    SEXP _V_ = PROTECT(NEW_OBJECT(PROTECT(MAKE_CLASS("mpfr1"))));	\
    SEXP prec_R = PROTECT(ALLOC_SLOT(_V_, Rmpfr_precSym, INTSXP, 1));	\
    SEXP sign_R = PROTECT(ALLOC_SLOT(_V_, Rmpfr_signSym, INTSXP, 1));	\
    SEXP exp_R  = PROTECT(ALLOC_SLOT(_V_, Rmpfr_expSym,  INTSXP, R_mpfr_exp_size)); \
    SEXP d_R    = PROTECT(ALLOC_SLOT(_V_, Rmpfr_d_Sym,   INTSXP, _d_length_)); \
    /* the integer vector which makes up the mantissa: */		\
    int *dd = INTEGER(d_R),						\
	*ex = INTEGER(exp_R) /* the one for the exponent */

/*------------------------*/
#if GMP_NUMB_BITS == 32
/*                  ---- easy : a gmp-limb is an int <--> R */

// This is only ok  if( mpfr_regular_p(.) ), i.e. not for {0, NaN, Inf}:
# define R_mpfr_FILL_DVEC(i)					\
    R_mpfr_dbg_printf(2,"r..d[i=%d] = 0x%lx\n",i,r->_mpfr_d[i]); \
    dd[i] = (int) r->_mpfr_d[i]

# define R_mpfr_GET_DVEC(i)					\
    r->_mpfr_d[i] = (mp_limb_t) dd[i];				\
    R_mpfr_dbg_printf(2,"dd[%d] = %10lu -> r..d[i=%d]= 0x%lx\n", \
		      i, dd[i], i,r->_mpfr_d[i])

# define R_mpfr_FILL_EXP ex[0] = (int)r->_mpfr_exp
# define R_mpfr_GET_EXP  r->_mpfr_exp = (mpfr_exp_t) ex[0]



/*------------------------*/
#elif GMP_NUMB_BITS == 64
/*                    ---- here a gmp-limb is 64-bit (long long):
 * ==> one limb  <---> 2 R int.s : */

# define RIGHT_HALF(_LONG_) ((long long)(_LONG_) & 0x00000000FFFFFFFF)
//                                                   1  4   8|  4   8
# define LEFT_SHIFT(_LONG_) (((unsigned long long)(_LONG_)) << 32)

// This is only ok  if( mpfr_regular_p(.) ), i.e. not for {0, NaN, Inf}:
# define R_mpfr_FILL_DVEC(i)						\
    R_mpfr_dbg_printf(2,"r..d[i=%d] = 0x%lx\n",i,r->_mpfr_d[i]);	\
    dd[2*i]  = (int) RIGHT_HALF(r->_mpfr_d[i]);				\
    dd[2*i+1]= (int) (r->_mpfr_d[i] >> 32)

# define R_mpfr_GET_DVEC(i)						\
    r->_mpfr_d[i] = (mp_limb_t)(RIGHT_HALF(dd[2*i]) | LEFT_SHIFT(dd[2*i+1])); \
    R_mpfr_dbg_printf(2,"dd[%d:%d]= (%10lu,%10lu) -> r..d[i=%d]= 0x%lx\n", \
	     2*i,2*i+1, dd[2*i],dd[2*i+1], i,r->_mpfr_d[i])

# define R_mpfr_FILL_EXP				\
    R_mpfr_dbg_printf(2,"_exp = 0x%lx\n",r->_mpfr_exp);	\
    ex[0] = (int) RIGHT_HALF(r->_mpfr_exp);		\
    ex[1] = (int) (r->_mpfr_exp >> 32)

# define R_mpfr_GET_EXP							\
    r->_mpfr_exp = (mpfr_exp_t) (RIGHT_HALF(ex[0]) | LEFT_SHIFT(ex1));	\
    R_mpfr_dbg_printf(2,"ex[0:1]= (%10lu,%10lu) -> _exp = 0x%lx\n",	\
		      ex[0],ex1, r->_mpfr_exp)

/*------------------------*/
#else
# error "will not happen"
#endif



#define R_mpfr_MPFR_2R_fill			\
    /* now fill the slots of val */		\
    INTEGER(prec_R)[0] = (int)r->_mpfr_prec;	\
    INTEGER(sign_R)[0] = (int)r->_mpfr_sign;	\
    R_mpfr_FILL_EXP;				\
    if(regular_p) {				\
	/* the full *vector* of limbs : */	\
	for(i=0; i < nr_limbs; i++) {		\
            R_mpfr_FILL_DVEC(i);		\
	}					\
    }

/* Return an R "mpfr1" object corresponding to mpfr input: */
SEXP MPFR_as_R(mpfr_t r) {

    int nr_limbs = R_mpfr_nr_limbs(r),
	regular_p = mpfr_regular_p(r), i;

    R_mpfr_MPFR_2R_init(val, (regular_p ? R_mpfr_nr_ints : 0));

    R_mpfr_MPFR_2R_fill;

    UNPROTECT(6);
    return val;
}

SEXP d2mpfr1_(double x, int i_prec, mpfr_rnd_t rnd)
{
    mpfr_t r;
    int nr_limbs = N_LIMBS(i_prec), regular_p, i;

    R_mpfr_check_prec(i_prec);

    mpfr_init2 (r, (mpfr_prec_t)i_prec);
    mpfr_set_d (r, x, rnd);

    regular_p = mpfr_regular_p(r);
    R_mpfr_MPFR_2R_init(val, (regular_p ? R_mpfr_nr_ints : 0));

    R_mpfr_MPFR_2R_fill;

    /* free space used by the MPFR variables */
    mpfr_clear (r);
    mpfr_free_cache(); /* <- Manual 4.8 "Memory Handling" strongly advises ...*/

    UNPROTECT(6);
    return val;
}/* d2mpfr1_ */

/**
 * Translate an "R rounding mode" into the correct MPFR one
 *
 * @param rnd_mode: an R character (string with nchar() == 1).
 *
 * @return one of the (currently 4) different  MPFR_RND[DNUZ] modes.
 */
mpfr_rnd_t R_rnd2MP(SEXP rnd_mode) {
    const char* r_ch = CHAR(asChar(rnd_mode));
    switch(r_ch[0]) {
    case 'D': return MPFR_RNDD;
    case 'N': return MPFR_RNDN;
    case 'U': return MPFR_RNDU;
    case 'Z': return MPFR_RNDZ;
    case 'A': return MPFR_RNDA; // since MPFR 3.0.0
    default:
	error(_("illegal rounding mode '%s'; must be one of {'D','N','U','Z','A'}"),
	      r_ch);
	/* Wall: */ return MPFR_RNDN;
    }
}

SEXP d2mpfr1(SEXP x, SEXP prec, SEXP rnd_mode) {
    if(LENGTH(x) != 1)
	error("length(x) (=%d) is not 1", LENGTH(x));
    return d2mpfr1_(asReal(x), asInteger(prec), R_rnd2MP(rnd_mode));
}

SEXP d2mpfr1_list(SEXP x, SEXP prec, SEXP rnd_mode)
{
    int nx = LENGTH(x), np = LENGTH(prec),
	n = (nx == 0 || np == 0) ? 0 : imax2(nx, np),
	nprot = 1;
    SEXP val = PROTECT(allocVector(VECSXP, n));
    if(nx > 0) {
	mpfr_rnd_t rnd = R_rnd2MP(rnd_mode);
	if(!isReal(x))       { PROTECT(x    = coerceVector(x,   REALSXP)); nprot++; }
	if(!isInteger(prec)) { PROTECT(prec = coerceVector(prec, INTSXP)); nprot++; }
	double *dx = REAL(x);
	int *iprec = INTEGER(prec);
	for(int i = 0; i < n; i++) {
	    /* FIXME: become more efficient by doing R_..._2R_init() only once*/
	    SET_VECTOR_ELT(val, i, d2mpfr1_(dx[i % nx], iprec[i % np], rnd));
	}
    }

    UNPROTECT(nprot);
    return val;
}

/*
  -- Function: int mpfr_set_z (mpfr_t ROP, mpz_t OP, mpfr_rnd_t RND)
  -- Function: int mpfr_set_q (mpfr_t ROP, mpq_t OP, mpfr_rnd_t RND)

--> would want functions
	SEXP mpz2mpfr1_(mpz_t x, int i_prec, mpfr_rnd_t rnd);
	SEXP mpz2mpfr1 (SEXP x, SEXP prec, SEXP rnd_mode);
	SEXP mpz2mpfr1_list(SEXP x, SEXP prec, SEXP rnd_mode);
    {and the same for 'q' instead of 'z'}

   completely parallel to the d2mpfr*() functions above

   *BUT* we cannot easily do the [R package gmp C++ code]-part of
   SEXP -> mpz !

   MM: still do it .. should not be so hard to "guess"
*/



/* From the MPFR (2.3.2, 2008) doc :
 -- Function:

 int mpfr_set_str (mpfr_t ROP, const char *S, int BASE, mpfr_rnd_t RND)

     Set ROP to the value of the whole string S in base BASE, rounded
     in the direction RND.  See the documentation of `mpfr_strtofr' for
     a detailed description of the valid string formats.  This function
     returns 0 if the entire string up to the final null character is a
     valid number in base BASE; otherwise it returns -1, and ROP may
     have changed.
*/
SEXP str2mpfr1_list(SEXP x, SEXP prec, SEXP base, SEXP rnd_mode)
{
/* NB: Both x and prec are "recycled" to the longer one if needed */
    int ibase = asInteger(base), *iprec,
	nx = LENGTH(x), np = LENGTH(prec),
	n = (nx == 0 || np == 0) ? 0 : imax2(nx, np),
	nprot = 1;
    SEXP val = PROTECT(allocVector(VECSXP, n));
    mpfr_rnd_t rnd = R_rnd2MP(rnd_mode);
    mpfr_t r_i;
    mpfr_init(r_i);

    if(!isString(x))     { PROTECT(x    = coerceVector(x,    STRSXP)); nprot++; }
    if(!isInteger(prec)) { PROTECT(prec = coerceVector(prec, INTSXP)); nprot++; }
    iprec = INTEGER(prec);

    for(int i = 0; i < n; i++) {
	int prec_i = iprec[i % np];
	R_mpfr_check_prec(prec_i);
	mpfr_set_prec(r_i, (mpfr_prec_t) prec_i);
	int ierr = mpfr_set_str(r_i, CHAR(STRING_ELT(x, i % nx)), ibase, rnd);
	if(ierr) {
	    if (!strcmp("NA", CHAR(STRING_ELT(x, i % nx))))
		mpfr_set_nan(r_i); // "NA" <=> "NaN" (which *are* treated well, by mpfr_set_str)
	    else
		error("str2mpfr1_list(x, *): x[%d] cannot be made into MPFR",
		      i+1);
	}
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


#ifdef _not_used_
/* This does *not* work: gives *empty* .Data slot [bug in NEW_OBJECT()? ] */
SEXP d2mpfr(SEXP x, SEXP prec)
{
    int i_prec = asInteger(prec),
	nx = LENGTH(x), np = LENGTH(prec),
	n = (nx == 0 || np == 0) ? 0 : imax2(nx, np),
	nprot = 1;
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("mpfr"))),
	lis = ALLOC_SLOT(val, Rmpfr_Data_Sym, VECSXP, n);
    double *dx;

    if(!isReal(x)) { PROTECT(x = coerceVector(x, REALSXP)); nprot++; }
    REprintf("d2mpfr(x, prec): length(x) = %d, prec = %d -> length(lis) = %d\n",
	     nx, i_prec, LENGTH(lis));
    dx = REAL(x);
    for(int i = 0; i < n; i++) {
	SET_VECTOR_ELT(lis, i, duplicate(d2mpfr1_(dx [i % nx],
						  i_prec [i % np])));
    }
    UNPROTECT(nprot);
    return val;
}
#endif

/* The inverse of  MPFR_as_R() :
 * From an R  "mpfr1" object, (re)build an mpfr one: */
void R_asMPFR(SEXP x, mpfr_ptr r)
{
    SEXP prec_R = GET_SLOT(x, Rmpfr_precSym);
    // SEXP sign_R = GET_SLOT(x, Rmpfr_signSym);// only used once
    SEXP exp_R  = GET_SLOT(x, Rmpfr_expSym);
    SEXP d_R    = GET_SLOT(x, Rmpfr_d_Sym);

    int x_prec = INTEGER(prec_R)[0],
	nr_limbs = N_LIMBS(x_prec), i;
    Rboolean regular_x = length(d_R) > 0;
    int *dd = INTEGER(d_R),/* the vector which makes up the mantissa */
	*ex = INTEGER(exp_R), ex1; /* the one for the exponent */

    if(regular_x && length(d_R) != R_mpfr_nr_ints)
	error("nr_limbs(x_prec)= nr_limbs(%d)= %d : length(<d>) == %d != R_mpfr_nr_ints == %d",
	      x_prec, nr_limbs, length(d_R), R_mpfr_nr_ints);
    if(length(exp_R) < R_mpfr_exp_size) {
	if(length(exp_R) == 0)
	    error("'exp' slot has length 0");
	/* else: we got a 32-bit one in a 64-bit system */
	ex1 = 0;
    } else ex1 = ex[1];

    mpfr_set_prec(r, (mpfr_prec_t) x_prec);
    r->_mpfr_sign = (mpfr_sign_t) INTEGER(GET_SLOT(x, Rmpfr_signSym))[0];
    R_mpfr_GET_EXP;
    if(regular_x)
	/* the full *vector* of limbs : */
	for(i=0; i < nr_limbs; i++) {
	    R_mpfr_GET_DVEC(i);
	}
    return;
}

#ifdef R_had_R_Outputfile_in_API
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
		  r, MPFR_RNDD);
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

/*	Rprintf */
/* #else  /\* requires R_Outputfile  from  R's Interfaces.h  ___Unix-alike only__ *\/ */
	mpfr_out_str (R_Outputfile, 10, use_x_digits ? 0 : asInteger(digits),
		      r, MPFR_RNDD);
/* #endif */
	Rprintf("\n");
    }

    mpfr_clear (r);
    mpfr_free_cache(); /* <- Manual 4.8 "Memory Handling" strongly advises ...*/
    return x;
}

#endif
/* ^^^ Unix-alike only */
#endif


/* Convert R "mpfr" object (list of "mpfr1")  to R "double" vector : */
SEXP mpfr2d(SEXP x, SEXP rnd_mode) {
    int n = length(x), i;
    SEXP val = PROTECT(allocVector(REALSXP, n));
    double *r = REAL(val);
    mpfr_t R_i;
    mpfr_init(R_i); /* with default precision */

    for(i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(x, i), R_i);
	r[i] = mpfr_get_d(R_i, R_rnd2MP(rnd_mode));
    }

    mpfr_clear (R_i);
    mpfr_free_cache();
    UNPROTECT(1);
    return val;
}

/* Convert R "mpfr" object (list of "mpfr1")  to R "integer" vector : */
SEXP mpfr2i(SEXP x, SEXP rnd_mode) {
    int n = length(x), i;
    SEXP val = PROTECT(allocVector(INTSXP, n));
    int *r = INTEGER(val);
    mpfr_t R_i;
    mpfr_init(R_i); /* with default precision */

    for(i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(x, i), R_i);
	if(!mpfr_fits_sint_p(R_i, R_rnd2MP(rnd_mode))) {
	    warning("NAs introduced by coercion from \"mpfr\" [%d]", i+1);
	    r[i] = NA_INTEGER;
	}
	else {
	    long lr = mpfr_get_si(R_i, R_rnd2MP(rnd_mode));
	    r[i] = (int) lr;
	}
    }
    mpfr_clear (R_i);
    mpfr_free_cache();
    UNPROTECT(1);
    return val;
}

/* Convert R "mpfr" object (list of "mpfr1")  to R "character" vector,
 * using 'digits' (or determinining it):
 *  1) digits = NULL , maybe_full = FALSE (<==> 'scientific = TRUE')
 *     --------------- -------------~~~~~ ==> set digits  <=>  getPrec(x)

 *  2) digits = NULL , maybe_full = TRUE (<=> 'scientific' = NA or FALSE  in the calling formatMpfr()
 *     -------------   ----------------- ==> set digits  <=>  max(getPrec(x), #{"digits left of '.'"})

 *  3) digits = <num>, maybe_full = TRUE (<=> 'scientific' = TRUE  in the calling formatMpfr()
 *     -------------   -----------------==> set digits  <=>  max(digit, getPrec(x), #{"digits left of '.'"}))

 */
SEXP mpfr2str(SEXP x, SEXP digits, SEXP maybeFull, SEXP base) {
    int n = length(x), i;
    int B = asInteger(base); // = base for output
    int n_dig = isNull(digits) ? 0 : asInteger(digits);
    if(n_dig < 0) error("'digits' must be NULL or a positive integer");
    Rboolean maybe_full = asLogical(maybeFull);
    if(maybe_full == NA_LOGICAL) // cannot happen when called "regularly"
	error("'maybe.full' must be TRUE or FALSE");

    R_mpfr_dbg_printf(1,"mpfr2str(*, digits=%d, maybeF=%s, base=%d): ",
		      n_dig, (maybe_full ? "TRUE" : "False"), B);

    /* int dig_n_max = -1; */
    /* SEXP val = PROTECT(allocVector(VECSXP, 4)), */
    /* 	nms, str, exp, fini, zero; */
    /* int *i_exp, *is_fin, *is_0; */
    char *ch = NULL;

    /* N_digits == 1 , for base = 2, 4, 8, 16, 32 (base <= 62 !) gives bad abort from MPFR:
       get_str.c:2306: MPFR assertion failed: m >= 2 || ((((b) & ((b) - 1)) == 0) == 0 && m >= 1)
       ...  Aborted ... (Speicherabzug geschrieben)

      the MPFR doc (see mpfr_get_str below) *says* N >= 2 is required,
      but we have used N = 1 for B = 10 a lot in the past ! */

    Rboolean base_is_2_power = (B == 2 || B == 4 || B == 8 || B == 16 || B == 32);
    Rboolean n_dig_1_problem = (n_dig == 1) && base_is_2_power;
    size_t N_digits = n_dig_1_problem ? 2 : n_dig;
    SEXP
	val = PROTECT(allocVector(VECSXP, 4)),
	nms = PROTECT(allocVector(STRSXP, 4)), str, exp, fini, zero;
    SET_VECTOR_ELT(val, 0, str = PROTECT(allocVector(STRSXP, n))); SET_STRING_ELT(nms, 0, mkChar("str"));
    SET_VECTOR_ELT(val, 1, exp = PROTECT(allocVector(INTSXP, n))); SET_STRING_ELT(nms, 1, mkChar("exp"));
    SET_VECTOR_ELT(val, 2, fini= PROTECT(allocVector(LGLSXP, n))); SET_STRING_ELT(nms, 2, mkChar("finite"));
    SET_VECTOR_ELT(val, 3, zero= PROTECT(allocVector(LGLSXP, n))); SET_STRING_ELT(nms, 3, mkChar("is.0"));
    setAttrib(val, R_NamesSymbol, nms);
    int *i_exp = INTEGER(exp),
	*is_fin= LOGICAL(fini),
	*is_0  = LOGICAL(zero);
    double p_fact = (B == 2) ? 1. : log(B) / M_LN2;
    int dig_n_max = -1; // := max_i { dig_needed[i] }
    mpfr_t R_i;
    mpfr_init(R_i); /* with default precision */

    for(i=0; i < n; i++) {
	mpfr_exp_t exp = (mpfr_exp_t) 0;
	mpfr_exp_t *exp_ptr = &exp;
	int dig_needed;

	R_asMPFR(VECTOR_ELT(x, i), R_i);

#ifdef __Rmpfr_FIRST_TRY_FAILS__
/* Observing memory problems, e.g., see ../tests/00-bug.R.~3~
 * Originally hoped it was solvable via  R_alloc() etc, but it seems the problem is
 * deeper and I currently suspect a problem/bug in MPFR library's  mpfr_get_str(..) */
	ch = mpfr_get_str(NULL, exp_ptr, B,
			  (size_t) N_digits, R_i, MPFR_RNDN);
#else
	if(n_dig) {/* use it as desired precision */
	    dig_needed = N_digits;
	    R_mpfr_dbg_printf(1," [i=%d]: ... -> dig.n = %d ", i, dig_needed);
	} else { /* N_digits = 0 --> string must use "enough" digits */
	    double need_dig =
		ceil(fmax2((double)R_i->_mpfr_prec,
			   // when prec is too small:
			   (double)mpfr_get_exp(R_i)) / p_fact);
	    if(need_dig > 268435456 /* = 2^28 */) // << FIXME, somewhat arbitrary
		error(_(".mpfr2str(): too large 'need_dig'; please set 'digits = <number>'"));
// FIXME: rather set   maybe_full = FALSE  (???)
	    dig_needed = (int) need_dig;
	    R_mpfr_dbg_printf(1," [i=%d]: prec=%ld, exp2=%ld -> (n.dig,dig.n)=(%g,%d) ",
			      i, R_i->_mpfr_prec, mpfr_get_exp(R_i),
			      need_dig, dig_needed);
	    if(dig_needed <= 1 && base_is_2_power) { // have n_dig_problem:
		R_mpfr_dbg_printf(1," [i=%d]: base_is_2_power & dig_needed=%d ==> fudge dig_n. := 2",
				  i, dig_needed);
		dig_needed = 2;
	    }
	}
	if (i == 0) { /* first time */
	    dig_n_max = dig_needed;
	    ch = (char *) R_alloc(imax2(dig_n_max + 2, 7), // 7 : '-@Inf@' (+ \0)n_str,
				  sizeof(char));
	}
	else if(!N_digits && dig_needed > dig_n_max) { // enlarge :
	    ch = (char *) S_realloc(ch,
				    imax2(dig_needed + 2, 7),
				    imax2(dig_n_max  + 2, 7),
				    sizeof(char));
	    dig_n_max = dig_needed;
	}

	/*  char* mpfr_get_str (char *STR, mpfr_exp_t *EXPPTR, int B,
	 *			size_t N, mpfr_t OP, mpfr_rnd_t RND)

	 Convert OP to a string of digits in base B, with rounding in the
	 direction RND, where N is either zero (see below) or the number of
	 significant digits output in the string; in the latter case, N must
	 be greater or equal to 2.  The base may vary from 2 to 62;
	 .........
	 .........  ==> MPFR info manual  "5.4 Conversion Functions"
	*/
	R_mpfr_dbg_printf(1," .. dig_n_max=%d\n", dig_n_max);

	mpfr_get_str(ch, exp_ptr, B,
		     (size_t) dig_n_max, R_i, MPFR_RNDN);
#endif
	SET_STRING_ELT(str, i, mkChar(ch));
	i_exp[i] = (int) exp_ptr[0];
	is_fin[i]= mpfr_number_p(R_i);
	is_0 [i] = mpfr_zero_p(R_i);
#ifdef __Rmpfr_FIRST_TRY_FAILS__
	mpfr_free_str(ch);
#endif
    }

    mpfr_clear (R_i);
    mpfr_free_cache();
    UNPROTECT(6);
    return val;
}
