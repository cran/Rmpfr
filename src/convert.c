/*
 * MPFR - Multiple Precision Floating-Point Reliable Library
 * ----   -        -         -              -
 */
#include <Rmath.h>
/* for imax2() */

#include "Rmpfr_utils.h"
extern
#include "Syms.h"

/*------------------------------------------------------------------------*/

/* NB:  int nr_limbs = R_mpfr_nr_limbs(r)  [ in MPFR_as_R() ]  or
 *                   = N_LIMBS(i_prec)     [ in d2mpfr1_()  ]
 *      R_mpfr_exp_size  is sufficient also for enlarged exponent range, as that is still < 2^62
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
    SEXP _V_ = PROTECT(R_do_new_object(PROTECT(R_do_MAKE_CLASS("mpfr1"))));	\
    SEXP prec_R = PROTECT(ALLOC_SLOT(_V_, Rmpfr_precSym, INTSXP, 1));	\
    SEXP sign_R = PROTECT(ALLOC_SLOT(_V_, Rmpfr_signSym, INTSXP, 1));	\
    SEXP exp_R  = PROTECT(ALLOC_SLOT(_V_, Rmpfr_expSym,  INTSXP, R_mpfr_exp_size)); \
    SEXP d_R    = PROTECT(ALLOC_SLOT(_V_, Rmpfr_d_Sym,   INTSXP, _d_length_)); \
    /* the integer vector which makes up the mantissa: */		\
    unsigned int *dd = (unsigned int *) INTEGER(d_R),			\
	*ex = (unsigned int *) INTEGER(exp_R) /* the one for the exponent */

/*------------------------*/

/* Convert integer 'u' of unsigned type 'utype' to corresponding signed
 * type 'stype' without relying on implementation-defined behaviour when
 * 'u' exceeds the maximum of 'stype'; see C99 6.3.1.3
 */
#define CAST_SIGNED(u, utype, stype) \
    (((u) <= ((utype) -1 >> 1)) ? (stype) u : -(stype) ~(u) - 1)

#if GMP_NUMB_BITS == 32
/*                  ---- easy : a gmp-limb is an int <--> R */

static R_INLINE void R_mpfr_FILL_DVEC(int i, mpfr_t r, unsigned int *dd) {
    mp_limb_t limb = r->_mpfr_d[i];
    dd[i] = (unsigned int) limb;
    R_mpfr_dbg_printf(2, "r..d[i=%d] = 0x%lx\n",
                      i, limb);
}

static R_INLINE void R_mpfr_GET_DVEC(int i, mpfr_t r, unsigned int *dd) {
    mp_limb_t limb = (mp_limb_t) dd[i];
    r->_mpfr_d[i] = limb;
    R_mpfr_dbg_printf(2, "dd[%d] = %10lu -> r..d[i=%d] = 0x%lx\n",
                      i, limb, i, limb);
}

// these work on (r, ex[0]) :
static R_INLINE void R_mpfr_FILL_EXP(mpfr_t r, unsigned int *ex) {
    ex[0] = (unsigned int) r->_mpfr_exp;
}
static R_INLINE void R_mpfr_GET_EXP(mpfr_t r, unsigned int *ex,
                                    unsigned int ex1) {
    r->_mpfr_exp = (mpft_exp_t) CAST_SIGNED(ex[0], unsigned int, int);
}


/*------------------------*/
#elif GMP_NUMB_BITS == 64
/*                    ---- here a gmp-limb is 64-bit (long long):
 * ==> one limb  <---> 2 R int.s : */

// This is only ok  if( mpfr_regular_p(.) ), i.e. not for {0, NaN, Inf}:
static R_INLINE void R_mpfr_FILL_DVEC(int i, mpfr_t r, unsigned int *dd) {
    mp_limb_t limb = r->_mpfr_d[i];
    dd[2*i  ] = (unsigned int) (limb & 0x00000000FFFFFFFFu);
    dd[2*i+1] = (unsigned int) (limb >> 32);
    R_mpfr_dbg_printf(2, "r..d[i=%d] = 0x%llx\n",
                      i, (unsigned long long) limb);
}

static R_INLINE void R_mpfr_GET_DVEC(int i, mpfr_t r, unsigned int *dd) {
    mp_limb_t limb = ((mp_limb_t) dd[2*i+1] << 32) | ((mp_limb_t) dd[2*i] & 0x00000000FFFFFFFFu);
    r->_mpfr_d[i] = limb;
    R_mpfr_dbg_printf(2, "dd[%d:%d] = (%10lu,%10lu) -> r..d[i=%d] = 0x%llx\n",
                      2*i, 2*i+1, dd[2*i], dd[2*i+1],
                      i, (unsigned long long) limb);
}

// these work on (r, ex[0], {ex[1] or ex1}) :
static R_INLINE void R_mpfr_FILL_EXP(mpfr_t r, unsigned int *ex) {
    mpfr_uexp_t exponent = (mpfr_uexp_t) r->_mpfr_exp;
    ex[0] = (unsigned int) (exponent & 0x00000000FFFFFFFFu);
    ex[1] = (unsigned int) (exponent >> 32);
    R_mpfr_dbg_printf(2, "_exp = 0x%llx\n",
                      (unsigned long long) exponent);
}

static R_INLINE void R_mpfr_GET_EXP(mpfr_t r, unsigned int *ex,
                                    unsigned int ex1) {
    mpfr_uexp_t exponent = ((mpfr_uexp_t) ex1 << 32) | ((mpfr_uexp_t) ex[0] & 0x00000000FFFFFFFFu);
    r->_mpfr_exp = CAST_SIGNED(exponent, mpfr_uexp_t, mpfr_exp_t);
    R_mpfr_dbg_printf(2, "ex[0:1] = (%10lu,%10lu) -> _exp = 0x%llx\n",
                      ex[0], ex1, (unsigned long long) exponent);
}

/*------------------------*/
#else
# error "will not happen"
#endif



static R_INLINE void
R_mpfr_MPFR_2R_fill(mpfr_t r, unsigned int *ex, int nr_limbs, int regular_p,
		    // Fill (i.e., modify) these :
		    unsigned int *dd, /* = INTEGER(d_R) , the vector which makes up the mantissa */
		    SEXP prec_R, SEXP sign_R)
{
    /* now fill the slots of val */
    INTEGER(prec_R)[0] = (int)r->_mpfr_prec;
    INTEGER(sign_R)[0] = (int)r->_mpfr_sign;
    R_mpfr_FILL_EXP(r, ex);
    if(regular_p) {
	/* the full *vector* of limbs : */
	for(int i=0; i < nr_limbs; i++) {
            R_mpfr_FILL_DVEC(i, r, dd);
	}
    }
}


/* Return an R "mpfr1" object corresponding to mpfr input: */
SEXP MPFR_as_R(mpfr_t r) {

    int nr_limbs = R_mpfr_nr_limbs(r),
	regular_p = mpfr_regular_p(r);

    R_mpfr_MPFR_2R_init(val, (regular_p ? R_mpfr_nr_ints : 0));

    R_mpfr_MPFR_2R_fill(r, ex, nr_limbs, regular_p,
			// Fill these :
			dd, /* = INTEGER(d_R) , the vector which makes up the mantissa */
			prec_R, sign_R);

    UNPROTECT(6);
    return val;
}

SEXP d2mpfr1_(double x, int i_prec, mpfr_rnd_t rnd)
{
    R_mpfr_check_prec(i_prec);

    mpfr_t r;
    mpfr_init2 (r, (mpfr_prec_t)i_prec);
    mpfr_set_d (r, x, rnd);

    int nr_limbs  = N_LIMBS(i_prec),
	regular_p = mpfr_regular_p(r);
    R_mpfr_MPFR_2R_init(val, (regular_p ? R_mpfr_nr_ints : 0));

    R_mpfr_MPFR_2R_fill(r, ex, N_LIMBS(i_prec), regular_p,
			// Fill these :
			dd, /* = INTEGER(d_R) , the vector which makes up the mantissa */
			prec_R, sign_R);

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
	    /* FIXME: become more efficient by doing R_mpfr_MPFR_2R_init() only once*/
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

   *BUT* we cannot easily do the [R package gmp C++ code]-part of  SEXP -> mpz !

   MM - FIXME: still do it .. should not be so hard to "guess"

... horrifically,
... for now, R's  ..bigq2mpfr(), ..bigz2mpfr()  "simply" go via character strings
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
	/* FIXME: become more efficient by doing R_mpfr_MPFR_2R_init() only once*/
	SET_VECTOR_ELT(val, i, MPFR_as_R(r_i));
    }
    mpfr_clear (r_i);
    mpfr_free_cache();
    UNPROTECT(nprot);
    return val;
}

#undef R_mpfr_MPFR_2R_init


#ifdef _not_used_
/* This does *not* work: gives *empty* .Data slot [bug in R_do_new_object()? ] */
SEXP d2mpfr(SEXP x, SEXP prec)
{
    int i_prec = asInteger(prec),
	nx = LENGTH(x), np = LENGTH(prec),
	n = (nx == 0 || np == 0) ? 0 : imax2(nx, np),
	nprot = 1;
    SEXP val = PROTECT(R_do_new_object(R_do_MAKE_CLASS("mpfr"))),
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
 * From an R  "mpfr1" object `x`,  create mpfr `r` (with correct prec): */
void R_asMPFR(SEXP x, mpfr_ptr r)
{
    SEXP prec_R = R_do_slot(x, Rmpfr_precSym);
    // SEXP sign_R = R_do_slot(x, Rmpfr_signSym);// only used once
    SEXP exp_R  = R_do_slot(x, Rmpfr_expSym);
    SEXP d_R    = R_do_slot(x, Rmpfr_d_Sym);

    int x_prec = INTEGER(prec_R)[0],
	nr_limbs = N_LIMBS(x_prec), i;
    Rboolean regular_x = length(d_R) > 0;
    /* the integer vector which makes up the mantissa: */
    unsigned int *dd = (unsigned int *) INTEGER(d_R),
	*ex = (unsigned int *) INTEGER(exp_R), ex1; /* the one for the exponent */

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
    r->_mpfr_sign = (mpfr_sign_t) INTEGER(R_do_slot(x, Rmpfr_signSym))[0];
    R_mpfr_GET_EXP(r, ex, ex1);
    if(regular_x)
	/* the full *vector* of limbs : */
	for(i=0; i < nr_limbs; i++) {
	    R_mpfr_GET_DVEC(i, r, dd);
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
    SEXP D = R_do_slot(x, Rmpfr_Data_Sym);/* an R list() of length n */
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
    mpfr_init(R_i); /* with default precision; set prec in R_asMPFR() */

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
    mpfr_init(R_i); /* with default precision; set prec in R_asMPFR() */

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


/* Get "format info" from  R "mpfr" object -- into list with (exp, finite, is.0),
 * a subset of mpfr2str() [below] : ---> see also R_mpfr_exp() in ./utils.c
 *
 */
SEXP R_mpfr_formatinfo(SEXP x) {
    int n = length(x);
    static const char *ans_nms[] = {"exp", "finite", "is.0", ""};
    SEXP val = PROTECT(mkNamed(VECSXP, ans_nms)), exp, fini, zero;
    int erange_is_int = mpfr_erange_int_p();
    SEXPTYPE exp_SXP = (erange_is_int ? INTSXP : REALSXP);
    SET_VECTOR_ELT(val, 0, exp = PROTECT(allocVector(exp_SXP,n)));
    SET_VECTOR_ELT(val, 1, fini= PROTECT(allocVector(LGLSXP, n)));
    SET_VECTOR_ELT(val, 2, zero= PROTECT(allocVector(LGLSXP, n)));
    int *is_fin= LOGICAL(fini),
	*is_0  = LOGICAL(zero);
    mpfr_t R_i;
    mpfr_init(R_i); /* with default precision; set prec in R_asMPFR() */
    if(erange_is_int) {
	int *exp_ = INTEGER(exp);
#define FOR_I_N_ASSIGN(exp_typ)				\
	for(int i=0; i < n; i++) {			\
	    R_asMPFR(VECTOR_ELT(x, i), R_i);		\
	    exp_  [i] = (exp_typ) mpfr_get_exp(R_i);	\
	    is_fin[i] = mpfr_number_p(R_i);		\
	    is_0  [i] = mpfr_zero_p(R_i);		\
	}
	FOR_I_N_ASSIGN(int)

    } else {/* 'exp' needs to use "double" as it may not fit into integer,
	       consistent with R_mpfr_get_erange(), or  R_mpfr_prec_range() : */
	double *exp_ = REAL(exp);
	FOR_I_N_ASSIGN(double)
    }
    mpfr_clear (R_i);
    mpfr_free_cache();
    UNPROTECT(4);
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
 *
 *  Rmpfr:::.mpfr_debug(1)  ==> to add debug output here
 *
 * mpfr2str() is the workhorse for R level  formatMpfr() , called from format(<mpfr>) <- print(<mpfr>)
 */
SEXP mpfr2str(SEXP x, SEXP digits, SEXP maybeFull, SEXP base) {
    int n = length(x), i;
    int B = asInteger(base); // = base for output
    int n_dig = isNull(digits) ? 0 : asInteger(digits);
    if(n_dig < 0) error("'digits' must be NULL or a positive integer");
    if(asLogical(maybeFull) == NA_LOGICAL) // cannot happen when called "regularly"
	error("'maybe.full' must be TRUE or FALSE");
    Rboolean maybe_full = asLogical(maybeFull);
    R_mpfr_dbg_printf(1,"mpfr2str(*, digits=%d, maybeF=%s, base=%d):\n",
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
    int N_digits = n_dig_1_problem ? 2 : n_dig;
    static const char *ans_nms[] = {"str", "exp", "finite", "is.0", ""};
    SEXP val = PROTECT(mkNamed(VECSXP, ans_nms)), str, exp, fini, zero;
    // NB: 'exp' may have to be 'double' instead of 'integer', when erange allows large exponents
    int erange_is_int = mpfr_erange_int_p();
    SEXPTYPE exp_SXP = (erange_is_int ? INTSXP : REALSXP);
    SET_VECTOR_ELT(val, 0, str = PROTECT(allocVector(STRSXP, n)));
    SET_VECTOR_ELT(val, 1, exp = PROTECT(allocVector(exp_SXP,n)));
    SET_VECTOR_ELT(val, 2, fini= PROTECT(allocVector(LGLSXP, n)));
    SET_VECTOR_ELT(val, 3, zero= PROTECT(allocVector(LGLSXP, n)));
    // depending on erange_is_int, only need one of  d_exp or i_exp  (but don't see a more elegant way):
    double *d_exp; // = REAL(exp);
    int    *i_exp; // = INTEGER(exp),
    int *is_fin= LOGICAL(fini),
	*is_0  = LOGICAL(zero);
    double p_fact = (B == 2) ? 1. : log(B) / M_LN2;// <==> P / p_fact == P *log(2)/log(B)
    int max_nchar = -1; // := max_i { dig_needed[i] }
    mpfr_t R_i;
    mpfr_init(R_i); /* with default precision; set prec in R_asMPFR() */
    if(erange_is_int) {
	i_exp = INTEGER(exp);
	d_exp = NULL;
    } else {
	i_exp = NULL;
	d_exp = REAL(exp);
    }

    for(i=0; i < n; i++) {
	mpfr_exp_t exp = (mpfr_exp_t) 0;
	mpfr_exp_t *exp_ptr = &exp;
	int nchar_i;
	Rboolean use_nchar = TRUE;

	R_asMPFR(VECTOR_ELT(x, i), R_i);

	int is0 = mpfr_zero_p(R_i);
	int isFin = mpfr_number_p(R_i);
	is_0  [i] = is0;
	is_fin[i] = isFin;

	if(N_digits) {/* use it as desired precision */
	    nchar_i = N_digits;
	    R_mpfr_dbg_printf(1,"N_digits: [i=%d]: ... -> dig.n = %d ", i, nchar_i);
	} else if(!isFin) {
	    nchar_i = 5; // @Inf@  @NaN@
	} else if(is0) {
	    nchar_i = 1 + base_is_2_power;
	} else { /* N_digits = 0 --> string must use "enough" digits */
	    // MPFR doc on mpfr_get_str(): use 'm + 1' where  m = 1+ceil(P * log(2)/log(B))
	    double P = (double)R_i->_mpfr_prec;
	    if(base_is_2_power) P--; // P := P-1  iff B is a power of 2
	    double m1 = 1 + ceil(P / p_fact) + 1;
	    double dchar_i = maybe_full ? // want all digits before "." :
		fmax2(m1, ceil((double)mpfr_get_exp(R_i) / p_fact)) : m1;
	    if(dchar_i > 536870912 /* = 2^29 */) // << somewhat arbitrary but < INT_MAX ~= 2^31-1
		error(_(".mpfr2str(): too large 'dchar_i = %g'; please set 'digits = <number>'"),
		      dchar_i);
	    nchar_i = (int) dchar_i;
	    R_mpfr_dbg_printf(1," [i=%d]: prec=%ld, exp2=%ld -> (nchar_i,dig.n)=(%g,%d) ",
			      i, R_i->_mpfr_prec, mpfr_get_exp(R_i),
			      dchar_i, nchar_i);
	    if(nchar_i <= 1 && base_is_2_power) { // have n_dig_problem:
		R_mpfr_dbg_printf_0(1," base_is_2_power & nchar_i=%d ==> fudge dig_n. := 2");
		nchar_i = 2;
	    }
	    use_nchar = FALSE;
	}

	if (i == 0) { /* first time */
	    max_nchar = nchar_i;
	    ch = (char *) R_alloc(imax2(max_nchar + 2, 7), // 7 : '-@Inf@' (+ \0)n_str,
			  sizeof(char));
	}
	else if(!N_digits && nchar_i > max_nchar) { // enlarge :
	    ch = (char *) S_realloc(ch,
				    imax2( nchar_i  + 2, 7),
				    imax2(max_nchar + 2, 7),
				    sizeof(char));
	    max_nchar = nchar_i;
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
	R_mpfr_dbg_printf_0(1," .. max_nchar=%d\n", max_nchar);

	/* // use nchar_i notably when that is smaller than max_nchar : */
	/* mpfr_get_str(ch, exp_ptr, B, (size_t) nchar_i, R_i, MPFR_RNDN);  */
	/* ---- alternatively,
	 * N = 0 : MPFR finds the number of digits needed : */
	mpfr_get_str(ch, exp_ptr, B, (size_t) (maybe_full || use_nchar) ? nchar_i : 0,
	//==========                                                               ---
		     R_i, MPFR_RNDN);

	SET_STRING_ELT(str, i, mkChar(ch));
	if(erange_is_int)
	    i_exp [i] =    (int) exp_ptr[0];
	else
	    d_exp [i] = (double) exp_ptr[0];
    }

    mpfr_clear (R_i);
    mpfr_free_cache();
    UNPROTECT(5);
    return val;
}

/* R  ldexpMpfr(f, E)
 *
 *  return "mpfr" x = f * 2^E   where 'f' is "mpfr" and  'E' is "integer".
 *  -------
 *  "Problem": here, the exponent is limited to  +/- 2^31-1 ("ok" with default erange;
 *             but the maximal erange is +- 2^61 which corresponds to *two* 32-bit integers
 */
SEXP R_mpfr_ldexp(SEXP f, SEXP E, SEXP rnd_mode) {

/* NB: Allow *recycling* for (f, E)
 * --  using 'mismatch' and the MISMATCH macros */
    mpfr_rnd_t rnd = R_rnd2MP(rnd_mode);
    int nprot = 0;
    if(!isInteger(E)) { PROTECT(E = coerceVector(E, INTSXP)); nprot++; }
    int *ee = INTEGER(E),
	nx = length(f), ny = length(E), // instead of 'nf, nE' for MISMATCH macros
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny),
	mismatch = 0;
    SEXP val = PROTECT(allocVector(VECSXP, n)); nprot++;
    mpfr_t x_i;
    mpfr_init(x_i); /* with default precision; set prec in R_asMPFR() */

    SET_MISMATCH;
    for(int i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(f, i % nx), x_i);
	mpfr_mul_2si(x_i, x_i, (long) ee[i % ny], rnd);
	/*
	  -- Function: int mpfr_mul_2si (mpfr_t ROP, mpfr_t OP1, long int OP2, mpfr_rnd_t RND)
	  Set ROP to OP1 times 2 raised to OP2 rounded in the direction RND.
	  Just increases the exponent by OP2 when ROP and OP1 are identical.
	*/
	SET_VECTOR_ELT(val, i, MPFR_as_R(x_i));
    }
    MISMATCH_WARN;

    mpfr_clear (x_i);
    mpfr_free_cache();
    UNPROTECT(nprot);
    return val;
}

#ifdef _not_yet_
/* For R functionality: from "mpfr" x, return  list(z, E),
 * z = "bigz", E = "integer" (or integer-valued double) such that   x = z * 2^E  exactly
 */
SEXP R_mpfr_get_2exp(SEXP x) {

    /*-- Function: mpfr_exp_t mpfr_get_z_2exp (mpz_t ROP, mpfr_t OP)

     Put the scaled significand of OP (regarded as an integer, with the
     precision of OP) into ROP, and return the exponent EXP (which may
     be outside the current exponent range) such that OP exactly equals
     ROP times 2 raised to the power EXP.  If OP is zero, the minimal
     exponent ‘emin’ is returned.  If OP is NaN or an infinity, the
     _erange_ flag is set, ROP is set to 0, and the the minimal exponent
     ‘emin’ is returned.  The returned exponent may be less than the
     minimal exponent ‘emin’ of MPFR numbers in the current exponent
     range; in case the exponent is not representable in the
     ‘mpfr_exp_t’ type, the _erange_ flag is set and the minimal value
     of the ‘mpfr_exp_t’ type is returned.
    */

    // placeholder -- FIXME !
    return R_Nilvalue;
}
#endif // _not_yet_

// R frexpMpfr(x)  |-> list(r, e)
SEXP R_mpfr_frexp(SEXP x, SEXP rnd_mode) {
    mpfr_rnd_t rnd = R_rnd2MP(rnd_mode);
    // NB: 'exp' may have to be 'double' instead of 'integer', when erange allows large exponents
    int erange_is_int = mpfr_erange_int_p();
    SEXPTYPE exp_SXP = (erange_is_int ? INTSXP : REALSXP);
    int n = length(x);
    static const char *ans_nms[] = {"r", "e", ""};
    SEXP ans = PROTECT(mkNamed(VECSXP, ans_nms)), r, e;
    // r: fractional parts; still "mpfr" numbers:
    SET_VECTOR_ELT(ans, 0, r = PROTECT(duplicate(x)));
    // e: 2-exponents (integer or double see aboe)
    SET_VECTOR_ELT(ans, 1, e = PROTECT(allocVector(exp_SXP, n)));
    int    *ei_;
    double *ed_;
    if(erange_is_int)
	ei_ = INTEGER(e);
    else
	ed_ = REAL(e);

    mpfr_t x_i, y_i;
    mpfr_init(x_i); /* with default precision; set prec in R_asMPFR() */
    mpfr_init(y_i);
    mpfr_exp_t Ex = (mpfr_exp_t)0; // create and initialize "valgrindably"
    // Rboolean int_ok = TRUE;
    for(int i=0; i < n; i++) {
	R_asMPFR(VECTOR_ELT(x, i), x_i);
	mpfr_set_prec(y_i, mpfr_get_prec(x_i));
	int ierr = mpfr_frexp(&Ex, y_i, x_i, rnd);
	/*-- Function: int mpfr_frexp (mpfr_exp_t *EXP, mpfr_t Y, mpfr_t X, mpfr_rnd_t RND)

	  Set EXP (formally, the value pointed to by EXP) and Y such that
	  0.5<=abs(Y)<1 and Y times 2 raised to EXP equals X rounded to the
	  precision of Y, using the given rounding mode.  If X is zero, then
	  Y is set to a zero of the same sign and EXP is set to 0.  If X is
	  NaN or an infinity, then Y is set to the same value and EXP is
	  undefined.
	*/
	if(ierr) {
#define MPFR_CLEAR mpfr_clear(x_i); mpfr_clear(y_i); mpfr_free_cache()
	    MPFR_CLEAR;
	    error("R_mpfr_frexp(): mpfr_frexp(x[%d]) gave error code %d\n", i+1, ierr);
	}
	if(erange_is_int)
	    ei_[i] = (int) Ex;
	else
	    ed_[i] = (double) Ex;

	SET_VECTOR_ELT(r, i, MPFR_as_R(y_i));
    }
    MPFR_CLEAR;
    UNPROTECT(3);
    return ans;
}


// TODO : get *exact* as.bigq.mpfr(), i.e, "mpfr" --> "bigq". (R's "bigq" is C 'mpq')
// ---- inside the GMP library, have
//                   _________
// -- Function: void mpq_set_f (mpq_t ROP, const mpf_t OP)
//                   ^^^^^^^^^
//       Set ROP to the value of OP.
//       There is no rounding, this conversion is **exact**.

/*
  ---not-yet---
#include <mpf2mpfr.h>
*/
/* maybe we only need very little of the above, such as
#define mpf_t mpfr_t
--------------------- ?
*/

