/*
 * MPFR - Multiple Precision Floating-Point Reliable Library
 * ----   -        -         -              -
 */
#include <Rmath.h>
/* imax2() */

#include "Rmpfr_utils.h"
extern
#include "Syms.h"

//Dbg: #define DEBUG_Rmpfr
#ifdef DEBUG_Rmpfr
/* ONLY for debugging  !! */
# include <R_ext/Utils.h>
//-> void R_CheckUserInterrupt(void);
# ifndef WIN32
#  include <Rinterface.h>
# endif
# define R_PRT(_X_) mpfr_out_str (R_Outputfile, 10, 0, _X_, MPFR_RNDD)
#endif

// Currently not in the API (hence "should be" (?) 'static') :
int my_mpfr_beta (mpfr_t ROP, mpfr_t X, mpfr_t Y, mpfr_rnd_t RND);
int my_mpfr_lbeta(mpfr_t ROP, mpfr_t X, mpfr_t Y, mpfr_rnd_t RND);

int my_mpfr_choose(mpfr_t ROP, long n,    mpfr_t X, mpfr_rnd_t RND);
int my_mpfr_poch  (mpfr_t ROP, long n,    mpfr_t X, mpfr_rnd_t RND);
int my_mpfr_round (mpfr_t ROP, long prec, mpfr_t X, mpfr_rnd_t RND);
/* argument order above must match the one of mpfr_jn() etc .. */

/* MM: for debugging, use
  gcc -I/u/maechler/R/D/r-patched/F19-64-inst/include  -I/usr/local/include -fpic  -g -O3 -pedantic -Wall --std=gnu99 -DDEBUG_Rmpfr -Wcast-align -Wclobbered  -c utils.c -o utils.o
*/
/*------------------------------------------------------------------------*/
int my_mpfr_beta (mpfr_t R, mpfr_t a, mpfr_t b, mpfr_rnd_t RND)
{
    mpfr_prec_t p_a = mpfr_get_prec(a), p_b = mpfr_get_prec(b);
    if(p_a < p_b) p_a = p_b;// p_a := max(p_a, p_b)
    if(mpfr_get_prec(R) < p_a)
	mpfr_prec_round(R, p_a, RND);// so prec(R) = max( prec(a), prec(b) )
    int ans;
    mpfr_t s; mpfr_init2(s, p_a);
#ifdef DEBUG_Rmpfr
    R_CheckUserInterrupt();
    int cc = 0;
#endif

    /* "FIXME": check each 'ans' below, and return when not ok ... */
    ans = mpfr_add(s, a, b, RND);

    if(mpfr_integer_p(s) && mpfr_sgn(s) <= 0) { // (a + b) is integer <= 0
	if(!mpfr_integer_p(a) && !mpfr_integer_p(b)) {
	    // but a,b not integer ==> R =  finite / +-Inf  = 0 :
	    mpfr_set_zero (R, +1);
	    mpfr_clear (s);
	    return ans;
	}// else: sum is integer; at least one {a,b} integer ==> both integer

	int sX = mpfr_sgn(a), sY = mpfr_sgn(b);
	if(sX * sY < 0) { // one negative, one positive integer
	    // ==> special treatment here :
	    if(sY < 0) // ==> sX > 0; swap the two
		mpfr_swap(a, b);
	    // now have --- a < 0 < b <= |a|  integer ------------------
	    /*              ================  and in this case:
	       B(a,b) = (-1)^b  B(1-a-b, b) = (-1)^b B(1-s, b)

		      = (1*2*..*b) / (-s-1)*(-s-2)*...*(-s-b)
	    */
	    /* where in the 2nd form, both numerator and denominator have exactly
	     * b integer factors. This is attractive {numerically & speed wise}
	     * for 'small' b */
#define b_large 100
#ifdef DEBUG_Rmpfr
	    Rprintf(" my_mpfr_beta(<neg int>): s = a+b= "); R_PRT(s);
	    Rprintf("\n   a = "); R_PRT(a);
	    Rprintf("\n   b = "); R_PRT(b); Rprintf("\n");
	    if(cc++ > 999) { mpfr_set_zero (R, +1); mpfr_clear (s); return ans; }
#endif
	    unsigned long b_ = 0;// -Wall
	    Rboolean
		b_fits_ulong = mpfr_fits_ulong_p(b, RND),
		small_b = b_fits_ulong &&  (b_ = mpfr_get_ui(b, RND)) < b_large;
	    if(small_b) {
#ifdef DEBUG_Rmpfr
		Rprintf("   b <= b_large = %d...\n", b_large);
#endif
		//----------------- small b ------------------
		// use GMP big integer arithmetic:
		mpz_t S; mpz_init(S); mpfr_get_z(S, s, RND); // S := s
		mpz_sub_ui (S, S, (unsigned long) 1); // S = s - 1 = (a+b-1)
		/* binomial coefficient choose(N, k) requires k a 'long int';
		 * here, b must fit into a long: */
		mpz_bin_ui (S, S, b_); // S = choose(S, b) = choose(a+b-1, b)
		mpz_mul_ui (S, S, b_); // S = S*b =  b * choose(a+b-1, b)
		// back to mpfr: R = 1 / S  = 1 / (b * choose(a+b-1, b))
		mpfr_set_ui(s, (unsigned long) 1, RND);
		mpfr_div_z(R, s, S, RND);
		mpz_clear(S);
	    }
	    else { // b is "large", use direct B(.,.) formula
#ifdef DEBUG_Rmpfr
		Rprintf("   b > b_large = %d...\n", b_large);
#endif
		// a := (-1)^b :
		// there is no  mpfr_si_pow(a, -1, b, RND);
		int neg; // := 1 ("TRUE") if (-1)^b = -1, i.e. iff  b is odd
		if(b_fits_ulong) { // (i.e. not very large)
		    neg = (b_ % 2); // 1 iff b_ is odd,  0 otherwise
		} else { // really large b; as we know it is integer, can still..
		    // b2 := b / 2
		    mpfr_t b2; mpfr_init2(b2, p_a);
		    mpfr_div_2ui(b2, b, 1, RND);
		    neg = !mpfr_integer_p(b2); // b is odd, if b/2 is *not* integer
#ifdef DEBUG_Rmpfr
		    Rprintf("   really large b; neg = ('b is odd') = %d\n", neg);
#endif
		}
		// s' := 1-s = 1-a-b
		mpfr_ui_sub(s, 1, s, RND);
#ifdef DEBUG_Rmpfr
		Rprintf("  neg = %d\n", neg);
		Rprintf("  s' = 1-a-b = "); R_PRT(s);
		Rprintf("\n  -> calling B(s',b)\n");
#endif
		// R := B(1-a-b, b) = B(s', b)
		if(small_b) {
		    my_mpfr_beta (R, s, b, RND);
		} else {
		    my_mpfr_lbeta (R, s, b, RND);
		    mpfr_exp(R, R, RND); // correct *if* beta() >= 0
		}
#ifdef DEBUG_Rmpfr
		Rprintf("  R' = beta(s',b) = "); R_PRT(R); Rprintf("\n");
#endif
		// Result = (-1)^b  B(1-a-b, b) = +/- s'
		if(neg) mpfr_neg(R, R, RND);
	    }
	    mpfr_clear(s);
	    return ans;
	}
   }

    ans = mpfr_gamma(s, s, RND);  /* s = gamma(a + b) */
#ifdef DEBUG_Rmpfr
    Rprintf("my_mpfr_beta(): s = gamma(a+b)= "); R_PRT(s);
    Rprintf("\n   a = "); R_PRT(a);
    Rprintf("\n   b = "); R_PRT(b);
#endif

    ans = mpfr_gamma(a, a, RND);
    ans = mpfr_gamma(b, b, RND);
    ans = mpfr_mul(b, b, a, RND); /* b' = gamma(a) * gamma(b) */

#ifdef DEBUG_Rmpfr
    Rprintf("\n    G(a) * G(b) = "); R_PRT(b); Rprintf("\n");
#endif

    ans = mpfr_div(R, b, s, RND);
    mpfr_clear (s);
    /* mpfr_free_cache() must be called in the caller !*/
    return ans;
}

int my_mpfr_lbeta(mpfr_t R, mpfr_t a, mpfr_t b, mpfr_rnd_t RND)
{
    mpfr_prec_t p_a = mpfr_get_prec(a), p_b = mpfr_get_prec(b);
    if(p_a < p_b) p_a = p_b;// p_a := max(p_a, p_b)
    if(mpfr_get_prec(R) < p_a)
	mpfr_prec_round(R, p_a, RND);// so prec(R) = max( prec(a), prec(b) )
    int ans;
    mpfr_t s;
    mpfr_init2(s, p_a);

    /* "FIXME": check each 'ans' below, and return when not ok ... */
    ans = mpfr_add(s, a, b, RND);

    if(mpfr_integer_p(s) && mpfr_sgn(s) <= 0) { // (a + b) is integer <= 0
	if(!mpfr_integer_p(a) && !mpfr_integer_p(b)) {
	    // but a,b not integer ==> R = ln(finite / +-Inf) = ln(0) = -Inf :
	    mpfr_set_inf (R, -1);
	    mpfr_clear (s);
	    return ans;
	}// else: sum is integer; at least one integer ==> both integer

	int sX = mpfr_sgn(a), sY = mpfr_sgn(b);
	if(sX * sY < 0) { // one negative, one positive integer
	    // ==> special treatment here :
	    if(sY < 0) // ==> sX > 0; swap the two
		mpfr_swap(a, b);
	    /* now have --- a < 0 < b <= |a|  integer ------------------
	     *              ================
	     * --> see my_mpfr_beta() above */
	    unsigned long b_ = 0;// -Wall
	    Rboolean
		b_fits_ulong = mpfr_fits_ulong_p(b, RND),
		small_b = b_fits_ulong &&  (b_ = mpfr_get_ui(b, RND)) < b_large;
	    if(small_b) {
		//----------------- small b ------------------
		// use GMP big integer arithmetic:
		mpz_t S; mpz_init(S); mpfr_get_z(S, s, RND); // S := s
		mpz_sub_ui (S, S, (unsigned long) 1); // S = s - 1 = (a+b-1)
		/* binomial coefficient choose(N, k) requires k a 'long int';
		 * here, b must fit into a long: */
		mpz_bin_ui (S, S, b_); // S = choose(S, b) = choose(a+b-1, b)
		mpz_mul_ui (S, S, b_); // S = S*b =  b * choose(a+b-1, b)

		// back to mpfr: R = log(|1 / S|) =  - log(|S|)
		mpz_abs(S, S);
		mpfr_set_z(s, S, RND); // <mpfr> s :=  |S|
		mpfr_log(R, s, RND);   // R := log(s) = log(|S|)
		mpfr_neg(R, R, RND);   // R = -R = -log(|S|)
		mpz_clear(S);
	    }
	    else { // b is "large", use direct B(.,.) formula
		// a := (-1)^b -- not needed here, neither 'neg': want log( |.| )
		// s' := 1-s = 1-a-b
		mpfr_ui_sub(s, 1, s, RND);
		// R := log(|B(1-a-b, b)|) = log(|B(s', b)|)
		my_mpfr_lbeta (R, s, b, RND);
	    }
	    mpfr_clear(s);
	    return ans;
	}
    }

    ans = mpfr_lngamma(s, s, RND); // s = lngamma(a + b)
    ans = mpfr_lngamma(a, a, RND);
    ans = mpfr_lngamma(b, b, RND);
    ans = mpfr_add (b, b, a, RND); // b' = lngamma(a) + lngamma(b)
    ans = mpfr_sub (R, b, s, RND);

    mpfr_clear (s);
    return ans;
}

/** Binomial Coefficient --
 * all initialization and cleanup is called in the caller
 * @result R = choose(X, n)
 */
int my_mpfr_choose (mpfr_t R, long n, mpfr_t X, mpfr_rnd_t RND)
{
    int ans;
    long i;
    mpfr_t r, x;
    mpfr_prec_t p_X = mpfr_get_prec(X);

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
int my_mpfr_poch (mpfr_t R, long n, mpfr_t X, mpfr_rnd_t RND)
{
    int ans;
    long i;
    mpfr_t r, x;
    mpfr_prec_t p_X = mpfr_get_prec(X);

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
int my_mpfr_round (mpfr_t R, long prec, mpfr_t X, mpfr_rnd_t RND)
{
    int ans;
    if(prec < MPFR_PREC_MIN)
	error("prec = %d < %d  is too small", prec, MPFR_PREC_MIN);
    if(prec > MPFR_PREC_MAX)
	error("prec = %d > %d  is too large", prec, MPFR_PREC_MAX);
    mpfr_set(R, X, RND);
    ans = mpfr_prec_round(R, (mpfr_prec_t) prec, RND);
    return ans;
}

/*------------------------------------------------------------------------*/

SEXP R_mpfr_get_version(void) {
    return mkString(mpfr_get_version());
}

SEXP R_mpfr_get_GMP_numb_bits(void) {// for diagnosing
    return ScalarInteger((int)GMP_NUMB_BITS);
}

/* Set or get the C-global debugging level --
 * currently only used in R_mpfr_dbg_printf() --> ./Rmpfr_utils.h
 *
 * Called from R  .mpfr_debug(i = NA)
*/
SEXP R_mpfr_set_debug(SEXP I)
{
    if(LENGTH(I) < 1 || INTEGER(I)[0] == NA_INTEGER)
	return ScalarInteger(R_mpfr_debug_);
    /* else : */
    int prev = R_mpfr_debug_;
    R_mpfr_debug_ = asInteger(I);
    return ScalarInteger(prev);
}

SEXP R_mpfr_get_default_prec(void) {
    return ScalarInteger((int) mpfr_get_default_prec());
}

SEXP R_mpfr_set_default_prec(SEXP prec) {
    // return the previous value
    int prev = (int) mpfr_get_default_prec();
    mpfr_set_default_prec((mpfr_prec_t) asInteger(prec));
    return ScalarInteger(prev);
}

// is MPFR's exponent range 'erange' representable as R's  (32 bit) integer  [INT_MIN not allowed] :
int mpfr_erange_int_p(void) {
    mpfr_exp_t r = mpfr_get_emin();
    int i_ok = (INT_MIN < r && r <= INT_MAX);
    if(i_ok) {
	r = mpfr_get_emax();
	i_ok = (INT_MIN < r && r <= INT_MAX);
    }
    return i_ok;
}
/** R's .mpfr_erange_is_int() - workhorse
 */
SEXP R_mpfr_erange_int_p(void) {
    return ScalarLogical(mpfr_erange_int_p());
}

/* MUST be sync'ed with  ../R/mpfr.R
 *                       ~~~~~~~~~~~ and its  .mpfr_erange_kinds
 */
typedef enum { E_min = 1, E_max,
	       min_emin, max_emin, min_emax, max_emax } erange_kind;

// Called from R's  .mpfr.erange(), now allows 'kind' to be a vector
SEXP R_mpfr_get_erange(SEXP kind_) {
    int k = LENGTH(kind_), nprot = 0;
    erange_kind *kind;
    if(TYPEOF(kind_) != INTSXP) {
	SEXP kk = PROTECT(coerceVector(kind_, INTSXP)); nprot++;
	kind = (erange_kind *) INTEGER(kk);
    } else {
	kind = (erange_kind *) INTEGER(kind_);
    }

    mpfr_exp_t *r = (mpfr_exp_t *) R_alloc(k, sizeof(mpfr_exp_t));
    Rboolean int_ok = TRUE;
    for(int j = 0; j < k; j++) {
	switch(kind[j]) { // keep the 'case' list in sync with 'erange_kind' enum above:
	case E_min:    r[j] = mpfr_get_emin();     if(int_ok && (r[j] <= INT_MIN || r[j] > INT_MAX)) int_ok=FALSE; break;
	case E_max:    r[j] = mpfr_get_emax();     if(int_ok && (r[j] <= INT_MIN || r[j] > INT_MAX)) int_ok=FALSE; break;
	case min_emin: r[j] = mpfr_get_emin_min(); if(int_ok) int_ok=FALSE; break;
	case max_emin: r[j] = mpfr_get_emin_max(); if(int_ok) int_ok=FALSE; break;
	case min_emax: r[j] = mpfr_get_emax_min(); if(int_ok) int_ok=FALSE; break;
	case max_emax: r[j] = mpfr_get_emax_max(); if(int_ok) int_ok=FALSE; break;
	default:
	    error("invalid kind[j(=%d)] (code = %d) in R_mpfr_get_erange()", j, kind);
	}
	R_mpfr_dbg_printf(1,"R_mpfr_get_erange(%d): %ld\n", kind[j], (long)r[j]);
    }
    SEXP ans;
    // int_ok: only now know if we can return integer or need double
    if(int_ok) {
	int* R = INTEGER(ans = allocVector(INTSXP, k));
	for(int j = 0; j < k; j++) R[j] = (int) r[j];
    } else {
	double* R = REAL(ans = allocVector(REALSXP, k));
	for(int j = 0; j < k; j++) R[j] = (double) r[j];
    }
    if(nprot) UNPROTECT(nprot);
    return ans;
}

SEXP R_mpfr_set_erange(SEXP kind_, SEXP val) {
    erange_kind kind = asInteger(kind_);
    mpfr_exp_t exp_val;
    if(isInteger(val))
	exp_val = asInteger(val);// assume this is always valid to set

    else { // we allow larger values from the R side
	PROTECT(val = coerceVector(val, REALSXP));
	exp_val = (mpfr_exp_t) asReal(val);
	UNPROTECT(1);
    }

    int i_err;
    switch(kind) {
    case E_min: i_err = mpfr_set_emin(exp_val); break;
    case E_max: i_err = mpfr_set_emax(exp_val); break;
    default:
	error("invalid kind (code = %d) in R_mpfr_set_erange()", kind);
    }
    if(i_err) warning("e%s exponent could not be set to %ld (code %d)",
		      (kind == E_min) ? "min" : "max", (long)exp_val, i_err);
    return ScalarInteger(i_err);
}

SEXP R_mpfr_prec_range(SEXP ind) {
    long r = (long) (
	(INTEGER(ind)[0] == 1)
	? MPFR_PREC_MIN
	: MPFR_PREC_MAX);
    R_mpfr_dbg_printf(1,"R_mpfr_prec_range(): %ld\n", r);
    // in 64 bit, int << long, so go "2nd best":
    return ScalarReal((double)r);
}

/** Get the  'base 2  exp slot' -- also in extended erange where it does not fit into integer
 *  Directly called from  R's  .mpfr2exp(x)
 */
SEXP R_mpfr_2exp(SEXP x) {
    int n = length(x);
    mpfr_t R_i; mpfr_init(R_i);
    SEXP val;
    if(mpfr_erange_int_p()) { // integer is ok: 'exp' values won't be too large
	val = PROTECT(allocVector(INTSXP,n));
	int *exp = INTEGER(val);
	for(int i=0; i < n; i++) {
	    R_asMPFR(VECTOR_ELT(x, i), R_i);
	    exp[i] = (int) mpfr_get_exp(R_i);
	}
    } else {
	val = PROTECT(allocVector(REALSXP,n));
	double *exp = REAL(val);
	for(int i=0; i < n; i++) {
	    R_asMPFR(VECTOR_ELT(x, i), R_i);
	    exp[i] = (double) mpfr_get_exp(R_i);
	}
    }
    mpfr_clear(R_i);
    mpfr_free_cache();
    UNPROTECT(1);
    return val;
}

//----------------------------------------------------------------------------

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

SEXP const_asMpfr(SEXP I, SEXP prec, SEXP rnd_mode)
{
    SEXP val;
    mpfr_t r;
    int i_p = asInteger(prec);
    R_mpfr_check_prec(i_p);
    mpfr_init2(r, i_p);

    switch(asInteger(I)) {
    case 1: mpfr_const_pi     (r, R_rnd2MP(rnd_mode)); break;
    case 2: mpfr_const_euler  (r, R_rnd2MP(rnd_mode)); break;
    case 3: mpfr_const_catalan(r, R_rnd2MP(rnd_mode)); break;
    case 4: mpfr_const_log2   (r, R_rnd2MP(rnd_mode)); break;
    default:
	error("invalid integer code {const_asMpfr()}"); /* -Wall */
    }

    FINISH_1_RETURN(r, val);
}

/** For functions     <logical>  <-  FUN(x = <mpfr>) :
 */
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

#define R_MPFRarray_Logic_Function(_FNAME, _MPFR_NAME)			\
SEXP _FNAME(SEXP x) {							\
    SEXP D = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym)),/* R list() */	\
       dim = PROTECT(GET_SLOT(x, Rmpfr_Dim_Sym)),			\
        dn = PROTECT(GET_SLOT(x, Rmpfr_Dimnames_Sym));			\
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
    setAttrib(val, R_DimSymbol,      duplicate(dim));			\
    setAttrib(val, R_DimNamesSymbol, duplicate(dn));			\
    UNPROTECT(4);							\
    return val;								\
}

R_MPFRarray_Logic_Function(R_mpfr_is_finite_A,   mpfr_number_p)
R_MPFRarray_Logic_Function(R_mpfr_is_infinite_A, mpfr_inf_p)
R_MPFRarray_Logic_Function(R_mpfr_is_integer_A,  mpfr_integer_p)
R_MPFRarray_Logic_Function(R_mpfr_is_na_A,       mpfr_nan_p)
R_MPFRarray_Logic_Function(R_mpfr_is_zero_A,     mpfr_zero_p)


SEXP R_mpfr_fac (SEXP n_, SEXP prec, SEXP rnd_mode)
{
    int n = length(n_), i, *nn;
    SEXP n_t, val = PROTECT(allocVector(VECSXP, n)); int nprot = 1;
    mpfr_rnd_t rnd = R_rnd2MP(rnd_mode);
    mpfr_t r_i;
    if(TYPEOF(n_) != INTSXP) {
	PROTECT(n_t = coerceVector(n_, INTSXP)); nprot++;/* or bail out*/
	nn = INTEGER(n_t);
    } else {
	nn = INTEGER(n_);
    }
    int i_p = asInteger(prec);
    R_mpfr_check_prec(i_p);
    mpfr_init2(r_i, i_p);
    for(i=0; i < n; i++) {
	// never happens when called from R:
	if(nn[i] < 0) error("R_mpfr_fac(%d): negative n.", nn[i]);
	mpfr_fac_ui(r_i, nn[i], rnd);
	SET_VECTOR_ELT(val, i, MPFR_as_R(r_i));
    }

    mpfr_clear(r_i);
    mpfr_free_cache();
    UNPROTECT(nprot);
    return val;
}


/** For functions   FUN(x = <mpfr>, y = <mpfr>) :
 */
#define R_MPFR_2_Numeric_Function(_FNAME, _MPFR_NAME)	\
SEXP _FNAME(SEXP x, SEXP y, SEXP rnd_mode) {		\
    SEXP xD = PROTECT(GET_SLOT(x, Rmpfr_Data_Sym));	\
    SEXP yD = PROTECT(GET_SLOT(y, Rmpfr_Data_Sym));	\
    mpfr_rnd_t rnd = R_rnd2MP(rnd_mode);		\
    int nx = length(xD), ny = length(yD), i,		\
	n = (nx == 0 || ny == 0) ? 0 : imax2(nx, ny);	\
    SEXP val = PROTECT(allocVector(VECSXP, n));		\
    mpfr_t R, x_i, y_i;					\
    mpfr_init(R); /* with default precision */		\
    mpfr_init(x_i); mpfr_init(y_i);			\
							\
    for(i=0; i < n; i++) {				\
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);		\
	R_asMPFR(VECTOR_ELT(yD, i % ny), y_i);		\
	_MPFR_NAME(R, x_i, y_i, rnd);			\
	SET_VECTOR_ELT(val, i, MPFR_as_R(R));		\
    }							\
							\
    mpfr_clear(R); mpfr_clear(x_i); mpfr_clear(y_i);	\
    mpfr_free_cache();					\
    UNPROTECT(3);					\
    return val;						\
}

R_MPFR_2_Numeric_Function(R_mpfr_atan2, mpfr_atan2)
R_MPFR_2_Numeric_Function(R_mpfr_hypot, mpfr_hypot)

#if (MPFR_VERSION >= MPFR_VERSION_NUM(3,2,0))
R_MPFR_2_Numeric_Function(R_mpfr_igamma, mpfr_gamma_inc)
#else
SEXP R_mpfr_igamma(SEXP x, SEXP y, SEXP rnd_mode) {
    error("mpfr_gamma_inc requires mpfr >= 3.2.0");
    return R_NilValue;
}
#endif

R_MPFR_2_Numeric_Function(R_mpfr_beta,  my_mpfr_beta)
R_MPFR_2_Numeric_Function(R_mpfr_lbeta, my_mpfr_lbeta)



/** For functions   FUN(x = <mpfr>, y = <integer>) :
 */
#define R_MPFR_2_Num_Long_Function(_FNAME, _MPFR_NAME)			\
SEXP _FNAME(SEXP x, SEXP y, SEXP rnd_mode) {				\
    SEXP xD, yt, val;							\
    mpfr_rnd_t rnd = R_rnd2MP(rnd_mode);				\
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
    PROTECT(val = allocVector(VECSXP, n));	nprot++;		\
    mpfr_init(x_i); /* with default precision */			\
									\
    for(i=0; i < n; i++) {						\
	R_asMPFR(VECTOR_ELT(xD, i % nx), x_i);				\
	_MPFR_NAME(x_i, (long) yy[i % ny], x_i, rnd);			\
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
