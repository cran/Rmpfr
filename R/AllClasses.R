#### All Class Definitions in package  "Rmpfr"

### Historically, we have used /usr/local/app/R/R_local/src/Brobdingnag/R/brob.R
### as partial role image                		    -----------

## NB: Most MPFR numbers are regular (C: mpfr_regular_p )
## ==  ---> then 'd' is well defined.
##     If they are not, i.e., it is in {0, NaN, +/- Inf}, then 'exp' shows this
##     and 'd' (in mpfr a pointer to the limbs) is *not* used and not defined.
## Since Jan.2018, Rmpfr 0.7-0, we reflect this by using a 0-length 'd' slot

setClass("mpfr1", ## a single Multi-precision float number
	 slots = c(prec = "integer", # precision in bits
                   exp = "integer",  # exponent
                   sign= "integer",  # signum
                   d = "integer"),   # the mantissa as a vector of (32 bit) integers
	 validity = function(object) {
	     gmp.numb <- .mpfr_gmp_numbbits() # 32 or 64
	     if(length(pr <- object@prec) != 1 || is.na(pr) || pr < 2)
		 "invalid 'prec' slot"
	     else if((lex <- length(ex <- object@exp)) != 2 && gmp.numb == 64)
		 "invalid 'exp' slot for 64-bit gmp.numbbits: must have length 2"
	     else if(lex != 1 && gmp.numb == 32)
		 "invalid 'exp' slot for 32-bit gmp.numbbits: must have length 1"
	     else if(length(sig <- object@sign) != 1 || is.na(sig) || abs(sig) > 1)
		 "'sign' slot not in {-1,1} is invalid"
	     else {
		 nd <- length(d <- object@d)
		 if(nd) { ## "regular"
		     need.d <- ceiling(pr / 32)
		     if((gmp.numb == 32 && nd != need.d) ||
			(gmp.numb == 64 && !any((nd - need.d) == 0:1)))
			 "length('d' slot) does not match 'prec'"
		     else TRUE
		 } else ## not regular: valid if exp slot shows so
		     if(gmp.numb == 64) {  ## ex of length 2
			 if((is.na(ex[2]) && any(ex[[1]] == (1:3))) ||  ## mpfr 3.1.5, Fedora 26
			    (ex[1] == ex[2] && any(ex[1]+2^31 == 1:3))) ## mpfr 3.1.3, Windows
			     TRUE
			 else
			     "'exp' slot invalid for non-regular number (64b, length(d) == 0)"
		     } else { ## gmp.numb == 32: 'exp' slot of length one
			 if(any(ex+2^31 == 1:3))
			     TRUE
			 else
			     "'exp' slot invalid for non-regular number (32b, length(d) == 0)"
		     }
	     }
	 })

setClass("mpfr", ## a *vector* of "mpfr1", i.e., multi-precision float numbers
	 contains = "list", ## of "mpfr1" entries:
	 validity = function(object) {
	     ## should be fast ( ==> not using	is(., "mpfr1") ) :
	     if(all(lengths(cls <- lapply(object@.Data, class)) == 1L) && 
		all(unlist(cls) == "mpfr1"))
		 return(TRUE)
	     ## else
		 "Not all components are of class 'mpfr1'"
	 })

setClass("mpfrArray", ## mpfr + "dim" + dimnames
	 contains = "mpfr",
	 slots = c(Dim = "integer", Dimnames = "list"),
	 prototype = prototype(new("mpfr"), Dim= 0L),
	 validity = function(object) {
	     if(length(object) != prod(D <- object@Dim))
		 "Dimension does not match length()"
	     else if(length(DN <- object@Dimnames) != length(D))
		 "Dimnames must have same length as 'Dim'"
	     else if(any(hasN <- !vapply(DN, is.null, NA)) &&
		     any(lengths(DN[hasN]) != D[hasN]))
		 "length of some 'Dimnames' do not match 'Dim'"
	     else
		 TRUE
	 })

setMethod("initialize", "mpfrArray", function(.Object, ..., Dim, Dimnames)
      {
	  if(!missing(Dim))
	      .Object@Dim <- as.integer(Dim)
	  k <- length(.Object@Dim)
	  if(missing(Dimnames))
	      .Object@Dimnames <- rep(list(NULL), k)
          else if(length(Dimnames) != k) {
              message(sprintf("in initialize: length(Dimnames) != k = %d;  setting to NULL",k))
	      .Object@Dimnames <- rep(list(NULL), k)
          } else
	      .Object@Dimnames <- Dimnames
	  callNextMethod()
      })


setClass("mpfrMatrix",
	 contains = "mpfrArray",
	 prototype = prototype(new("mpfrArray"),
			       Dim= c(0L,0L),
			       Dimnames = list(NULL, NULL)),
	 validity = function(object) {
	     if(length(object@Dim) != 2L)
		 "'Dim' is not of length 2"
	     else TRUE
	 })

## "atomic vectors" (-> ?is.atomic ) -- exactly as in "Matrix":
## ---------------
setClassUnion("atomicVector", ## "double" is not needed, and not liked by some
	      members = c("logical", "integer", "numeric",
			  "complex", "raw", "character"))

## This is tricky ...
## With the following class,  arrays/matrices  are covered as
## they are also with "vector" already. *However*, they are
## *not* made into vectors in method dispatch,
## which they would be if we used simply "vector"
setClassUnion("array_or_vector",
	      members = c("array", "matrix", "atomicVector"))
## However (FIXME?), the above is too large: "matrix" extends "vector"
## and that has "character", "list", ...

## For this class, we want to define  '...' methods for cbind & rbind :
## FIXME(?): "array_or_vector" also contains "character"
##         (and even if it wouldn't, a "matrix" could have "character" entries!)
setClassUnion("Mnumber",
	      members = c("array_or_vector", # *but* must be numeric-like
	      "mpfr", "mpfrArray", "mpfrMatrix",
	      ## from package 'gmp' :
	      "bigz", "bigq"))

if(FALSE) { ## if we did this, then ... {see below}
    setValidity("Mnumber",
		function(object) {
		    if(is.numeric(object) ||
		       is.logical(object) ||
		       is.mpfr(object)) return(TRUE)
		    ## else
		    "Not a valid 'Mnumber' class object"
		})

    ## ...., then, the following would fail (!)
    validObject( new("character", LETTERS) )
}

###----- Simpler {without 'matrix' -> 'character' ...} -------------------------
###
setClassUnion("numericVector", members = c("logical", "integer", "numeric"))

setClassUnion("mNumber",
	      members = c("numericVector",
	      "mpfr", "mpfrArray", "mpfrMatrix",
	      ## from package 'gmp' :
	      "bigz", "bigq"))
setValidity("mNumber",
            function(object) {
                if(is.numeric(object) ||
                   is.logical(object) ||
                   is.mpfr(object)) return(TRUE)
                ## else
                "Not a valid 'mNumber' class object"
		})
