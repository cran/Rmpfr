#### All Class Definitions in package  "Rmpfr"

### NB:	 Use   /usr/local/app/R/R_local/src/Brobdingnag/R/brob.R
###					    -----------
### as a partial role image

setClass("mpfr1", ## a single Multi-precision float number
	 representation(prec = "integer", # precision in bits
			exp = "integer",  # exponent
			sign= "integer",  # signum
			d = "integer"),	  # the mantissa as a vector of long
	 validity = function(object) {
	     if(length(pr <- object@prec) != 1 || is.na(pr) || pr < 2)
		 "invalid 'prec' slot"
	     else if(length(ex <- object@prec) != 1)
		 "invalid 'exp' slot"
	     else if(length(sig <- object@sign) != 1 || is.na(sig) ||
		     abs(sig) > 1)
		 "invalid 'sign' slot"
	     else if(length(d <- object@d) != ceiling(pr / 32))
		 "length('d' slot) does not match 'prec'"
	     else TRUE
	 })

setClass("mpfr", ## a *vector* of "mpfr1", i.e., multi-precision float numbers
	 contains = "list", ## of "mpfr1" entries:
	 validity = function(object) {
	     ## should be fast ( ==> not using	is(., "mpfr1") ) :
	     if(all(vapply(object@.Data, class, "") == "mpfr1"))
		 return(TRUE)
	     ## else
		 "Not all components are of class 'mpfr1'"
	 })

setClass("mpfrArray", ## mpfr + "dim" + dimnames
	 contains = "mpfr",
	 representation = list(Dim = "integer", Dimnames = "list"),
	 prototype = prototype(new("mpfr"), Dim= 0L,
			       Dimnames = list(NULL)),
	 validity = function(object) {
	     if(length(object) != prod(D <- object@Dim))
		 "Dimension does not match length()"
	     else if(length(DN <- object@Dimnames) != length(D))
		 "Dimnames must have same length as 'Dim'"
	     else if(any(hasN <- !vapply(DN, is.null, NA)) &&
		     any((lDN <- vapply(DN[hasN], length, 1L)) != D[hasN]))
		 "length of some 'Dimnames' do not match 'Dim'"
	     else
		 TRUE
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

## S3 classes from package gmp --- to be used in signatures {together with "mpfr"}:
setOldClass("bigz")
setOldClass("bigq")

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
		       is(object, "mpfr")) return(TRUE)
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
                   is(object, "mpfr")) return(TRUE)
                ## else
                "Not a valid 'mNumber' class object"
		})
