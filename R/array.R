## From an "mpfr" object make an mpfr(Array|Matrix) :

setMethod("dim", "mpfrArray", function(x) x@Dim)
setMethod("dimnames", "mpfrArray", function(x) x@Dimnames)

## 2 basic methods to construct "mpfr - arrays" ( mpfrArray | mpfrMatrix ) :

##' "mpfr" --> "mpfrArray"  --- basically  dim(<mpfr>) <- dd
mpfr2array <- function(x, dim, dimnames=NULL, check=FALSE) {
    if(check) stopifnot(is(x,"mpfr"))
    if(is.numeric(dim) && all(dim == (iv <- as.integer(dim)))) {
	if(check) {
	    cl <- if(length(iv) == 2) "mpfrMatrix" else "mpfrArray"
	    if(is.null(dimnames))
		new(cl, x, Dim = iv)
	    else new(cl, x, Dim = iv, Dimnames = dimnames)
	} else { ## faster, non-checking
	    r <- setDataPart(new(if(length(iv) == 2) "mpfrMatrix" else "mpfrArray"),
			     x, check=FALSE)
	    r@Dim <- iv
	    if(!is.null(dimnames)) r@Dimnames <- dimnames
	    r
	}
    }
    else if(is.null(dim))
	as.vector(x)
    else
	stop("invalid dimension specified")
}
setMethod("dim<-", signature(x = "mpfr", value = "ANY"),
	  function(x, value) mpfr2array(x, value))


mpfrArray <- function(x, precBits, dim = length(x), dimnames = NULL,
                      rnd.mode = c('N','D','U','Z'))
{
    dim <- as.integer(dim)
    rnd.mode <- toupper(rnd.mode)
    rnd.mode <- match.arg(rnd.mode)

    ml <- .Call(d2mpfr1_list, x, precBits, rnd.mode)
    vl <- prod(dim)
    if (length(x) != vl) {
        if (vl > .Machine$integer.max)
            stop("'dim' specifies too large an array")
        ml <- rep(ml, length.out = vl)
    }
    new(if(length(dim) == 2) "mpfrMatrix" else "mpfrArray",
        ml, Dim = dim,
	Dimnames = if(is.null(dimnames)) vector("list", length(dim))
		   else dimnames)
}

setAs("array", "mpfr", function(from) mpfr(from, 128L))

setMethod("dimnames<-", signature(x = "mpfrArray", value = "ANY"),
	  function(x, value) {
	      if(!is.list(value)) stop("non-list RHS")
	      if(length(value) != length(x@Dim))
		  stop("RHS (new dimnames) differs in length from dim(.)")
	      x@Dimnames <- value
	      x
	  })

setMethod("t", "mpfrMatrix",
	  function(x) {
	      d <- x@Dim; n <- d[1]; m <- d[2]
	      ## These are the indices to get the transpose of m {n x m} :
	      ## ind.t <- function(n,m)rep.int(1:n, rep(m,n)) + n*(0:(m-1))
	      x@Dim <- c(m,n)
	      x@Dimnames <- x@Dimnames[2:1]
	      ## faster than { x@.Data <- x@.Data[rep.int(1:n, rep(m,n)) + n*(0:(m-1))] ; x } :
	      setDataPart(x, getD(x)[rep.int(1:n, rep(m,n)) + n*(0:(m-1))], check=FALSE)
	  })
setMethod("t", "mpfr",
	  function(x) { # t(<n-vector>) |-->  {1 x n} matrix
	      r <- new("mpfrMatrix")
	      r@Dim <- c(1L, length(x))
              ## faster than  { r@.Data <- x@.Data ; r } :
	      setDataPart(r, getD(x), check=FALSE)
	  })

setMethod("aperm", signature(a="mpfrArray"),
	  function(a, perm, resize=TRUE) {
	      stopifnot(1 <= (k <- length(d <- a@Dim)))
	      if(missing(perm)) perm <- k:1
	      else stopifnot(length(perm <- as.integer(perm)) == k, 1 <= perm, perm <= k)
	      if(!resize)
		  stop("'resize != TRUE is not (yet) implemented for 'mpfrArray'")
	      a@Dim <- d[perm]
	      a@Dimnames <- a@Dimnames[perm]
	      ii <- c(aperm(array(1:prod(d), dim=d), perm=perm, resize=FALSE))
              ## faster than  { a@.Data <- a@.Data[ ii ] ; a } :
	      setDataPart(a, getD(a)[ ii ], check=FALSE)
	  })


## `` drop the  dim() part '' :
setMethod("as.vector", "mpfrArray", function(x) as(x, "mpfr"))
## a "vector" in  *one* sense at least, and "mpfr" does extend "vector":
setAs("mpfrArray", "vector", function(from) as(from, "mpfr"))

toNum <- function(from) {
    structure(.Call(mpfr2d, from),
	      dim = dim(from),
	      dimnames = dimnames(from))
}

setAs("mpfrArray", "array", toNum)

setAs("mpfrMatrix", "matrix", toNum)

setAs("mpfrArray", "matrix", function(from) {
    if(length(dim(from)) != 2)
	stop("dim(.) != 2  ==> cannot be coerced to 'matrix'")
    toNum(from)
})

print.mpfrArray <- function(x, digits = NULL, drop0trailing = FALSE, ...) {
##                                                            -----
## would like 'drop0... = TRUE', but that's only ok once we have a
## format() allowing to "jointly format a column"

    stopifnot(is(x, "mpfrArray"), is.null(digits) || digits >= 2)
    ## digits = NA --> the inherent precision of x will be used
    n <- length(x)
    ch.prec <-
	if(n >= 1) {
	    rpr <- range(.getPrec(x))
	    paste("of precision ", rpr[1],
		   if(rpr[1] != rpr[2]) paste("..",rpr[2]), " bits")
	}
    cl <- class(x)
    p0 <- function(...) paste(..., sep="")
    cat(p0("'",cl,"'"), "of dim(.) = ",
        p0("(",paste(x@Dim, collapse=", "),")"),
        ch.prec, "\n")
    if(n >= 1) {
        ## FIXME: really need a 'format' method for mpfrArrays
        ## -----  which properly alings columns !!
        fx <- format(x, digits=digits, drop0trailing=drop0trailing)
        dim(fx) <- dim(x)
        dimnames(fx) <- dimnames(x)
	print(fx, ..., quote = FALSE)
    }
    invisible(x)
}
setMethod(show, "mpfrArray", function(object) print.mpfrArray(object))

## FIXME : should happen in C, where we could "cut & paste" much of
## -----  do_matprod() and matprod() from ~/R/D/r-devel/R/src/main/array.c
##/* "%*%" (op = 0), crossprod (op = 1) or tcrossprod (op = 2) */
.matmult.R <- function(x,y, op = 0)
{
    if(!(is.numeric(x) || is(x,"mpfr")))
	stop("'x' must be numeric or mpfr(Matrix)")
    sym <- missing(y)
    if (sym && (op > 0)) y <- x
    else if(!(is.numeric(y) || is(y,"mpfr")))
	stop("'y' must be numeric or mpfr(Matrix)")
    ldx <- length(dx <- dim(x))
    ldy <- length(dy <- dim(y))
    ## "copy, paste & modify" from  do_matprod():
    if (ldx != 2 && ldy != 2) {		#* x and y non-matrices */
	if (op == 0) {
	    nrx <- 1L; ncx <- length(x)
	} else {
	    nrx <- length(x); ncx <- 1L
	}
	nry <- length(y)
	ncy <- 1L
    }
    else if (ldx != 2) {		#* x not a matrix */
	nry <- dy[1]
	ncy <- dy[2]
	nrx <- ncx <- 0L
	if (op == 0) {
	    if (length(x) == nry) {	#* x as row vector */
		nrx <- 1L
		ncx <- nry # == length(x)
	    }
	    else if (nry == 1) {	#* x as col vector */
		nrx <- length(x)
		ncx <- 1L # == nry
	    }
	}
	else if (op == 1) { #* crossprod
	    if (length(x) == nry) {	#* x is a col vector */
		nrx <- nry # = length(x)
		ncx <- 1L
	    }
	}
	else { # op == 2: tcrossprod
	    if (length(x) == ncy) {	#* x as row vector */
		nrx <- 1L
		ncx <- ncy # == length(x)
	    }
	    else if (ncy == 1) {	#* x as col vector */
		nrx <- length(x)
		ncx <- 1L # == ncy
	    }
	}
    }
    else if (ldy != 2) {		#* y not a matrix */
	nrx <- dx[1]
	ncx <- dx[2]
	nry <- ncy <- 0L
	if (op == 0) {
	    if (length(y) == ncx) {	#* y as col vector */
		nry <- ncx # = length(y)
		ncy <- 1L
	    }
	    else if (ncx == 1) {	#* y as row vector */
		nry <- 1L # = ncx
		ncy <- length(y)
	    }
	}
	else if (op == 1) { #* crossprod
	    if (length(y) == nrx) {	#* y is a col vector */
		nry <- nrx # = length(y)
		ncy <- 1L
	    }
	}
	else { # op == 2: tcrossprod	y is a col vector
	    nry <- length(y)
	    ncy <- 1L
	}
    }
    else {				#* x and y matrices */
	nrx <- dx[1]
	ncx <- dx[2]
	nry <- dy[1]
	ncy <- dy[2]
    }
    ##* nr[ow](.) and nc[ol](.) are now defined for x and y */

    z <- new("mpfrMatrix")
    z0 <- as(0, "mpfr")

    if (op == 0) { ## %*%
	if (ncx != nry) stop("non-conformable arguments")

	z@Dim <- c(nrx, ncy)
	z@.Data <- vector("list", nrx*ncy)
	if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
	    for(i in 1:nrx) {
		j <- 0L:(ncx - 1L)
		for (k in 0L:(ncy - 1L))
		    z[i + k * nrx] <- sum(x[i + j * nrx] * y[1L+ j + k * nry])
	    }
	else	    #/* zero-extent operations should return zeroes */
	    for(i in seq_len(nrx*ncy)) z[i] <- z0
    }
    else if (op == 1) { ## crossprod() :  x' %*% y
	if (nrx != nry) stop("non-conformable arguments")

	z@Dim <- c(ncx, ncy)
	z@.Data <- vector("list", ncx*ncy)
	if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
	    for(i in 0L:(ncx - 1L)) {
		j <- 1L:nrx
		for (k in 0L:(ncy - 1L))
		    z[1L +i + k * ncx] <- sum(x[j + i * nrx] * y[j + k * nry])
	    }
	else
	    for(i in seq_len(ncx*ncy)) z[i] <- z0
    }
    else { ## op == 2 :	 tcrossprod() :	 x %*% y'
	if (ncx != ncy) stop("non-conformable arguments")

	z@Dim <- c(nrx, nry)
	z@.Data <- vector("list", nrx*nry)
	if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
	    for(i in seq_len(nrx)) {
		j <- 0L:(ncx - 1L)
		for (k in 0L:(nry - 1L))
		    z[i + k * nrx] <- sum(x[i + j * nrx] * y[1L +k + j * nry])
	    }
	else
	    for(i in seq_len(nrx*nry)) z[i] <- z0
    }
    z
} ## .matmult.R()

## "FIXME"?  make working also with "Matrix" class matrices ..
##                             ----------------------------


setMethod("%*%", signature(x = "mpfrMatrix", y = "mpfrMatrix"),
	  function(x,y) .matmult.R(x,y, op= 0))
setMethod("%*%", signature(x = "mpfrMatrix", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 0))
setMethod("%*%", signature(x = "mpfr", y = "mpfrMatrix"),
	  function(x,y) .matmult.R(x,y, op= 0))
setMethod("%*%", signature(x = "mpfr", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 0))
## These cover vectors, etc (!) :
setMethod("%*%", signature(x = "mpfr", y = "Mnumber"),
	  function(x,y) .matmult.R(x,y, op= 0))
setMethod("%*%", signature(x = "Mnumber", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 0))


setMethod("crossprod", signature(x = "mpfrMatrix", y = "mpfrMatrix"),
	  function(x,y) .matmult.R(x,y, op= 1))
setMethod("crossprod", signature(x = "mpfrMatrix", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 1))
setMethod("crossprod", signature(x = "mpfr", y = "mpfrMatrix"),
	  function(x,y) .matmult.R(x,y, op= 1))
setMethod("crossprod", signature(x = "mpfr", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 1))
setMethod("crossprod", signature(x = "mpfr", y = "Mnumber"),
	  function(x,y) .matmult.R(x,y, op= 1))
setMethod("crossprod", signature(x = "Mnumber", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 1))
## one argument-case:
setMethod("crossprod", signature(x = "mpfr", y = "missing"),
	  function(x,y) .matmult.R(x,x, op= 1))

setMethod("tcrossprod", signature(x = "mpfrMatrix", y = "mpfrMatrix"),
	  function(x,y) .matmult.R(x,y, op= 2))
setMethod("tcrossprod", signature(x = "mpfrMatrix", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 2))
setMethod("tcrossprod", signature(x = "mpfr", y = "mpfrMatrix"),
	  function(x,y) .matmult.R(x,y, op= 2))
setMethod("tcrossprod", signature(x = "mpfr", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 2))
setMethod("tcrossprod", signature(x = "mpfr", y = "Mnumber"),
	  function(x,y) .matmult.R(x,y, op= 2))
setMethod("tcrossprod", signature(x = "Mnumber", y = "mpfr"),
	  function(x,y) .matmult.R(x,y, op= 2))
## one argument-case:
setMethod("tcrossprod", signature(x = "mpfr", y = "missing"),
	  function(x,y) .matmult.R(x,x, op= 2))


.mpfrA.subset <- function(x,i,j, ..., drop) {
    nA <- nargs()
    if(getOption("verbose"))
        message(sprintf("nargs() == %d  mpfrArray indexing ... ", nA))

    r <- getD(x)
    if(nA == 2) ## A[i]
        return(new("mpfr", r[i]))
    ## else: nA != 2 : nA > 2 -
    dim(r) <- (dx <- dim(x))
    dimnames(r) <- dimnames(x)
    r <- r[i,j, ..., drop=drop]
    if(drop & is.null(dim(r)))
        new("mpfr", r)
    else {
        D <- if(is.null(dr <- dim(r))) # ==> drop is FALSE; can this happen?
            rep.int(1L, length(r)) else dr
        x@Dim <- D
        x@Dimnames <- if(is.null(dn <- dimnames(r)))
            vector("list", length(D)) else dn
	if(length(D) == 2 && class(x) != "mpfrMatrix")
	    ## low-level "coercion" from mpfrArray to *Matrix :
	    attr(x,"class") <- getClass("mpfrMatrix")@className
        attributes(r) <- NULL
        setDataPart(x, r, check=FALSE)
    }
}

## "["
setMethod("[", signature(x = "mpfrArray", i = "ANY", j = "ANY", drop = "ANY"),
          .mpfrA.subset)

## this signature needs a method here, or it triggers the one for "mpfr"
setMethod("[", signature(x = "mpfrArray", i = "ANY", j = "missing",
                         drop = "missing"),
          .mpfrA.subset)


.mA.subAssign <- function(x,i,j,..., value, n.a, isMpfr)
{
    ## n.a :=== nargs() -- in the calling "[<-" method --
    r <- getD(x)
    if(n.a >= 4) {
	## A[i,j]  /  A[i,]  /	A[,j]	but not A[i]
	## A[i,j,k] <- v : n.a == 5
	dim(r) <- dim(x)
	dimnames(r) <- dimnames(x)
	if(!isMpfr)
	   value <- mpfr(value, precBits =
			 pmax(getPrec(value),
			      .getPrec(if(n.a == 4) r[i,j] else r[i,j, ...]))
			 )
	vD <- getD(value)
	if(n.a == 4) {
	    r[i,j] <- vD
	} else { ## n.a >= 5
	    r[i, j, ...] <- vD
	}
	attributes(r) <- NULL
    }
    else if(n.a == 3) { ##  A [ i ] <- v
	if(!isMpfr)
	    value <- mpfr(value, precBits = pmax(getPrec(value), .getPrec(r[i])))

	r[i] <- value
    } else { ## n.a <= 2
	stop(sprintf("nargs() == %d  mpfrArray[i,j] <- value  IMPOSSIBLE?",
		     n.a))
    }
    setDataPart(x, r, check=FALSE)
}## .mA.subAssign

## "[<-" :
## -------
## E.g., for A[1,,2] <- V
## these are to trigger before the  ("mpfr", i,j, "mpfr")  [ ./mpfr.R ] does
for(it in c("ANY", "missing"))
    for(jt in c("ANY", "missing"))
    setReplaceMethod("[", signature(x = "mpfrArray", i = it, j = jt, value = "mpfr"),
		     function(x,i,j,..., value)
		     .mA.subAssign(x,i=i,j=j,...,value=value,
				   n.a=nargs(), isMpfr = TRUE))
## non-"mpfr" value
for(it in c("ANY", "missing"))
    for(jt in c("ANY", "missing"))
    setReplaceMethod("[", signature(x = "mpfrArray", i = it, j = jt, value = "ANY"),
		     function(x,i,j, ..., value)
		     .mA.subAssign(x,i=i,j=j,...,value=value,
				   n.a=nargs(), isMpfr = FALSE))
rm(it,jt)

###-----------

setGeneric("cbind", signature = "...")# -> message about override & deparse.level
setGeneric("rbind", signature = "...")

setMethod("cbind", "Mnumber",
	  function(..., deparse.level = 1) {
	      args <- list(...)
	      if(all(sapply(args, is.atomic)))
		  return( base::cbind(..., deparse.level = deparse.level) )
	      ## else: at least one is "mpfr(Matrix/Array)"

	      if(any(sapply(args, is.character))) {
		  ## result will be  <character> matrix !
		  isM <- sapply(args, is, class2 = "mpfr")
		  args[isM] <- lapply(args[isM], as, Class = "character")
		  return(do.call(base::cbind,
				 c(args, list(deparse.level=deparse.level))))

	      } else if(any(sapply(args, is.complex))) {
		  ## result will be  <complex> matrix;
		  ## in the future <complex_mpfr>  ???

		  stop("cbind(...) of 'complex' and 'mpfr' objects is not implemented")
		  ## give at least warning !!
              }
              ## else
	      L <- function(a) if(is.numeric(n <- nrow(a))) n else length(a)
	      W <- function(a) if(is.numeric(n <- ncol(a))) n else 1L
	      ## the number of rows of the result :
	      NR <- max(lengths <- sapply(args, L))
	      NC <- sum(widths	<- sapply(args, W))
	      r <- setDataPart(new("mpfrMatrix"), vector("list", NR*NC))
	      r@Dim <- as.integer(c(NR, NC))
	      if(deparse.level >= 1 && !is.null(nms <- names(widths)))
		  r@Dimnames[[2]] <- nms
	      j <- 0
	      prec <- .Machine$double.digits
	      for(ia in seq_along(args)) {
		  w <- widths[ia]
		  a <- args[[ia]]
		  if(is(a,"mpfr")) {
		      prec <- max(prec, .getPrec(a))
		  } else { ## not "mpfr"
		      a <- mpfr(a, prec)
		  }
		  if((li <- lengths[ia]) != 1 && li != NR) { ## recycle
		      if(!is.null(dim(a)))
			  stop("number of rows of matrices must match")
		      ## else
		      if(NR %% li)
			  warning("number of rows of result is not a multiple of vector length")
		      a <- a[rep(seq_len(li), length.out = NR)]
		  }
		  r[, j+ 1:w] <- a
		  j <- j + w
	      }
	      r
	  })

setMethod("rbind", "Mnumber",
	  function(..., deparse.level = 1) {
	      args <- list(...)
	      if(all(sapply(args, is.atomic)))
		  return( base::rbind(..., deparse.level = deparse.level) )
	      ## else: at least one is "mpfr(Matrix/Array)"

	      if(any(sapply(args, is.character))) {
		  ## result will be  <character> matrix !
		  isM <- sapply(args, is, class2 = "mpfr")
		  args[isM] <- lapply(args[isM], as, Class = "character")
		  return(do.call(base::rbind,
				 c(args, list(deparse.level=deparse.level))))

	      } else if(any(sapply(args, is.complex))) {
		  ## result will be  <complex> matrix;
		  ## in the future <complex_mpfr>  ???

		  stop("rbind(...) of 'complex' and 'mpfr' objects is not implemented")
		  ## give at least warning !!
	      }
              ## else
 	      L <- function(a) if(is.numeric(n <- nrow(a))) n else 1L
	      W <- function(a) if(is.numeric(n <- ncol(a))) n else length(a)
	      ## the number of rows of the result :
	      NR <- sum(lengths <- sapply(args, L))
	      NC <- max(widths	<- sapply(args, W))

	      r <- setDataPart(new("mpfrMatrix"), vector("list", NR*NC))
	      r@Dim <- as.integer(c(NR, NC))
	      if(deparse.level >= 1 && !is.null(nms <- names(widths)))
		  r@Dimnames[[1]] <- nms
	      i <- 0
	      prec <- .Machine$double.digits
	      for(ia in seq_along(args)) {
		  le <- lengths[ia]
		  a <- args[[ia]]
		  if(is(a,"mpfr")) {
		      prec <- max(prec, .getPrec(a))
		  } else { ## not "mpfr"
		      a <- mpfr(a, prec)
		  }

		  if((wi <- widths[ia]) != 1 && wi != NC) { ## recycle
		      if(!is.null(dim(a)))
			  stop("number of rows of matrices must match")
		      ## else
		      if(NC %% wi)
			  warning("number of columns of result is not a multiple of vector length")
		      a <- a[rep(seq_len(wi), length.out = NC)]
		  }
		  r[i+ 1:le, ] <- a
		  i <- i + le
	      }
	      r
	  })
