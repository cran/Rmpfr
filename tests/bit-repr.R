stopifnot(require("Rmpfr"))

### Try to look at the internal bit-representation of the limbs

.limbs <- function(x) {
    stopifnot(is(x, "mpfr"))
    lapply(x, slot, "d") # not sapply() each can have different prec. & #{limbs}
}
.expo <- function(x) {
    stopifnot(is(x, "mpfr"))
    sapply(x, slot, "exp")
}

Bits <- function(x) {
    L <- .limbs(x)# list(length n) each of "k(prec)" 32-bit ints
    hasNA <- any(iNA <- sapply(lapply(L, is.na), any)) # iNA: TRUE if there's an NA
    ## need to catch them later
    CC <- function(ch) paste(ch, collapse="")
    hex <- sapply(L, function(I) CC(sprintf("%x", rev(I))))
    if(hasNA) hex[iNA] <- NA_character_
    hex <- strsplit(hex, NULL)

    db <- t(expand.grid(0:1,0:1,0:1,0:1, KEEP.OUT.ATTRS=FALSE)[,4:1])
    storage.mode(db) <- "character" # "0" or "1"
    dimnames(db) <- list(NULL, c(paste(0:9), letters[1:6]))
    ## db is  4 x 16  matrix  with col.names "0" "1" .. "9" "a" "b" ... "f"

    ex <- .expo(x)
    if(is.matrix(ex))
        ## 64-bit case: exponent is long == two ints
        ## -----------  but currently the 2nd int is always 0 (or NA for '0')
        ex <- ex[1,]
    pat <- paste("(", sapply(pmax(0, ex),
                             function(n) CC(rep.int(".", n))),
                 ")0+$", sep="")
    pat <- ifelse(iNA, NA_character_, pat)

    if(hasNA) {
	r <- as.list(iNA)
	r[!iNA] <- lapply(hex[!iNA], function(ch) CC(as.vector(db[,ch])))
	r[iNA ] <- NA_character_
        ## now keep correct number of trailing zeros :
        r[!iNA] <- lapply(which(!iNA), function(i) sub(pat[i], "\\1", r[[i]]))
        unlist(r)
    }
    else {
	r <- lapply(hex, function(ch) CC(as.vector(db[,ch])))
        ## now keep correct number of trailing zeros :
        sapply(seq_along(r), function(i) sub(pat[i], "\\1", r[[i]]))
    }

}

x <- mpfr(c(3:5,11:15, 59, 125:127, 1025), 64)
x
data.frame(x= as.numeric(x), I(Bits(x)))

x <- mpfr(c(-20:30),64)
x <- x[x != 0] # mpfr(0, *) has "random" bits
data.frame(x= as.numeric(x), I(Bits(x)))

## pi, in varying number of bits :
p. <- round(pi* 2^c(10,16,5*(4:8)))
dput(p.)
p <- c(mpfr(c(3217, 205887, 3294199, 105414357,
              3373259426, 107944301636, 3454217652358), 64),
       Const("pi", 64))
Bits(p)


cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()
