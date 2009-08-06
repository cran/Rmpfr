require("Rmpfr")

### Simple basic examples of creation of   "mpfr"  objects

pi. <- Const("pi", prec = 260)
pi. # nicely prints 80 digits [260 * log10(2) ~= 78.3 ~ 80]

## This is TRUE for 0 and -0 :
Zero <- mpfr(c(0,1/-Inf), 20)
stopifnot(mpfr.is.0(Zero))
Zero

d.spec <- c(0,NA,NaN,Inf,-Inf)
(spec <- mpfr(d.spec, 3))
stopifnot(identical(is.na(spec), is.na(d.spec)),
          identical(is.finite(spec), is.finite(d.spec)),
          identical(is.infinite(spec), is.infinite(d.spec)),
          identical(format(spec),
                    c("0.00", "NaN", "NaN", "Inf", "-Inf")))

x <- c(-12, 1:3 * pi)
sss <- mpfr(x, 100)
validObject(sss)
sss
sss2 <- sss * sss
stopifnot(identical(sss2, sss * x),
          identical(sss2, x * sss),
          sss ^ 2 == sss2)
## and go back {not sure if  identical() is guaranteed here, but it seems...}:
stopifnot(identical(x, as(sss, "numeric")))

(cs <- as(sss, "character"))

y <- c(0, 100,-10, 1.25, -2.5,
       x * c(1,100,1e5,1e20),
       x / 100^(1:4))
(Y <- mpfr(y, 100))
cbind(y, as.data.frame(.mpfr2str(Y, 20))[,c("exp","str")])

eps8 <- 8 * .Machine$double.eps
## checking  mpfr -> character -> mpfr:
stopifnot(all.equal(y, as.numeric(format(Y, digits=20)), tol= eps8),
          all.equal(Y, as(format(Y), "mpfr"), tol= eps8))

## More  character -> mpfr  checking :
## from   echo 'scale=200; 4*a(1)' | bc -l :
cpi <- "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196"
pi. <- Const("pi", prec=667)
stopifnot(cpi == format(mpfr(cpi, prec=667), digits=201),
          all.equal(pi., as(cpi, "mpfr")),
          all.equal(pi., as(cpi, "mpfr"), tol = 1e-200))

## Check double -> mpfr -> character -> double :
rSign <- function(n) sample(c(-1,1), size = n, replace=TRUE)
N <- function(x) as.numeric(x)
for(n in 1:40) {
    cat(if(n %% 10)"." else n)
    x. <- rSign(100) * rlnorm(100)
    X. <- mpfr(x., precBits = 120L)
    stopifnot(all.equal(x., N(format(X., digits=20)), tol = eps8)
              , all.equal(x., N(log(exp(X.))), tol = 32*eps8)
    )
}; cat("\n")

X. <- X.[!mpfr.is.0(X.)]
stopifnot(all( X./X. == 1)) # TRUE

u <- mpfr(as.raw(0:100))
stopifnot(0:100 == u,
	  all.equal(u, mpfr(0:100, prec = 8), tol = 0),
	  0:1 == mpfr(1:2 %% 2 == 0))
