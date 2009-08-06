Const <- function(name = c("pi", "gamma", "catalan"), prec = 120L)
{
    stopifnot(is.numeric(prec))
    i <- pmatch(name, eval(formals()$name))
    new("mpfr", list(.Call("const_asMpfr", i, prec, PACKAGE="Rmpfr")))
}
## fails here; must happen *after* dyn.load ... : Pi <- Const("pi")
