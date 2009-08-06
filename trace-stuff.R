trace("[<-", browser, exit=browser,
      signature=c(x="mpfrArray",i="numeric",j="numeric",value="mpfr"))
trace("[<-", browser, exit=browser,
      signature=c(x="mpfrArray",i="missing",j="numeric",value="mpfr"))
trace("[<-", browser, exit=browser,
      signature=c(x="mpfrArray",i="numeric",j="missing",value="mpfr"))
trace("[<-", browser, exit=browser,
      signature=c(x="mpfrArray",i="missing",j="missing",value="mpfr"))

## Revert:

untrace("[<-",
      signature=c(x="mpfrArray",i="numeric",j="numeric",value="mpfr"))
untrace("[<-",
      signature=c(x="mpfrArray",i="missing",j="numeric",value="mpfr"))
untrace("[<-",
      signature=c(x="mpfrArray",i="numeric",j="missing",value="mpfr"))
untrace("[<-",
      signature=c(x="mpfrArray",i="missing",j="missing",value="mpfr"))

## --- same with s/Array/Matrix : -------------
trace("[<-", browser, exit=browser,
      signature=c(x="mpfrMatrix",i="numeric",j="numeric",value="mpfr"))
trace("[<-", browser, exit=browser,
      signature=c(x="mpfrMatrix",i="missing",j="numeric",value="mpfr"))
trace("[<-", browser, exit=browser,
      signature=c(x="mpfrMatrix",i="numeric",j="missing",value="mpfr"))
trace("[<-", browser, exit=browser,
      signature=c(x="mpfrMatrix",i="missing",j="missing",value="mpfr"))

## Revert:

untrace("[<-",
      signature=c(x="mpfrMatrix",i="numeric",j="numeric",value="mpfr"))
untrace("[<-",
      signature=c(x="mpfrMatrix",i="missing",j="numeric",value="mpfr"))
untrace("[<-",
      signature=c(x="mpfrMatrix",i="numeric",j="missing",value="mpfr"))
untrace("[<-",
      signature=c(x="mpfrMatrix",i="missing",j="missing",value="mpfr"))
