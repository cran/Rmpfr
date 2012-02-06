stopifnot(require("Rmpfr"))

stopifnot(chooseMpfr(1:10, 0) == 1,# failed earlier
	  chooseMpfr(20, 0:20) == choose(20, 0:20),
	  chooseMpfr(19, 0:20) == choose(19, 0:20),
	  chooseMpfr	(30, 4:30) * (-1)^(4:30) ==
	  chooseMpfr.all(30, k0=4, alternating=TRUE)
          )


cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()
