(TeX-add-style-hook "Maechler_Rmpfr"
 (lambda ()
    (LaTeX-add-labels
     "sec:intro"
     "sec:BinCoef")
    (TeX-add-symbols
     '("Ew" 1)
     '("fn" 1)
     "Rp"
     "CRAN"
     "W"
     "Ip")
    (TeX-run-style-hooks
     "MM-colors"
     "MM-slides"
     "mmVignette"
     "relsize"
     "inputenc"
     "utf8"
     "babel"
     "english"
     "latex2e"
     "beamer10"
     "dvipsnames"
     "pdflatex"
     "beamer"
     "log1exp"
     "BinCoef"
     "sumBinC"
     "RmpfrArith")))

