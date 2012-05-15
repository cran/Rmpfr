R_PROG = R
RNWSRC = log1exp BinCoef loglikFrank sumBinC RmpfrArith integrateBeta
DRIVERS= Maechler_Rmpfr
TEXSRC =
SWEAVE = SweaveSlides.sty # SweaveSlides.cfg

XtrTARGETS = $(DRIVERS:=-4.pdf) $(DRIVERS:=-4a4.pdf) $(RNWSRC:_R.txt) allR

## 'make' should make the  slides file only
# slides: $(DRIVERS:=D.pdf)
slides: $(DRIVERS:=.pdf)

all: slides $(XtrTARGETS)

allR: $(RNWSRC:=.R)
	touch allR
allTeX: $(RNWSRC:=.tex)
	touch allTeX
%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex
	${R_PROG} CMD texi2pdf $<

$(DRIVERS:=.pdf): $(DRIVERS:=.tex) allTeX
	TEXINPUTS=${TPUTS} texi2dvi --pdf $<

dist: clean
## ^ currently not more
clean:
	$(RM) $(RNWSRC:=.tex) \
          Rplots.ps  Rplots.pdf $(DRIVERS:=*.out) allR \
	  figs/*.pdf *-*-*.pdf \
          *.aux *.log *.toc *.nav *.snm *.vrb allTeX

realclean: clean
	$(RM) $(DRIVERS:=.pdf) $(XtrTARGETS)



%.pdf: %.tex
	TEXINPUTS=${TPUTS} texi2dvi --pdf $<
#	texclean $<
%H.pdf: %H.tex
	TEXINPUTS=${TPUTS} texi2dvi --pdf $<
#	texclean $<

%D.pdf: %D.tex
	TEXINPUTS=${TPUTS} texi2dvi --pdf $<
#	texclean $<

%-4.tex: %.pdf
	@(echo '#documentclass[letterpaper,landscape]{article}' ; \
	  echo '#usepackage{pdfpages}' ; \
	  echo '#begin{document}' ; \
	  echo '#includepdf[nup=2x2,pages=-]{'$*'.pdf}' ; \
	  echo '#end{document}' ) | sed 's/^#/\\/' > $@

%-4.pdf: %-4.tex
	texi2dvi --pdf $<

%-2.tex: %.pdf
	@echo '#documentclass[letterpaper]{article}' ; \
	  echo '#usepackage{pdfpages}' ; \
	  echo '#begin{document}' ; \
	  echo '#includepdf[nup=1x2,pages=-]{'$*'.pdf}' ; \
	  echo '#end{document}' ) | sed 's/^#/\\/' > $@


%-2.pdf: %-2.tex
	texi2dvi --pdf $<
#	texclean $<


# A4 paper version: is   echo '-4' | sed s/letter/a4/
%-4a4.tex: %.pdf
	@(echo '#documentclass[a4paper,landscape]{article}' ; \
	  echo '#usepackage{pdfpages}' ; \
	  echo '#begin{document}' ; \
	  echo '#includepdf[nup=2x2,pages=-]{'$*'.pdf}' ; \
	  echo '#end{document}' ) | sed 's/^#/\\/' > $@


%-4a4.pdf: %-4a4.tex
	texi2dvi --pdf $<
