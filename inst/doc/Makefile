### Makefile (adapted from copula,... original from Douglas Bates)

DRIVERS = Rmpfr-pkg log1mexp-note

#R_PROG = R-devel
R_PROG = R
R_HOME = $(shell ${R_PROG} RHOME)

## "Writing R Extensions" requires: 1st target must also clean up
## NB: This is *BAD* for submitting to JSS (which needs  $(DRIVERS).tex !)
all: All
	${MAKE} clean
## use this for "testing"
All: $(DRIVERS:=.pdf)

$(DRIVERS:=.zip): $(DRIVERS:=.Rnw) $(DRIVERS:=.pdf) $(DRIVERS:=.R) *.bib
	zip $@ $^
	ls -l $@ ; zipinfo $@

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<
%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex *.bib
	${R_PROG} CMD texi2pdf $<
	$(R_HOME)/bin/Rscript -e 'tools::compactPDF("$@", gs_quality = "ebook")'
#	qpdf --stream-data=compress --object-streams=generate $@ temp$@


purge:
	$(RM) Rplots.pdf $(DRIVERS:=.out) *.aux *.log *.blg
clean: purge
	$(RM) $(DRIVERS:=.tex) $(DRIVERS:=.bbl) $(DRIVERS:=-*.pdf)

realclean: clean
	$(RM) $(DRIVERS:=.pdf)
