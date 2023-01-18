.PHONY: all clean distclean

all: lab-07.pdf lab-07-solution.pdf

lab-07.pdf: lab-07.tex lab-07-common.tex ../../common/preamble.tex
	latexmk -pdf $<

lab-07-solution.pdf: lab-07-solution.tex lab-07-common.tex ../../common/preamble.tex ../../common/preamble-solution.tex
	latexmk -pdf $<

clean:
	rm -f *.aux *.fls *.tdo *.blg *.bbl *.out *.log *.fdb_latexmk

distclean: clean
	rm *.pdf
