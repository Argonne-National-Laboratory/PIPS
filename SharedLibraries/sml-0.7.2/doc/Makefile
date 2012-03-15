pdflatex=pdflatex
bibtex=bibtex
xpdf=acroread

all: userguide.pdf

userguide.pdf: userguide.tex userguide.bib
	$(pdflatex) userguide.tex
	$(bibtex) userguide
	$(pdflatex) userguide.tex
	$(pdflatex) userguide.tex

view: userguide.pdf
	$(xpdf) userguide.pdf

clean:
	-rm userguide.blg userguide.aux userguide.log userguide.bbl userguide.pdf 2> /dev/null
