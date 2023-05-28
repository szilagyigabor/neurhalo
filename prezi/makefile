all: beamer

beamer: beamer.tex
	pdflatex -interaction=nonstopmode -halt-on-error beamer.tex > output.txt
	pdflatex -interaction=nonstopmode -halt-on-error beamer.tex >> output.txt
	rm -f beamer.aux beamer.bbl beamer.blg beamer.out beamer.log beamer.toc beamer.snm beamer.nav beamer.xmpdata pdfa.xmpi texput.log

