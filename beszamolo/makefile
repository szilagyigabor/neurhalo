all: document

document: beszamolo.tex
	pdflatex -interaction=nonstopmode -halt-on-error beszamolo.tex > output.txt
	pdflatex -interaction=nonstopmode -halt-on-error beszamolo.tex >> output.txt
	rm -f beszamolo.aux beszamolo.bbl beszamolo.blg beszamolo.out beszamolo.log beszamolo.toc beszamolo.snm beszamolo.nav beszamolo.xmpdata pdfa.xmpi texput.log

