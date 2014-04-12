sa1_az325.pdf: sa1_az325.Rnw
	R --vanilla -e "library(knitr); knit2pdf('sa1_az325.Rnw');"
	pdflatex sa1_az325.tex
	bibtex sa1_az325.aux
	pdflatex sa1_az325.tex

clean:
	rm -f sa1_az325.aux sa1_az325.log sa1_az325.out sa1_az325.tex
	rm -rf figure
	rm -f *~
