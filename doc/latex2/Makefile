main.dvi : main.aux
	latex main

main.aux : main.tex
	latex main

html : dummy
	rm -rf ../html/*
	latex2html -antialias_text -no_math -html_version 3.2,math \
           -show_section_numbers \
           -dir ../html -split +0 -toc_depth 3 -no_navigation main

clean :
	rm -rf html/* *.log *~

dummy :




