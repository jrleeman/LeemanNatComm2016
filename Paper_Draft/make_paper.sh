#!/bin/tcsh

# Remove all of the old aux type files
rm Leeman_et_al.aux
rm Leeman_et_al.bbl
rm Leeman_et_al.dvi
rm Leeman_et_al.log
rm Leeman_et_al.pdf
rm Leeman_et_al.blg

# Make the paper, references, etc
pdflatex Leeman_et_al.tex
bibtex Leeman_et_al.aux
pdflatex Leeman_et_al.tex
pdflatex Leeman_et_al.tex
