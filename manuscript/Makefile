all: betternorm.pdf suppmat.pdf Figure1.pdf Figure2.pdf Figure3.pdf Figure4.pdf Figure5.pdf

betternorm.pdf: betternorm.tex references.bib 
	pdflatex betternorm.tex
	bibtex betternorm.aux
	pdflatex betternorm.tex
	pdflatex betternorm.tex

suppmat.pdf: betternorm.pdf suppmat.tex ../simulations/results/sizeP_*.pdf pics/library_ring.pdf ../realdata/Zeisel_*.pdf
	pdflatex suppmat.tex

pics/library_ring.pdf: pics/library_ring.svg
	inkscape -D -A $@ $<

Figure1.pdf: figure_existing_sim.tex ../simulations/results/size_* ../simulations/results/lib_*.pdf ../simulations/results/TMM_*.pdf
	pdflatex figure_existing_sim.tex
	pdfcrop figure_existing_sim.pdf $@

Figure2.pdf: figure_zerosim.tex pics/small_cell.pdf pics/large_cell.pdf
	pdflatex figure_zerosim.tex
	pdfcrop figure_zerosim.pdf $@

pics/small_cell.pdf: pics/zero_effect.R
	cd pics; Rscript zero_effect.R

Figure3.pdf: pics/deconvolution.svg
	inkscape -D -A $@ $<

Figure4.pdf: figure_deconv_sim.tex ../simulations/results/sumClust_*.pdf
	pdflatex figure_deconv_sim.tex
	pdfcrop figure_deconv_sim.pdf $@

Figure5.pdf: figure_realdata.tex ../realdata/Zeisel_*.pdf ../realdata/Klein_*.pdf
	pdflatex figure_realdata.tex
	pdfcrop figure_realdata.pdf $@
