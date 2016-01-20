all: betternorm.pdf AdditionalFile1.pdf Figure1.pdf Figure2.pdf Figure3.pdf Figure4.pdf Figure5.pdf AdditionalFile2.tsv AdditionalFile3.tsv

clean:
	rm *.pdf
	rm *.tsv

betternorm.pdf: betternorm.tex references.bib 
	pdflatex betternorm.tex
	bibtex betternorm.aux
	pdflatex betternorm.tex
	pdflatex betternorm.tex

AdditionalFile1.pdf: betternorm.pdf suppmat.tex ../simulations/results/sizeP_*.pdf pics/library_ring.pdf ../realdata/Zeisel_*.pdf
	pdflatex suppmat.tex 
	mv suppmat.pdf $@

pics/library_ring.pdf: pics/library_ring.svg
	inkscape -D -A $@ $<

Figure1.pdf: pics/figure_existing_sim.tex ../simulations/results/size_* ../simulations/results/lib_*.pdf ../simulations/results/TMM_*.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_existing_sim.pdf $@

Figure2.pdf: pics/figure_zerosim.tex pics/small_cell.pdf pics/large_cell.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_zerosim.pdf $@

pics/small_cell.pdf: pics/zero_effect.R
	cd pics; Rscript zero_effect.R

Figure3.pdf: pics/deconvolution.svg
	inkscape -D -A $@ $<

Figure4.pdf: pics/figure_deconv_sim.tex ../simulations/results/sumClust_*.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_deconv_sim.pdf $@

Figure5.pdf: pics/figure_realdata.tex ../realdata/Zeisel_*.pdf ../realdata/Klein_*.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_realdata.pdf $@

AdditionalFile2.tsv: ../realdata/GOresults/GOoutputUQinDecon_BP.tsv
	cp $< $@

AdditionalFile3.tsv: ../realdata/GOresults/GOoutputUQinLib_BP.tsv
	cp $< $@
