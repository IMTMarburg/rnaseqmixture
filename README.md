RNAseqmixture
===
A python module to correct (ajust/deconvolute) for contamination of primary cell RNAseq with 
other cell types (eg. Macrophages).

No RNAseq of primary cells is pure - there are always other cell types sequenced within your sample.
Given your purest sample(s) (=reference) and samples consisting just of the likely contaminations, this module
attempts to correct for the contamination.

How much contamination is present is assessed based on a number of genes that are specific
for your contamination - they should not show any expression in the pure reference.


