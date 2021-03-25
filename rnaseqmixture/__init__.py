# Copyright (c) 2015, Florian Finkernagel, Philipps University Marburg, Germany.
#  All rights reserved.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""
No RNAseq of primary cells is pure - there are always other cell types sequenced within your sample.
Given your purest sample(s) (=reference) and samples consisting just of the likely contaminations, this module
attempts to correct for the contamination.

How much contamination is present is assessed based on a number of genes that are specific
for your contamination - they should not show any expression in the pure reference.

You first need to define a set of marker genes with find_good_genes_top_k_skip(),
and then correct your samples with correct().

Both take a Pandas DataFrame with samples in columns, genes in rows,
all values in Transcripts Per Million (TPM).
"""


from correction import correct, calculate_percentage, correct_with_percentage, find_good_genes_top_k, find_good_genes_top_k_skip, find_good_genes_threeway_top_k_skip, calculate_percentage_threeway


all = [
    correct, correct_with_percentage,
    calculate_percentage, calculate_percentage_threeway,
    find_good_genes_top_k, find_good_genes_top_k_skip, find_good_genes_threeway_top_k_skip,

]

__version__ = '0.1.0'