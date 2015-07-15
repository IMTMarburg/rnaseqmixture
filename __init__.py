## Copyright (c) 2015, Florian Finkernagel, Philipps University Marburg, Germany. 
#  All rights reserved.

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.
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

import pandas as pd
import numpy as np
import scipy
import scipy.optimize
try:
    #optional imports, optional code for genomics pipeline
    import pypipegraph as ppg
    import genomics
    import pydataframe
    make_annotators = True
except ImportError:  
    make_annotators = False
    pass


def correct(df, sample_to_correct, reference_columns, contamination_columns, 
        genes_for_percentage_calculation):
    """
    Correct a column of RNAseq expression by removing contamination expression.
    Needs a Pandas DataFrame @df, containing a column with the transcripts per million for the sample 
    to correct, any number (>=1) of as pure as possible reference and contamination columns,
    and a list of genes to use to judge the percentage of contamination.

    @genes_for_percentage_calculation elements must be in the index of df (can be found via find_good_genes_top_k)

    Returns contamination percentage, corrected Series
    
    """
    genes = df.ix[genes_for_percentage_calculation]
    reference = genes[reference_columns].median(axis=1)
    contamination = genes[contamination_columns].median(axis=1)
    sample = genes[sample_to_correct] 
    percentage = (sample - reference) / contamination
    percentage_std = percentage.std()
    percentage = percentage.median() 
    if percentage > 1:
        raise ValueError("Calculated percentage above 1.0")

    if len(contamination_columns) > 1:
        cont = df[contamination_columns].median(axis=1)
    else:
        cont = df[contamination_columns[0]]
    corrected = (df[sample_to_correct] -  cont * percentage) / (1.0 - percentage)
    corrected[corrected < 0] = 0
    return percentage, percentage_std, corrected


def calculate_percentage(df, sample_to_correct, reference_columns, contamination_columns, 
        genes_for_percentage_calculation):
    """
    Calculate contamination percentage for an RNAseq sample
    Needs a pandas dataframe @df, containing a column with the transcripts per million for the sample 
    to correct, any number (>=1) of as pure as possible reference and contamination columns,
    and a list of genes to use to judge the percentage of contamination.

    @genes_for_percentage_calculation elements must be in the index of df (can be found via find_good_genes_top_k)

    Returns contamination percentage
    
    """
    genes = df.ix[genes_for_percentage_calculation]
    reference = genes[reference_columns].median(axis=1)
    contamination = genes[contamination_columns].median(axis=1)
    sample = genes[sample_to_correct] 
    percentage = (sample - reference) / contamination
    percentage = percentage.median() 
    if percentage > 1:
        raise ValueError("Calculated percentage above 1.0")
    return percentage


def correct_with_percentage(df, sample_to_correct, reference_columns, contamination_columns, 
        percentage):
    corrected = (df[sample_to_correct] - df[contamination_columns].median(axis=1) * percentage) / (1.0 - percentage)
    corrected[corrected < 0] = 0
    return percentage, corrected


def find_good_genes_top_k(df, reference_columns, contamination_columns, k = 5):
    """Find decent candidate genes for the correction"""
    reference_mean = df[reference_columns].max(axis=1) #worst case here
    contamination_mean = df[contamination_columns].max(axis=1) #best case here
    contamination_std = df[contamination_columns].std(axis=1)
    contamination_std[np.isnan(contamination_std)] = 0
    enrichment = (contamination_mean + 1) / (reference_mean + 1)

    ok = (
        (reference_mean < 10)   &
        #(contamination_mean > min_tpm_threshold) & 
        (enrichment > 3) 
        #((contamination_std / contamination_mean).abs() < 0.1)
        )
    if not ok.any():
        raise ValueError("none passed filters")
    
    selected = df[ok]
    selected.insert(0, 'enrichment', enrichment[ok])
    selected = selected.sort('enrichment', ascending=False)
    #selected = selected.drop(['enrichment'], axis=1)
    return selected.index[:k]

def find_good_genes_top_k_skip(df, reference_columns, contamination_columns, k = 5, skip=5):
    """Find decent candidate genes for the correction, bug ignore the very best ones
    (they are often highly enriched, but not good markers).

    @df Pandas DataFrame
    @reference_columns list of columns to use as pure samples
    @contamination_columns list of columns to use as pure contamination
    @k how many genes to pick
    @skip how many genes to skip
    """
    reference_mean = df[reference_columns].max(axis=1) #worst case here
    contamination_mean = df[contamination_columns].max(axis=1) #best case here
    contamination_std = df[contamination_columns].std(axis=1)
    contamination_std[np.isnan(contamination_std)] = 0
    enrichment = (contamination_mean + 1) / (reference_mean + 1)

    ok = (
        (reference_mean < 10)   &
        #(contamination_mean > min_tpm_threshold) & 
        (enrichment > 3) 
        #((contamination_std / contamination_mean).abs() < 0.1)
        )
    if not ok.any():
        raise ValueError("No gene passed the necessary filters")
    
    selected = df[ok]
    selected.insert(0, 'enrichment', enrichment[ok])
    selected = selected.sort('enrichment', ascending=False)
    #selected = selected.drop(['enrichment'], axis=1)
    return selected.index[skip:skip + k]

if make_annotators:  # support code for IMT pipeline - useless otherwise
    class MixtureCorrector(genomics.annotators.Annotator):

        def __init__(self, tpm_to_correct, tpm_a, tpm_b, markers_for_b):
            """@tpm* are annotators.
            Markers for b may be alist, an integer (use topK), or a tuple (k, j) , use top j to j+k"""
            self.tpm_to_correct = tpm_to_correct
            self.tpm_a = tpm_a
            self.tpm_b = tpm_b
            if not hasattr(markers_for_b, '__iter__'):
                markers_for_b = [markers_for_b]
            self.genes_to_consider = markers_for_b
            self.column_names = ['corrected %s' % tpm_to_correct.column_name]
            genomics.annotators.Annotator.__init__(self)

        def annotate(self, genes):
            #corrected_tpm = (tpm_to_correct - tpm_b * percentage /) (1 - percentage)
            #calculate percentage.
            df = pd.DataFrame(genes.df[:, ['stable_id', self.tpm_to_correct.column_name, self.tpm_a.column_name, self.tpm_b.column_name]].value_dict)
            df = df.set_index('stable_id')
            
            marker_genes = set()
            if isinstance(self.genes_to_consider, list):
                for name_or_stable_id in self.genes_to_consider:
                    stable_ids = genes.genome.name_to_gene_ids_hard(name_or_stable_id)
                    if len(stable_ids) > 1:
                        raise ValueError("Multiple stable ids for %s - pick one: %s " % (name_or_stable_id, stable_ids))
                    stable_id = list(stable_ids)[0]
                    marker_genes.update(stable_id)
            elif isinstance(self.genes_to_consider, int):
                marker_genes = find_good_genes_top_k(df, [self.tpm_a.column_name], [self.tpm_b.column_name], self.genes_to_consider)
            elif isinstance(self.genes_to_consider, tuple):
                marker_genes = find_good_genes_top_k_skip(df, [self.tpm_a.column_name], [self.tpm_b.column_name], self.genes_to_consider[0], self.genes_to_consider[1])
            else:
                raise ValueError("could not understand genes_to_consider")
            
            percentage, percentage_std, column = correct(df, self.tpm_to_correct.column_name, [self.tpm_a.column_name], [self.tpm_b.column_name], marker_genes)
            with open(genes.write().job_id + ' %s.txt' % self.name, 'wb') as op:
                op.write("Genes selected: %s\n" % ", ".join(marker_genes))
                p.write("Contamination percentage: %f\n" % percentage)
            return pydataframe.DataFrame({self.column_name: column})

        def get_dependencies(self, genes):
            return [genes.add_annotator(self.tpm_to_correct), genes.add_annotator(self.tpm_a), genes.add_annotator(self.tpm_b),
                    ppg.ParameterInvariant('rnaseqmixture.correct', correct)
                    
                    ]
                

