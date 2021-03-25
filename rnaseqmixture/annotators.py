import pypipegraph as ppg
import genomics
import pydataframe
import pandas as pd
from correction import correct, find_good_genes_top_k, find_good_genes_top_k_skip, find_good_genes_threeway_top_k_skip, correct_threeway


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
        # corrected_tpm = (tpm_to_correct - tpm_b * percentage /) (1 - percentage)
        # calculate percentage.
        df = pd.DataFrame(genes.df[:, ['stable_id', self.tpm_to_correct.column_name,
                                       self.tpm_a.column_name, self.tpm_b.column_name]].value_dict)
        df = df.set_index('stable_id')

        marker_genes = set()
        if isinstance(self.genes_to_consider, list):
            for name_or_stable_id in self.genes_to_consider:
                stable_ids = genes.genome.name_to_gene_ids_hard(
                    name_or_stable_id)
                if len(stable_ids) > 1:
                    raise ValueError(
                        "Multiple stable ids for %s - pick one: %s " % (name_or_stable_id, stable_ids))
                stable_id = list(stable_ids)[0]
                marker_genes.update(stable_id)
        elif isinstance(self.genes_to_consider, int):
            marker_genes = find_good_genes_top_k(df, [self.tpm_a.column_name], [
                                                 self.tpm_b.column_name], self.genes_to_consider)
        elif isinstance(self.genes_to_consider, tuple):
            marker_genes = find_good_genes_top_k_skip(df, [self.tpm_a.column_name], [
                                                      self.tpm_b.column_name], self.genes_to_consider[0], self.genes_to_consider[1])
        else:
            raise ValueError("could not understand genes_to_consider")

        percentage, percentage_std, column = correct(df, self.tpm_to_correct.column_name, [
                                                     self.tpm_a.column_name], [self.tpm_b.column_name], marker_genes)
        with open(genes.write().job_id + ' %s.txt' % self.name, 'wb') as op:
            op.write("Genes selected: %s\n" % ", ".join(marker_genes))
            op.write("Contamination percentage: %f\n" % percentage)
        return pydataframe.DataFrame({self.column_name: column})

    def get_dependencies(self, genes):
        return [genes.add_annotator(self.tpm_to_correct), genes.add_annotator(self.tpm_a), genes.add_annotator(self.tpm_b),
                ppg.ParameterInvariant('rnaseqmixture.correct', correct)

                ]


class MixtureCorrectorThreeway(genomics.annotators.Annotator):

    def __init__(self, tpm_to_correct, reference_tpm_target_tissue, reference_tpm_contamination1, reference_tpm_contamination2,
             top_k, skip_j, filter_to_biotypes = None
             ):
            """@tpm* are annotators.
            Markers for b may be alist, an integer (use topK), or a tuple (k, j) , use top j to j+k.

            If filter_to_biotypes is set, only these genes are considered for the correction
            """
            self.tpm_to_correct = tpm_to_correct
            self.reference_tpm_target_tissue = reference_tpm_target_tissue
            self.reference_tpm_contamination1 = reference_tpm_contamination1
            self.reference_tpm_contamination2 = reference_tpm_contamination2
            self.top_k = top_k
            self.skip_j = skip_j
            self.filter_to_biotypes = set(filter_to_biotypes)
            self.column_names = ['corrected threeway %s' % tpm_to_correct.column_name]
            genomics.annotators.Annotator.__init__(self)

    def annotate(self, genes):
        # corrected_tpm = (tpm_to_correct - tpm_b * percentage /) (1 - percentage)
        # calculate percentage.
        df = pd.DataFrame(genes.df[:, ['stable_id', 'biotype', self.tpm_to_correct.column_name,
                                       self.reference_tpm_target_tissue.column_name,
                                       self.reference_tpm_contamination1.column_name,
                                       self.reference_tpm_contamination2.column_name,
                                       ]].value_dict)
        df = df.set_index('stable_id')
        search_df = df
        if self.filter_to_biotypes is not None:
            search_df = df[df['biotype'].isin(self.filter_to_biotypes)]

        marker_genes = find_good_genes_threeway_top_k_skip(search_df,
               reference_columns = [self.reference_tpm_target_tissue.column_name],
               contamination1_columns = [self.reference_tpm_contamination1.column_name],
               contamination2_columns = [self.reference_tpm_contamination2.column_name],
               k=self.top_k, skip=self.skip_j)

        percentages, column = correct_threeway(df,
               self.tpm_to_correct.column_name,
               reference_columns = [self.reference_tpm_target_tissue.column_name],
               contamination1_columns = [self.reference_tpm_contamination1.column_name],
               contamination2_columns = [self.reference_tpm_contamination2.column_name],
               marker_genes1 = marker_genes[0],
               marker_genes2 = marker_genes[1]
        )
        with open(genes.write().job_id + ' %s.txt' % self.column_name, 'wb') as op:
            op.write("Genes for %s: %s\n" % (self.reference_tpm_contamination1.column_name, ", ".join(marker_genes[0])))
            op.write("Genes for %s: %s\n" % (self.reference_tpm_contamination2.column_name, ", ".join(marker_genes[1])))
            op.write("Contamination percentags: %s\n" % percentages)
        return pydataframe.DataFrame({self.column_name: column})

    def get_dependencies(self, genes):
        return [genes.add_annotator(self.tpm_to_correct),
                genes.add_annotator(self.reference_tpm_target_tissue),
                genes.add_annotator(self.reference_tpm_contamination1),
                genes.add_annotator(self.reference_tpm_contamination2),
                ppg.ParameterInvariant('rnaseqmixture.correct', (correct, self.filter_to_biotypes))
                ]
