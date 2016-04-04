import rpy2.robjects
import pyggplot as gg
import numpy as np
from ...correction import correct_with_percentage


def Percentage_UNDO():
    import UNDO
    UNDO.load_r()

    def calc(df_tpm, genes_for_estimation):
        rpy2.robjects.r("""call_undo = function(x) {
            library("Biobase")
            y = as.matrix(x)
            colnames(y) = c(1,2)
            rownames(y) = rownames(x)
            a = unlist(gene_expression_input(y))
            markergene <- marker_gene_selection(a,lowper=0.4,highper=0.1,epsilon1=0.01,epsilon2=0.01)
            two_source_deconv(y,lowper=0.4,highper=0.1,epsilon1=0.01, epsilon2=0.01,return=0)
            }
            """
                        )
        try:
            query_df = df_tpm[['observed', 'observed_1']].copy()
            # UNDO can't handle zeros in the dataset...
            query_df['observed'] += 1
            query_df['observed_1'] += 1

            x = rpy2.robjects.r('call_undo')(
                gg.convert_dataframe_to_r(query_df, True),
            )
            calc_percentage = np.array(x)[0][0]
            percentage, corrected = correct_with_percentage(
                df_tpm,
                'observed',
                ['reference_sample'],
                ['reference_contamination'],
                calc_percentage)
        except rpy2.rinterface.RRuntimeError as e:
            if 'infinite or missing values in' in str(e):
                calc_percentage = np.nan
                corrected = np.nan
            else:
                raise

        return {
            'percentage': calc_percentage,
            # 'percentage_std': calc_percentage_std,
            'corrected': corrected,
            # 'p-value': p_value
        }
    return 'UNDO', calc


def GeneChoice_Unused():
    def select_genes(df_tpm, reference_information):
        return []
    return 'unused', select_genes
