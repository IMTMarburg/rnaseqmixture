import DeconRNASeq
import rpy2.robjects
import rpy2.robjects.pandas2ri
import numpy as np
import pyggplot as gg
from ...correction import correct_with_percentage
import pandas as pd


def Percentage_DeconRNASeq():
    DeconRNASeq.load_r()

    def calc(df_tpm, genes_for_estimation):
        if isinstance(genes_for_estimation, list):
            signature_df = df_tpm.ix[genes_for_estimation, [
                'reference_sample', 'reference_contamination']]
        else:
            signature_df = genes_for_estimation
        if not isinstance(signature_df, rpy2.robjects.vectors.DataFrame):
            signature_df = gg.convert_dataframe_to_r(signature_df, True)
        if signature_df is None:
            print genes_for_estimation
            print type(genes_for_estimation)
            raise ValueError("Here it is")

        mixture_df = df_tpm[['observed']]
        # for some arcane reason it won't work if you run it with one sample
        # only...
        mixture_df.insert(1, 'duplicate', mixture_df['observed'])
        mixture_df = gg.convert_dataframe_to_r(mixture_df, True)
        result = rpy2.robjects.r('DeconRNASeq')(
            mixture_df,
            signature_df,
        )
        calc_percentage = result[0][1]  # the contamination percentage
        calc_percentage_std = np.nan
        percentage, corrected = correct_with_percentage(
            df_tpm,
            'observed',
            ['reference_sample'],
            ['reference_contamination'],
            calc_percentage)
        return {
            'percentage': calc_percentage,
            'percentage_std': calc_percentage_std,
            'corrected': corrected,
        }
    return 'DeconRNAseq', calc


def GeneChoice_TopAndBottomNSkip_X(no_of_genes, skip_x_genes):
    """Choose the most regulated genes from up and down (fold change, no t-test, no min expression selection"""
    def select_genes(df_tpm, reference_information):
        df = df_tpm
        reference_mean = df[reference_information['sample']].max(
            axis=1)  # worst case here
        contamination_mean = df[reference_information['contamination']].max(
            axis=1)  # best case here
        contamination_std = df[reference_information['contamination']].std(axis=1)
        contamination_std[np.isnan(contamination_std)] = 0
        enrichment = (contamination_mean + 1) / (reference_mean + 1)

        ok = ((
            (enrichment > 3) |
            (enrichment < 1. / 3)
        ) & (
            (reference_mean > 10) |
            (contamination_mean > 10)
        )
        )

        if not ok.any():
            raise ValueError("none passed filters")

        selected = df[ok]
        selected.insert(0, 'enrichment', enrichment[ok])
        selected = selected.sort('enrichment', ascending=False)
        result = list(selected.index[skip_x_genes: skip_x_genes + no_of_genes])
        result.extend(
            selected.index[-1 * (skip_x_genes + no_of_genes): -1 * skip_x_genes])
        return result
    return 'TopAndBottomNSkip_X_%i_%i' % (no_of_genes, skip_x_genes), select_genes


def GSE19830FixedGenes():
    """Fixed set of genes from the deconrnaseq example"""
    def select_genes(df_tpm, reference_information):
        DeconRNASeq.load_r()
        rpy2.robjects.r("data(rat_liver_brain)")
        x = rpy2.robjects.r('array.signatures')
        return list(rpy2.robjects.r('rownames')(x))
    return 'GSE19830FixedGenes', select_genes


def GSE19830FixedSet():
    """Iterate over the fixed set of mixed samples from GSE19830 / the DeconRNASeq paper"""
    def generate_sample(sample_info, sample_data_raw, list_of_sample_sources, ii, target_percentage):
        samples = [
            ('d15285bm', 0.75,),
            ('d15286bm', 0.75,),
            ('d15287bm', 0.75,),
            ('d15288bm', 0.50,),
            ('d15289bm', 0.50,),
            ('d15290bm', 0.50,),
            ('d15291bm', 0.25,),
            ('d15292bm', 0.25,),
            ('d15293bm', 0.25,),
        ]
        reference_samples = ['d15282bm', 'd15283bm', 'd15284bm']
        reference_contaminations = ['d15294bm', 'd15295bm', 'd15296bm']

        chosen = ii
        if chosen >= len(samples):
            raise ValueError('skipped')

        sub_df = pd.DataFrame(
            {'observed': sample_data_raw[samples[chosen][0]]})
        if ii == (0, 3, 6):
            sub_df.insert(0, 'sample', sample_data_raw['d15282bm'])
            sub_df.insert(0, 'contamination', sample_data_raw['d15294bm'])
        elif ii == (1, 4, 7):
            sub_df.insert(0, 'sample', sample_data_raw['d15283bm'])
            sub_df.insert(0, 'contamination', sample_data_raw['d15295bm'])
        else:
            sub_df.insert(0, 'sample', sample_data_raw['d15284bm'])
            sub_df.insert(0, 'contamination', sample_data_raw['d15296bm'])
        actual_p = samples[chosen][1]
        # can't do that, they mixed the actual samples :).
        # if scipy.stats.pearsonr(sub_df['sample'] * actual_p + sub_df['contamination'] * (1 - actual_p), sub_df['observed'])[0] < 0.9999:
        #     raise ValueError("Sample mixing is not correctly replicated")
        mixture_information = {
                'sample': ['d15282bm'],
                'contamination': ['d15294bm']
        }
        reference_information = {
                'sample': reference_samples,
                'contamination': reference_contaminations
        }

        return sub_df, actual_p, mixture_information, reference_information

    return 'FixedGSE19830', generate_sample, 1


def LiverKidneyFixedSet():
    """Iterate over the fixed set of mixed samples from the 'liver_kidney' data from the package / supplementare fig 1 (deconrnaseq)"""
    def generate_sample(sample_info, sample_data_raw, list_of_sample_sources, ii, target_percentage):
        samples = [
            ('reads.1', 0.33333,),
            ('reads.2', 0.666667,),
            ('reads.3', 0.5,),
            ('reads.4', 0.25,),
            ('reads.5', 0.75,),
            ('reads.6', 0.1,),
            ('reads.7', 0.9,),
        ]
        # all the reference and sample stuff is fake - can't do MAE on this
        # dataset
        reference_samples = ['reads.1']
        reference_contaminations = ['reads.2']

        chosen = ii
        if chosen >= len(samples):
            raise ValueError('skipped')

        sub_df = pd.DataFrame(
            {'observed': sample_data_raw[samples[chosen][0]]})
        sub_df.insert(0, 'sample', [np.nan] * len(sample_data_raw))
        sub_df.insert(0, 'contamination', [np.nan] * len(sample_data_raw))
        actual_p = samples[chosen][1]

        mixture_information = {
                'sample': ['reads.1'],
                'contamination': ['reads.2']
        }
        reference_information = {
                'sample': reference_samples,
                'contamination': reference_contaminations
        }

        return sub_df, actual_p, mixture_information, reference_information

    return 'LiverKidney', generate_sample, 1


def LiverKidneyFixedGenes():
    def select_genes(df_tpm, reference_information):
        DeconRNASeq.load_r()
        rpy2.robjects.r("data(liver_kidney)")
        x = rpy2.robjects.r('signatures')
        return x
    return 'LiverKidneyFixedGenes', select_genes


# question: Does DeconRNASeq work with only differential genes in one direction
def LiverKidneyKidneyOnlyFixedGenes():
    def select_genes(df_tpm, reference_information):
        DeconRNASeq.load_r()
        rpy2.robjects.r("data(liver_kidney)")
        x = rpy2.robjects.r('signatures')
        xdf = rpy2.robjects.pandas2ri.ri2py(x)
        xdf.index = rpy2.robjects.r("rownames")(x)
        xdf = xdf[xdf.kidney > xdf.liver]
        return xdf  # gg.convert_dataframe_to_r(x, True)
    return 'LiverKidneyKidneyOnlyFixedGenes', select_genes


def GSE19830Source():
    """A source for the example dataset of DeconRNASeq"""
    def load():
        DeconRNASeq.load_r()
        rpy2.robjects.r("data(rat_liver_brain)")
        x = rpy2.robjects.r('all.datasets')
        df = rpy2.robjects.pandas2ri.ri2py(x)
        df.index = rpy2.robjects.r('rownames')(x)

        sample_info = {
            'd15282bm': 'liver',
            'd15283bm': 'liver',
            'd15284bm': 'liver',
            'd15285bm': 'mix',
            'd15286bm': 'mix',
            'd15287bm': 'mix',
            'd15288bm': 'mix',
            'd15289bm': 'mix',
            'd15290bm': 'mix',
            'd15291bm': 'mix',
            'd15292bm': 'mix',
            'd15293bm': 'mix',
            'd15294bm': 'brain',
            'd15295bm': 'brain',
            'd15296bm': 'brain',
        }
        return sample_info, df
    return 'GSE19830', load, []


def LiverKidneySource():
    """A source for the example dataset of DeconRNASeq"""
    def load():
        DeconRNASeq.load_r()
        rpy2.robjects.r("data(liver_kidney)")
        x = rpy2.robjects.r('datasets')
        df = rpy2.robjects.pandas2ri.ri2py(x)
        df.index = rpy2.robjects.r('rownames')(x)

        sample_info = {}
        for x in df.columns:
            sample_info[x] = np.nan
        return sample_info, df
    return 'LiverKidney', load, []
