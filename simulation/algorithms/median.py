from ...correction import find_good_genes_top_k, find_good_genes_top_k_skip, correct, find_good_genes_threeway_top_k_skip, correct_threeway
import numpy as np


def GeneChoice_NGenes(no_of_genes):
    """Chose the top N genes"""
    def select_genes(df_tpm, reference_information):
        genes_for_estimation = ([x for x in find_good_genes_top_k(
            df_tpm, reference_information['sample'], reference_information['contamination'], no_of_genes)])
        return genes_for_estimation
    return 'NGenes_%i' % no_of_genes, select_genes


def GeneChoice_NGenesSkipTopX(no_of_genes, skip_x_genes):
    """Choose the top N genes, skipping the first X"""
    def select_genes(df_tpm, reference_information):
        genes_for_estimation = ([x for x in find_good_genes_top_k_skip(
            df_tpm, reference_information['sample'], reference_information['contamination'], no_of_genes, skip_x_genes)])
        return genes_for_estimation
    return 'NGenesSkipTopX_%i_%i' % (no_of_genes, skip_x_genes), select_genes


def Percentage_Median():
    """Our method of calculating the contamination percentage"""
    def calc(df_tpm, genes_for_estimation):
        try:
            calc_percentage, calc_percentage_std, corrected = correct(
                df_tpm,
                'observed',
                ['reference_sample'],
                ['reference_contamination'],
                genes_for_estimation)
        except ValueError:
            calc_percentage = np.nan
            calc_percentage_std = np.nan
            corrected = df_tpm['observed']
        return {
            'percentage': calc_percentage,
            'percentage_std': calc_percentage_std,
            'corrected': corrected,
        }
    return 'Median', calc


def GeneChoice_Threeway(no_of_genes, skip_x_genes):
    def select_genes(df_tpm, reference_information):
        return find_good_genes_threeway_top_k_skip(df_tpm,
                                                   reference_information['sample'],
                                                   reference_information['contamination1'],
                                                   reference_information['contamination2'],
                                                   no_of_genes, skip_x_genes,
                                                   )
    return 'NGenesSkipTopXThreeway_%i_%i' % (no_of_genes, skip_x_genes), select_genes


def Percentage_Median_Threeway():
    def calc(df_tpm, genes_for_estimation):
        try:
            calc_percentages, corrected = correct_threeway(
                    df_tpm,
                    'observed',
                    ['reference_sample'],
                    ['reference_contamination1'],
                    ['reference_contamination2'],
                    genes_for_estimation[0],
                    genes_for_estimation[1])
        except ValueError:
            calc_percentages = [np.nan, np.nan]
            corrected = df_tpm['observed']
        return {
            'percentage': calc_percentages,
            'corrected': corrected,
        }

    return 'Median3w', calc
