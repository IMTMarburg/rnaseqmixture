import numpy as np


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
    corrected = (df[sample_to_correct] - cont *
                 percentage) / (1.0 - percentage)
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
    percentage[percentage < 0] = 0
    percentage = percentage.median()
    if percentage > 1:
        raise ValueError("Calculated percentage above 1.0")
    return percentage


def correct_with_percentage(df, sample_to_correct, reference_columns, contamination_columns,
                            percentage):
    corrected = (df[sample_to_correct] - df[contamination_columns].median(axis=1)
                 * percentage) / (1.0 - percentage)
    corrected[corrected < 0] = 0
    return percentage, corrected


def find_good_genes_top_k(df, reference_columns, contamination_columns, k=5):
    """Find decent candidate genes for the correction, bug ignore the most enriched ones
     (they are often not good markers).

    @df Pandas DataFrame
    @reference_columns list of columns to use as pure samples
    @contamination_columns list of columns to use as pure contamination
    @k how many genes to pick
    @skip how many genes to skip
    """
    reference_mean = df[reference_columns].max(axis=1)  # worst case here
    contamination_mean = df[contamination_columns].max(
        axis=1)  # best case here
    contamination_std = df[contamination_columns].std(axis=1)
    contamination_std[np.isnan(contamination_std)] = 0
    enrichment = (contamination_mean + 1) / (reference_mean + 1)

    ok = (
        (reference_mean < 10) &
        # (contamination_mean > min_tpm_threshold) &
        (enrichment > 3)
        # ((contamination_std / contamination_mean).abs() < 0.1)
    )
    if not ok.any():
        raise ValueError("No gene passed the necessary filters")

    selected = df[ok]
    selected.insert(0, 'enrichment', enrichment[ok])
    selected = selected.sort('enrichment', ascending=False)
    # selected = selected.drop(['enrichment'], axis=1)
    return selected.index[:k]


def find_good_genes_top_k_skip(df, reference_columns, contamination_columns, k=5, skip=5):
    """Find decent candidate genes for the correction, bug ignore the very best ones"""
    reference_mean = df[reference_columns].max(axis=1)  # worst case here
    contamination_mean = df[contamination_columns].max(
        axis=1)  # best case here
    contamination_std = df[contamination_columns].std(axis=1)
    contamination_std[np.isnan(contamination_std)] = 0
    enrichment = (contamination_mean + 1) / (reference_mean + 1)

    ok = (
        (reference_mean < 10) &
        # (contamination_mean > min_tpm_threshold) &
        (enrichment > 3)
        # ((contamination_std / contamination_mean).abs() < 0.1)
    )
    if not ok.any():
        raise ValueError("none passed filters")

    selected = df[ok]
    selected.insert(0, 'enrichment', enrichment[ok])
    selected = selected.sort('enrichment', ascending=False)
    return selected.index[skip:skip + k]


def find_good_genes_threeway_top_k_skip(df, reference_columns, contamination1_columns, contamination2_columns, k, skip):
    reference_mean = df[reference_columns].max(axis=1)
    contamination1_mean = df[contamination1_columns].max(axis=1)
    contamination2_mean = df[contamination2_columns].max(axis=1)

    enrichment1 = (contamination1_mean + 1) / (reference_mean + 1)
    enrichment2 = (contamination2_mean + 1) / (reference_mean + 1)

    ok1 = ((reference_mean < 10) & (enrichment1 > 3) & (contamination2_mean < 10))
    ok2 = ((reference_mean < 10) & (enrichment2 > 3) & (contamination1_mean < 10))
    if not ok1.any():
        raise ValueError("None passed filters for contamination 1")
    if not ok2.any():
        raise ValueError("None passed filters for contamination 2")

    selected1 = df[ok1]
    selected1.insert(0, 'enrichment', enrichment1[ok1])
    selected1 = selected1.sort('enrichment', ascending=False)
    marker_genes1 = selected1.index[skip: skip + k]

    selected2 = df[ok2]
    selected2.insert(0, 'enrichment', enrichment2[ok2])
    selected2 = selected2.sort('enrichment', ascending=False)
    marker_genes2 = selected2.index[skip: skip + k]

    return marker_genes1, marker_genes2


def calculate_percentage_threeway(df, sample_to_correct, reference_columns, contamination1_columns, contamination2_columns,
                                  marker_genes1, marker_genes2):
    percentage1 = calculate_percentage(
        df, sample_to_correct, reference_columns, contamination1_columns, marker_genes1)
    percentage2 = calculate_percentage(
        df, sample_to_correct, reference_columns, contamination2_columns, marker_genes2)
    if percentage1 + percentage2 > 1.0:
        raise ValueError("Combined contamination percentage above 1.0 (%.2f, %.2f)" % (
            percentage1, percentage2))
    return percentage1, percentage2


def correct_threeway(df, sample_to_correct, reference_columns, contamination1_columns, contamination2_columns, marker_genes1, marker_genes2):
    print marker_genes1, marker_genes2
    print reference_columns, contamination1_columns, contamination2_columns
    percentage1, percentage2 = calculate_percentage_threeway(
        df, sample_to_correct, reference_columns, contamination1_columns, contamination2_columns, marker_genes1, marker_genes2)
    corrected = (df[sample_to_correct] - df[contamination1_columns].median(axis=1) * percentage1 -
                 df[contamination2_columns].median(axis=1) * percentage2
                 ) / (1.0 - percentage1 - percentage2)
    corrected[corrected < 0] = 0
    return [percentage1, percentage2], corrected
