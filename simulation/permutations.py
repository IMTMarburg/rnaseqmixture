import random
import pandas as pd


def PermutationMonteCarlo(no_of_reference_samples=1):
    """Pick one reference sample for contamination and target each,
    choose (different) samples to mix, build mixture"""

    def generate_sample(sample_info, sample_data_raw, list_of_sample_sources, list_of_contamination_sources, ii, target_percentage):
        # what samples are available and in the list of sample sources?
        samples_available = [k for (k, v) in sample_info.items(
        ) if v in list_of_sample_sources and k in sample_data_raw.columns]
        contamination_available = [k for (k, v) in sample_info.items(
        ) if v in list_of_contamination_sources and k in sample_data_raw.columns]

        # decide on which samples to use for the mixture
        if not samples_available:
            raise ValueError("no samples founde")
        sample = random.choice(samples_available)
        contamination = random.choice(contamination_available)
        sample_tissue = sample_info[sample]
        contamination_tissue = sample_info[contamination]
        # if you explicitly order identical tissues for sample and
        # contamination, we accept that
        if list_of_sample_sources[0] != list_of_sample_sources[1]:
            # make sure, we don't select the same tissue for reference and
            # contamination
            while contamination == sample or contamination_tissue == sample_tissue:
                contamination = random.choice(samples_available)
                contamination_tissue = sample_info[contamination]

        # select which samples we'll use for the reference - excluding the
        # reference
        foreground_samples_available = [
            x for x in samples_available if sample_info[x] == sample_tissue and x != sample]
        contamination_samples_available = [x for x in samples_available if sample_info[
            x] == contamination_tissue and x != contamination]
        reference_samples = random.sample(
            foreground_samples_available, no_of_reference_samples)
        reference_contaminations = random.sample(
            contamination_samples_available, no_of_reference_samples)

        # create dataframe with known column names and combine reference
        # samples
        # we need the raw data for the mixture
        sub_df = sample_data_raw[[sample, contamination]]
        sub_df.columns = ['sample', 'contamination']
        # simulate a mixture
        count_sample = sub_df['sample'].sum()
        count_contamination = sub_df['contamination'].sum()
        # we simulate a lane with the same sequencing depth as the smaller lane
        total = min(count_sample, count_contamination)
        from_sample = (1 - target_percentage) * total
        from_contamination = target_percentage * total
        # how many do we need to pick from each condition?
        ratio_sample = from_sample / count_sample
        ratio_contamination = from_contamination / count_contamination
        # we truncate to integers - we simulate complete reads
        expression_sample = (sub_df['sample'] * ratio_sample).astype(int)
        expression_contamination = (
            sub_df['contamination'] * ratio_contamination).astype(int)
        combined_expression = expression_sample + expression_contamination
        actual_p = float(expression_contamination.sum()) / \
            (expression_sample.sum() + expression_contamination.sum())
        sub_df.insert(0, 'observed', combined_expression)
        mixture_information = {
                'sample': [sample],
                'contamination': [contamination]
        }
        reference_information = {
                'sample': reference_samples,
                'contamination': reference_contaminations
        }
        return sub_df, actual_p, mixture_information, reference_information

    return 'permutate-%i' % no_of_reference_samples, generate_sample, 1


def PermutationMonteCarloMultipleSamples(no_of_mixed_samples=2, no_of_reference_samples=1):
    """Pick one reference sample for contamination and target each,
    (repeatedly) choose (different) samples to mix, build mixtures.
    After the second sample, mixture percentages alter between 1-target_percentage and target_percentage

    """

    def generate_sample(sample_info, sample_data_raw, list_of_sample_sources, list_of_contamination_sources, ii, target_percentage):
        # what samples are available and in the list of sample sources?
        samples_available = [k for (k, v) in sample_info.items(
        ) if v in list_of_sample_sources and k in sample_data_raw.columns]
        contamination_available = [k for (k, v) in sample_info.items(
        ) if v in list_of_contamination_sources and k in sample_data_raw.columns]

        # decide on which samples to use for the mixture
        if not samples_available:
            raise ValueError("no samples founde")
        # pick two tissues
        sample = random.choice(samples_available)
        contamination = random.choice(contamination_available)
        sample_tissue = sample_info[sample]
        contamination_tissue = sample_info[contamination]
        # if you explicitly order identical tissues for sample and
        # contamination, we accept that
        if list_of_sample_sources[0] != list_of_sample_sources[1]:
            # make sure, we don't select the same tissue for reference and
            # contamination
            while contamination == sample or contamination_tissue == sample_tissue:
                contamination = random.choice(samples_available)
                contamination_tissue = sample_info[contamination]

        # now pick no_of_mixed_samples + no_of_reference_samples from each
        foreground_samples_available = [
            x for x in samples_available if sample_info[x] == sample_tissue and x != sample]
        contamination_samples_available = [x for x in samples_available if sample_info[
            x] == contamination_tissue and x != contamination]
        samples = random.sample(
            foreground_samples_available, no_of_mixed_samples)
        reference_samples = random.sample(
            foreground_samples_available, no_of_reference_samples)
        contaminations = random.sample(
            contamination_samples_available, no_of_mixed_samples)
        reference_contaminations = random.sample(
            contamination_samples_available, no_of_reference_samples)

        # create dataframe with known column names and combine reference
        # samples
        # we need the raw data for the mixture
        results = {}
        for ii in xrange(0, no_of_mixed_samples):
            if ii > 0:
                target_percentage = 1 - target_percentage  # alternate for undo
            sub_df = sample_data_raw[[samples[ii], contaminations[ii]]]
            sub_df.columns = ['sample', 'contamination']
            # simulate a mixture
            count_sample = sub_df['sample'].sum()
            count_contamination = sub_df['contamination'].sum()
            # we simulate a lane with the same sequencing depth as the smaller
            # lane
            total = min(count_sample, count_contamination)
            from_sample = (1 - target_percentage) * total
            from_contamination = target_percentage * total
            # how many do we need to pick from each condition?
            ratio_sample = from_sample / count_sample
            ratio_contamination = from_contamination / count_contamination
            # we truncate to integers - we simulate complete reads
            expression_sample = (sub_df['sample'] * ratio_sample).astype(int)
            expression_contamination = (
                sub_df['contamination'] * ratio_contamination).astype(int)
            combined_expression = expression_sample + expression_contamination
            if ii == 0:
                actual_p = float(expression_contamination.sum()) / \
                    (expression_sample.sum() + expression_contamination.sum())
            if not results:
                results['index'] = sub_df.index
            if ii == 0:
                key = ''
            else:
                key = '_%i' % ii
            results['observed' + key] = combined_expression
            results['sample' + key] = sample_data_raw[samples[ii]]
            results['contamination' + key] = sample_data_raw[contaminations[ii]]
        result = pd.DataFrame(results).set_index('index')

        mixture_information = {
                'sample': samples,
                'contamination': contaminations
        }
        reference_information = {
                'sample': reference_samples,
                'contamination': reference_contaminations
        }
        return result, actual_p, mixture_information, reference_information

    return 'permutate-ms-%i-%i' % (no_of_mixed_samples, no_of_reference_samples), generate_sample, 1


def PermutationThreewayMonteCarlo(no_of_reference_samples=1):
    """Pick one reference sample for each of (target tissue, contamination 1, contamination 2)
    choose (different) samples to mix, build mixture"""

    def generate_sample(sample_info, sample_data_raw, list_of_sample_sources, list_of_contamination_sources, ii, target_percentage):
        """ @list_of_sample_sources is a list of tissues we accept for this experiment"""
        if len(list_of_sample_sources) != len(set(list_of_sample_sources)):
            raise ValueError("Threeway does not accept duplicates in sample sources""")
        # what samples are available and in the list of sample sources?
        samples_available = [k for (k, v) in sample_info.items(
        ) if v in list_of_sample_sources and k in sample_data_raw.columns]
        contamination_available = [k for (k, v) in sample_info.items(
        ) if v in list_of_contamination_sources and k in sample_data_raw.columns]

        # decide on which samples to use for the mixture
        if not samples_available:
            raise ValueError("no samples founde")
        sample_tissue = None
        contamination1_tissue = None
        contamination2_tissue = None
        while sample_tissue == contamination1_tissue or contamination1_tissue == contamination2_tissue:
            sample = random.choice(samples_available)
            contamination1 = random.choice(contamination_available)
            contamination2 = random.choice(contamination_available)
            sample_tissue = sample_info[sample]
            contamination1_tissue = sample_info[contamination1]
            contamination2_tissue = sample_info[contamination2]

        # select which samples we'll use for the reference - excluding the
        # reference
        foreground_samples_available = [
            x for x in samples_available if sample_info[x] == sample_tissue and x != sample]
        contamination1_samples_available = [x for x in samples_available if sample_info[
            x] == contamination1_tissue and x != contamination1]
        contamination2_samples_available = [x for x in samples_available if sample_info[
            x] == contamination2_tissue and x != contamination2]

        reference_samples = random.sample(
            foreground_samples_available, no_of_reference_samples)
        reference_contaminations1 = random.sample(
            contamination1_samples_available, no_of_reference_samples)
        reference_contaminations2 = random.sample(
            contamination2_samples_available, no_of_reference_samples)

        # create dataframe with known column names and combine reference
        # samples
        # we need the raw data for the mixture
        sub_df = sample_data_raw[[sample, contamination1, contamination2]]
        sub_df.columns = ['sample', 'contamination1', 'contamination2']
        # simulate a mixture
        count_sample = sub_df['sample'].sum()
        count_contamination1 = sub_df['contamination1'].sum()
        count_contamination2 = sub_df['contamination2'].sum()
        # we simulate a lane with the same sequencing depth as the smallest lane

        total = min(count_sample, count_contamination1, count_contamination2)
        from_sample = (1 - target_percentage[0] - target_percentage[1]) * total
        from_contamination1 = target_percentage[0] * total
        from_contamination2 = target_percentage[1] * total
        # how many do we need to pick from each condition?
        ratio_sample = from_sample / count_sample
        ratio_contamination1 = from_contamination1 / count_contamination1
        ratio_contamination2 = from_contamination2 / count_contamination2
        # we truncate to integers - we simulate complete reads
        expression_sample = (sub_df['sample'] * ratio_sample).astype(int)
        expression_contamination1 = (
            sub_df['contamination1'] * ratio_contamination1).astype(int)
        expression_contamination2 = (
            sub_df['contamination2'] * ratio_contamination2).astype(int)

        combined_expression = expression_sample + expression_contamination1 + expression_contamination2
        actual_p = [
                float(expression_contamination1.sum()) / (expression_sample.sum() + expression_contamination1.sum() + expression_contamination2.sum()),
                float(expression_contamination2.sum()) / (expression_sample.sum() + expression_contamination1.sum() + expression_contamination2.sum()),
        ]
        sub_df.insert(0, 'observed', combined_expression)
        mixture_information = {
                'sample': [sample],
                'contamination1': [contamination1],
                'contamination2': [contamination2]
        }
        reference_information = {
                'sample': reference_samples,
                'contamination1': reference_contaminations1,
                'contamination2': reference_contaminations2
        }
        return sub_df, actual_p, mixture_information, reference_information

    return 'perm_three-%i' % no_of_reference_samples, generate_sample, 2


