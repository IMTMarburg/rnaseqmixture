import os
import cPickle
import pypipegraph as ppg
import numpy as np
import random
import pyggplot as gg
import scipy.stats
import pandas as pd
import collections


def ensure_path(path):
    result = True
    try:
        if not os.path.exists(path):
            os.makedirs(path)
            result = False
    except OSError, e:
        print e
        if str(e).find('File exists') != -1:
            result = False
        else:
            raise e
    if not os.path.exists(path) and os.path.isdir(path):
        raise ValueError("Failed to create %s" % path)
    return result

_source_loading_jobs = {}
_source_data = {}


def load_source(name, callback, dependencies):
    cache_dir = 'cache/mixture_simulation'
    if name not in _source_loading_jobs:
        def calc():
            sample_info, sample_data_raw = callback()
            return sample_info, sample_data_raw

        def load(tup):
            sample_info, df = tup
            factors = 1e6 / df.sum(axis=0)
            df_tpm = df.mul(factors, axis=1)
            _source_data[name] = (sample_info, df, df_tpm)

        _source_loading_jobs[name] = ppg.CachedDataLoadingJob(cache_dir + '/source_%s' % name,
                calc, load).depends_on(dependencies).depends_on(
                        ppg.FunctionInvariant(cache_dir + '/source_func_%s' % name, callback))
    return _source_loading_jobs[name]


class MixtureSimulation:

    def __init__(self, data_loading_strategy, genes_selection_strategy, percentage_calculation_strategy, permutation_strategy):
        self.name = "%s_%s_%s_%s" % (permutation_strategy[0], data_loading_strategy[
                                     0], genes_selection_strategy[0], percentage_calculation_strategy[0])
        self.cache_dir = os.path.join('cache', 'mixture_simulation', self.name)
        self.result_dir = os.path.join(
            'results', 'mixture_simulation', self.name)
        self.data_loading_strategy_name = data_loading_strategy[0]
        self.data_loading_strategy = data_loading_strategy[1]
        self.data_loading_strategy_dependencies = data_loading_strategy[2]
        self.genes_selection_strategy = genes_selection_strategy[1]
        self.percentage_calculation_strategy = percentage_calculation_strategy[
            1]
        self.permutation_strategy = permutation_strategy[1]
        self.percentages_needed_per_permutation = permutation_strategy[2]

    def load(self):
        def load():
            name = self.data_loading_strategy_name
            self.sample_info = _source_data[name][0]
            self.sample_data_raw = _source_data[name][1]
            self.sample_data_tpm = _source_data[name][2]
        return ppg.DataLoadingJob(self.cache_dir + '/load', load).depends_on(load_source(self.data_loading_strategy_name, self.data_loading_strategy, self.data_loading_strategy_dependencies))

    def simulate(self, list_of_sample_sources, no_of_repetitions, list_of_contamination_sources = None, dump_all=False):
        if list_of_contamination_sources is None:
            key = str(list_of_sample_sources)
            list_of_contamination_sources = list_of_sample_sources
        else:
            key = str(list_of_sample_sources) + '-vs-' + str(list_of_contamination_sources)

        ensure_path(os.path.join(self.cache_dir, key))
        np.random.seed(1234 * hash(self.name + key) % 4294967295)
        percentages = self.create_target_percentages(no_of_repetitions)
        jobs = []
        for ii in xrange(no_of_repetitions):
            p = percentages[ii]
            output_filename = os.path.join(
                self.cache_dir, key, '%i.pickle' % ii)

            def calc(output_filename=output_filename, job_id=ii, p=p):
                random.seed((job_id + 30000 * hash(self.name + key)) %
                            4294967295)

                try:
                    sub_df, actual_p, mixture_information, reference_information = self.permutation_strategy(
                        self.sample_info, self.sample_data_raw, list_of_sample_sources, list_of_contamination_sources, job_id, target_percentage=p)
                except ValueError as e:
                    if 'skipped' in e:
                        result = {'skipped': True}
                        with open(output_filename, 'wb') as op:
                            cPickle.dump(result, op, cPickle.HIGHEST_PROTOCOL)
                        return
                    else:
                        raise

                try:
                    genes_for_estimation = self.genes_selection_strategy(
                        self.sample_data_tpm, reference_information)  # a list of stable_ids or gene snames
                except ValueError as e:
                    result = {"error": str(e)}
                    result['mixture_info'] = mixture_information
                    result['reference_info'] = reference_information
                    with open(output_filename, 'wb') as op:
                        cPickle.dump(result, op, cPickle.HIGHEST_PROTOCOL)
                    return
                # add references - median between lanes, based on TPM!
                for group in reference_information:
                    if len(reference_information[group]) > 1:
                        sub_df.insert(0, 'reference_%s' % group, self.sample_data_tpm[
                                      reference_information[group]].median(axis=1))
                    else:
                        sub_df.insert(0, 'reference_%s' % group,
                                      self.sample_data_tpm[reference_information[group][0]])

                sub_df = sub_df.copy()
                # (re)normalize to TPM
                for x in sub_df.columns:
                    if sub_df.dtypes[x] == np.int64:
                        sub_df[x] = sub_df[x].astype(
                            'float') * 1e6 / float(sub_df[x].sum())

                algo_result = self.percentage_calculation_strategy(
                    sub_df, genes_for_estimation)
                calc_percentage = algo_result['percentage']
                corrected = algo_result['corrected']
                del algo_result['corrected']  # don't store that...
                del algo_result['percentage']  # don't store that...
                percentage_ok = True
                if hasattr(calc_percentage, '__iter__'):
                    if np.isnan(calc_percentage).all():
                        percentage_ok = False
                else:
                    if np.isnan(calc_percentage):
                        percentage_ok = False
                if percentage_ok:
                    # sample is the original, observed the mixture value
                    # and corrected is after the algorithm massaged it
                    r_squared_observed = np.sqrt(
                        (sub_df['sample'] - sub_df['observed']) ** 2)
                    r_squared_corrected = np.sqrt(
                        (corrected - sub_df['sample']) ** 2)

                    absolute_error_observed = np.abs(
                        sub_df['sample'] - sub_df['observed'])
                    absolute_error_corrected = np.abs(corrected - sub_df['sample'])
                else:
                    r_squared_observed = np.array([np.nan])
                    r_squared_corrected = np.array([np.nan])
                    absolute_error_observed = np.array([np.nan])
                    absolute_error_corrected = np.array([np.nan])

                result = {
                    'mixture_info': mixture_information,
                    'reference_info': reference_information,
                    'percentage': p,  # that's what should have been mixe
                    'actual percentage': actual_p,  # that's what was mixed
                    'calculated percentage': calc_percentage,  # that's what the algorithm estimated
                    'r_squared_observed_sum': r_squared_observed.sum(),
                    'r_squared_corrected_sum': r_squared_corrected.sum(),
                    'r_squared_observed_mean': r_squared_observed.mean(),
                    'r_squared_corrected_mean': r_squared_corrected.mean(),
                    'absolute_error_observed_mean': absolute_error_observed.mean(),
                    'absolute_error_corrected_mean': absolute_error_corrected.mean(),
                }

                for k in algo_result:
                    result['algo_' + k] = algo_result[k]

                if dump_all:
                    result.update({
                        'r_squared_observed': r_squared_observed,
                        'r_squared_corrected': r_squared_corrected,
                        'marker_genes': genes_for_estimation,
                        'absolute_error_observed': absolute_error_observed,
                        'absolute_error_corrected': absolute_error_corrected,
                    })
                    if isinstance(genes_for_estimation, list):
                        columns = []
                        for group in reference_information:
                            columns.extend(reference_information[group])
                        result['signature_matrix'] = self.sample_data_tpm[
                            columns].ix[genes_for_estimation]
                    else:
                        result['signature_matrix'] = genes_for_estimation
                with open(output_filename, 'wb') as op:
                    cPickle.dump(result, op, cPickle.HIGHEST_PROTOCOL)
            job = ppg.FileGeneratingJob(output_filename, calc).depends_on(
                self.load()).depends_on(ppg.ParameterInvariant(output_filename, dump_all))
            job.ignore_code_changes()
            jobs.append(job)
        if self.percentages_needed_per_permutation == 1:
            self._plot_actual_vs_calculated_percentage(jobs, key)
            self._plot_mae_improvement(jobs, key)
            self._plot_mae_improvement_calc_over_30(jobs, key)
        elif self.percentages_needed_per_permutation > 1:
            self._plot_actual_vs_calculated_percentage_threeway(jobs, key)
        self._plot_sum_vs_delta(jobs, key)
        self._plot_delta_by_contamination(jobs, key)
        self._plot_delta_by_sample(jobs, key)
        self._plot_delta_by_pair(jobs, key)
        self.tabularize_errors(jobs, key)

    def create_target_percentages(self, no_of_repetitions):
        if self.percentages_needed_per_permutation == 1:
            # return np.random.beta(1.5, 5, no_of_repetitions)
            return np.random.uniform(0, 0.5, no_of_repetitions)
        elif self.percentages_needed_per_permutation == 2:
            # p1 = np.random.beta(1.5, 5, no_of_repetitions)
            # p2 = np.random.beta(1.5, 5, no_of_repetitions)
            p1 = np.random.uniform(0, 0.5, no_of_repetitions)
            p2 = np.random.uniform(0, 0.5, no_of_repetitions)
            while (p1 + p2 > 1).any():
                p2[p1 + p2 > 1] /= 2.0
            return zip(p1, p2)
        else:
            raise NotImplementedError()

    def _load_simulation_results(self, simulation_jobs, key):
        def calc():
            simulation_runs = {
                'sample_tissue': [],
                'percentage_ok': [],
                # 'calculated_percentage': [],
                # 'r2_observed_sum': [],
                # 'r2_corrected_sum': [],
                # 'r2_observed_mean': [],
                # 'r2_corrected_mean': [],
                # 'calculated percentage_std': [],
                # 'absolute_error_observed_mean': [],
                # 'absolute_error_corrected_mean': [],
            }
            path = os.path.join(self.cache_dir, key)
            for fn in os.listdir(path):
                if fn.endswith('.pickle'):
                    with open(os.path.join(path, fn)) as op:
                        d = cPickle.load(op)
                        if 'skipped' in d and d['skipped']:
                            continue
                        elif 'error' in d:
                            continue
                        else:
                            simulation_runs['sample_tissue'].append(
                                self.sample_info[d['mixture_info']['sample'][0]])
                            for x in ['contamination', 'contamination1', 'contamination2']:
                                if x in d['mixture_info']:
                                    if not x + '_tissue' in simulation_runs:
                                        simulation_runs[x + '_tissue'] = []
                                    simulation_runs[x + '_tissue'].append(
                                        self.sample_info[d['mixture_info'][x][0]])
                            simulation_runs['percentage_ok'].append(
                                d['calculated percentage'] < .3)
                            for k in d:
                                if k not in ['r_squared_observed', 'r_squared_corrected', 'marker_genes', 'absolute_error_observed',
                                        'absolute_error_corrected']:
                                    if k not in simulation_runs:
                                        simulation_runs[k] = []
                                    simulation_runs[k].append(d[k])
                            # simulation_runs['r2_observed_sum'].append(d['r_squared_observed_sum'])
                            # simulation_runs['r2_corrected_sum'].append(d['r_squared_corrected_sum'])
                            # simulation_runs['r2_observed_mean'].append(d['r_squared_observed_mean'])
                            # simulation_runs['r2_corrected_mean'].append(d['r_squared_corrected_mean'])
                            # simulation_runs['calculated_percentage'].append(d['calculated percentage'])
                            # simulation_runs['actual_percentage'].append(d['actual percentage'])
                            # simulation_runs['calculated percentage_std'].append(d['calculated percentage_std'])
                            # simulation_runs['absolute_error_observed_mean'].append(d['absolute_error_observed_mean'])
                            # simulation_runs['absolute_error_corrected_mean'].append(d['absolute_error_corrected_mean'])
            for k in simulation_runs:
                print k, len(simulation_runs[k])
            simulation_runs = pd.DataFrame(simulation_runs)
            return simulation_runs
        return ppg.CachedAttributeLoadingJob(os.path.join(self.cache_dir, key, 'simulation_results'), self, 'simulation_results', calc).depends_on(simulation_jobs).depends_on(self.load())

    def tabularize_errors(self, simulate_jobs, key):
        of = os.path.join(self.result_dir, key, 'errors.txt')
        ensure_path(os.path.join(self.result_dir, key))

        def write():
            errors = collections.Counter()
            path = os.path.join(self.cache_dir, key)
            for fn in os.listdir(path):
                if fn.endswith('.pickle'):
                    with open(os.path.join(path, fn)) as op:
                        d = cPickle.load(op)
                        if 'error' in d:
                            errors[d['error']] += 1
            with open(of, 'wb') as op:
                op.write("error count:\n")
                for count, error in reversed(sorted([(y, x) for (x, y) in errors.items()])):
                    op.write("%i %s\n" % (count, error))
        return ppg.FileGeneratingJob(of, write).depends_on(simulate_jobs)

    def _plot_actual_vs_calculated_percentage(self, simulate_jobs, key):
        ensure_path(os.path.join(self.result_dir, key))
        output_filename = os.path.join(
            self.result_dir, key, 'actual vs calculated percentage.svg')

        def plot():
            df = self.simulation_results
            df = df[~np.isnan(df['calculated percentage'])]
            p = gg.Plot(df)
            p.theme_bw()
            if 'algo_p-value' in df.columns:
                p.add_scatter('actual percentage', 'calculated percentage',
                              alpha=0.3, size=1, color='algo_p-value')
                p.scale_color_gradient('red', 'blue', name='p-value')
            elif 'algo_calculated_percentage_std' in df.columns:
                p.add_scatter('actual percentage', 'calculated percentage',
                              alpha=0.3, size=1, color='algo_calculated_percentage_std')
                p.scale_color_gradient('red', 'blue', name='% std')
            else:
                p.add_scatter('actual percentage',
                              'calculated percentage', alpha=0.3, size=1)
            p.add_ab_line(0, 1, alpha=0.5)
            r = scipy.stats.linregress(df['actual percentage'], df['calculated percentage'])
            p.add_ab_line(r[1], r[0], color='blue')

            p.add_density_2d()
            p.scale_x_continuous(limits=[0, 1])
            p.scale_y_continuous(limits=[0, 1])
            r = scipy.stats.pearsonr(df['actual percentage'], df[
                                     'calculated percentage'])[0]
            delta = df['actual percentage'] - df['calculated percentage']
            rmse = np.sqrt(np.mean(2 * delta ** 2))
            p.title(self.name + ' ' + key + "\nr=%.2f, rmse=%.2f" % (r, rmse))
            p.render(output_filename, width=11.69, height=8.27)

        return ppg.FileGeneratingJob(output_filename, plot).depends_on(self._load_simulation_results(simulate_jobs, key)).depends_on(self.load())

    def _plot_actual_vs_calculated_percentage_threeway(self, simulate_jobs, key):
        ensure_path(os.path.join(self.result_dir, key))
        output_filename = os.path.join(
            self.result_dir, key, 'actual vs calculated percentage threeway.svg')

        def plot(df):
            df = self.simulation_results
            actual = []
            calculated = []
            no = []
            for values in df['actual percentage']:
                actual.extend(values)
                no.extend(xrange(len(values)))
            for values in df['calculated percentage']:
                calculated.extend(values)
            plot_df = pd.DataFrame({'actual percentage': actual, 'calculated percentage': calculated,
                'contamination number': no})
            plot_df = plot_df[~np.isnan(plot_df['calculated percentage'])]
            plot_df['contamination number'] = pd.Categorical(plot_df['contamination number'])
            p = gg.Plot(plot_df)
            p.theme_bw()
            p.add_scatter('actual percentage',
                          'calculated percentage', alpha=0.3, size=1, color='contamination number')
            p.add_ab_line(0, 1, alpha=0.5)
            p.add_density_2d()
            r = scipy.stats.linregress(plot_df['actual percentage'], plot_df['calculated percentage'])
            p.add_ab_line(r[1], r[0], color='blue', alpha=0.5)
            
            p.scale_x_continuous(limits=[0, 1])
            p.scale_y_continuous(limits=[0, 1])
            p.legend_position('bottom')
            p.scale_color_manual(['red', 'black'])
            r = scipy.stats.pearsonr(plot_df['actual percentage'], plot_df[
                                     'calculated percentage'])[0]
            delta = plot_df['actual percentage'] - plot_df['calculated percentage']
            rmse = np.sqrt(np.mean(2 * delta ** 2))
            p.title(self.name + ' ' + key + "\nr=%.2f, rmse=%.2f" % (r, rmse))
            p.render(output_filename, width=11.69, height=8.27)
            plot_df.to_csv(output_filename + '.tsv', sep="\t")

        return ppg.FileGeneratingJob(output_filename, plot).depends_on(self._load_simulation_results(simulate_jobs, key)).depends_on(self.load())

    def _plot_delta_by_contamination(self, simulate_jobs, key):
        ensure_path(os.path.join(self.result_dir, key))
        output_filename = os.path.join(
            self.result_dir, key, 'delta by contamination sample.svg')

        def plot():
            df = self.simulation_results
            deltas = []
            contaminations = []
            if 'contamination' in df.ix[0, 'mixture_info']:
                deltas.extend((df['actual percentage'] - df['calculated percentage']))
                print df.head()
                contaminations.extend([x['contamination'][0] for x in df['mixture_info']])
            else:
                for idx, row in df.iterrows():
                    deltas.append((row['actual percentage'][0] - row['calculated percentage'][0]))
                    deltas.append((row['actual percentage'][1] - row['calculated percentage'][1]))
                    contaminations.append(row['mixture_info']['contamination1'][0])
                    contaminations.append(row['mixture_info']['contamination2'][0])
            plot_df = pd.DataFrame({"delta": deltas, 'contamination sample': contaminations})
            medians = plot_df.groupby('contamination sample').aggregate(np.min)['delta'].copy()
            medians.sort()
            plot_df['contamination sample'] = pd.Categorical(plot_df['contamination sample'], list(reversed(medians.index)))

            p = gg.Plot(plot_df)
            p.theme_bw()
            # p.add_boxplot('contamination sample', 'delta', )
            p.add_jitter('contamination sample', 'delta', jitter_x=True, jitter_y=False, alpha=0.3, size=1, color='blue')
            p.add_horizontal_bar(0)
            p.turn_x_axis_labels(90)
            p.width = len(medians) * 0.1 + 2
            p.render(output_filename, width=11.69, height=8.27)
        return ppg.FileGeneratingJob(output_filename, plot).depends_on(self._load_simulation_results(simulate_jobs, key)).depends_on(self.load())

    def _plot_delta_by_sample(self, simulate_jobs, key):
        ensure_path(os.path.join(self.result_dir, key))
        output_filename = os.path.join(
            self.result_dir, key, 'delta by reference sample.svg')

        def plot():
            df = self.simulation_results
            deltas = []
            contaminations = []
            if 'contamination' in df.ix[0, 'mixture_info']:
                deltas.extend((df['actual percentage'] - df['calculated percentage']))
                contaminations.extend([x['sample'][0] for x in df['mixture_info']])
            else:
                for idx, row in df.iterrows():
                    deltas.append((row['actual percentage'][0] - row['calculated percentage'][0]))
                    deltas.append((row['actual percentage'][1] - row['calculated percentage'][1]))
                    contaminations.append(row['mixture_info']['sample'][0])
                    contaminations.append(row['mixture_info']['sample'][0])
            plot_df = pd.DataFrame({"delta": deltas, 'reference sample': contaminations})
            medians = plot_df.groupby('reference sample').aggregate(np.min)['delta'].copy()
            medians.sort()
            plot_df['reference sample'] = pd.Categorical(plot_df['reference sample'], list(reversed(medians.index)))

            p = gg.Plot(plot_df)
            p.theme_bw()
            # p.add_boxplot('reference sample', 'delta', )
            p.add_jitter('reference sample', 'delta', jitter_x=True, jitter_y=False, alpha=0.3, size=1, color='blue')
            p.add_horizontal_bar(0)
            p.turn_x_axis_labels(90)
            p.width = len(medians) * 0.1 + 2
            p.render(output_filename, width=11.69, height=8.27)
        return ppg.FileGeneratingJob(output_filename, plot).depends_on(self._load_simulation_results(simulate_jobs, key)).depends_on(self.load())

    def _plot_delta_by_pair(self, simulate_jobs, key):
        ensure_path(os.path.join(self.result_dir, key))
        output_filename = os.path.join(
            self.result_dir, key, 'delta by reference pair.png')

        def plot():
            df = self.simulation_results
            deltas = []
            pairs = []
            for idx, row in df.iterrows():
                sample = row['mixture_info']['sample'][0]
                if 'contamination' in row['mixture_info']:
                    contamination = row['mixture_info']['contamination'][0]
                    deltas.append((row['actual percentage'] - row['calculated percentage']))
                    pairs.append("%s %s" % (sample, contamination))
                else:
                    contamination = row['mixture_info']['contamination1'][0] + ' ' + row['mixture_info']['contamination2'][0]
                    deltas.append((row['actual percentage'][0] - row['calculated percentage'][0]))
                    deltas.append((row['actual percentage'][1] - row['calculated percentage'][1]))
                    pairs.append("%s %s" % (sample, contamination))
                    pairs.append("%s %s" % (sample, contamination))

            plot_df = pd.DataFrame({"delta": deltas, 'pair': pairs})
            if False:
                p = gg.Plot(pd.DataFrame({"x": [0], 'y': [0], 'text': ['no_data']}))
                p.add_text('x', 'y', 'text')
                p.render(output_filename, width=11.69, height=8.27)
                return
            if len(plot_df['pair'].unique()) >= 1:
                medians = plot_df.groupby('pair').aggregate(np.min)['delta'].copy()
                medians.sort()
                plot_df['pair'] = pd.Categorical(plot_df['pair'], list(reversed(medians.index)))
            print plot_df.head()

            p = gg.Plot(plot_df)
            p.theme_bw()
            # p.add_boxplot('reference sample', 'delta', )
            p.add_jitter('pair', 'delta', jitter_x=False, jitter_y=False, alpha=0.3, size=1, color='pair')
            p.scale_color_manual(['red', 'black'] * len(medians))
            p.add_alternating_background('pair')
            p.add_horizontal_bar(0)
            p.turn_x_axis_labels(90)
            p.hide_x_axis_labels()
            p.width = len(medians) * 0.1 + 2
            p.render(output_filename, width=11.69, height=8.27)
        return ppg.FileGeneratingJob(output_filename, plot).depends_on(self._load_simulation_results(simulate_jobs, key)).depends_on(self.load())

    def _plot_sum_vs_delta(self, simulate_jobs, key):
        ensure_path(os.path.join(self.result_dir, key))
        output_filename = os.path.join(
            self.result_dir, key, 'actual vs calculated percentage sum vs delta.svg')

        def plot():
            df = self.simulation_results
            actual = []
            calculated = []
            deltas = []
            for a, c in zip(df['actual percentage'], df['calculated percentage']):
                if hasattr(c, '__iter__'):
                    actual.extend([np.sum(a)] * len(c))
                    calculated.extend([np.sum(c)] * len(c))
                    deltas.extend((np.array(a) - np.array(c)))
                else:
                    actual.append(a)
                    calculated.append(c)
                    deltas.append((a - c))
            plot_df = pd.DataFrame({'calculated percentage sum': calculated, 'delta': deltas, 'actual percentage sum': actual})
            plot_df = plot_df[~np.isnan(plot_df['calculated percentage sum'])]
            p = gg.Plot(plot_df)
            p.theme_bw()
            p.add_scatter('actual percentage sum',
                          'delta', alpha=0.3, size=1)
            p.add_ab_line(0, 0, alpha=0.5)
            p.scale_x_continuous(limits=[0, 1])
            p.scale_y_continuous(limits=[-1, 1])
            p.title(self.name + ' ' + key)
            p.render(output_filename, width=11.69, height=8.27)
        return ppg.FileGeneratingJob(output_filename, plot).depends_on(self._load_simulation_results(simulate_jobs, key)).depends_on(self.load())

    def _plot_mae_improvement(self, simulate_jobs, key):
        ensure_path(os.path.join(self.result_dir, key))
        output_filename = os.path.join(
            self.result_dir, key, 'mae improvement.svg')

        def plot():
            df = self.simulation_results
            df.insert(0, 'deltaMAE', df[
                      'absolute_error_corrected_mean'] - df['absolute_error_observed_mean'])
            if np.isnan(df['deltaMAE']).all():
                p = gg.Plot(pd.DataFrame({"x": [0], 'y': [0], 'text': ['no_data']}))
                p.add_text('x', 'y', 'text')
                p.render(output_filename, width=11.69, height=8.27)
                return
            p = gg.Plot(df)
            p.add_box_plot('sample_tissue', 'deltaMAE')
            p.add_horizontal_bar(0)
            p.facet('contamination_tissue')
            p.render(output_filename, width=11.69, height=8.27)
        return ppg.FileGeneratingJob(output_filename, plot).depends_on(self._load_simulation_results(simulate_jobs, key)).depends_on(self.load())

    def _plot_mae_improvement_calc_over_30(self, simulate_jobs, key):
        ensure_path(os.path.join(self.result_dir, key))
        output_filename = os.path.join(
            self.result_dir, key, 'mae improvement, calc percentage below 0.3.svg')

        def plot():
            df = self.simulation_results
            df = df[df['calculated percentage'] < 0.3]
            df.insert(0, 'deltaMAE', df[
                      'absolute_error_corrected_mean'] - df['absolute_error_observed_mean'])
            if np.isnan(df['deltaMAE']).all():
                p = gg.Plot(pd.DataFrame({"x": [0], 'y': [0], 'text': ['no_data']}))
                p.add_text('x', 'y', 'text')
                p.render(output_filename, width=11.69, height=8.27)
                return

            p = gg.Plot(df)
            p.add_box_plot('sample_tissue', 'deltaMAE')
            p.add_horizontal_bar(0)
            p.facet('contamination_tissue')
            p.render(output_filename, width=11.69, height=8.27)
        return ppg.FileGeneratingJob(output_filename, plot).depends_on(self._load_simulation_results(simulate_jobs, key)).depends_on(self.load())
