import numpy as np
import argparse
import pandas as pd
from multiprocessing.dummy import Pool as ThreadPool
from collections import defaultdict
import sys
from scipy.stats import mannwhitneyu, fisher_exact


def mann_whitney_u(pool_args):
    """Conduct the Mann Whitney statistical test, returning the significance and directionality of the comparison"""
    gene_or_interval = pool_args.get('gene_or_interval')
    case_results = np.array(pool_args.get('case_results'))
    control_results = np.array(pool_args.get('control_results'))

    try:
        case_median = np.median(case_results)
    except ValueError:
        case_median = 0
    try:
        control_median = np.median(control_results)
    except ValueError:
        control_median = 0
    try:
        case_mean = np.mean(case_results)
    except ValueError:
        case_mean = 0
    try:
        control_mean = np.mean(control_results)
    except ValueError:
        control_mean = 0

    # Fisher's exact on the means
    try:
        mean_fisher_OR, mean_fisher_p = fisher_exact([[control_mean * len(control_results), len(control_results)],
                                             [case_mean * len(case_results), len(case_results)]])
        if mean_fisher_p > 0.05:
            mean_fisher_directionality = 'same'
        else:
            if mean_fisher_OR > 1:
                mean_fisher_directionality = 'controls_higher'
            else:
                mean_fisher_directionality = 'cases_higher'
    except ValueError:
        mean_fisher_OR, mean_fisher_p, mean_fisher_directionality = 0 , 0, 0

    # Fisher's exact on the medians instead of the means? This ensures that outliers don't skew the results
    try:
        fisher_OR, fisher_p = fisher_exact([[control_median * len(control_results), len(control_results)],
                                             [case_median * len(case_results), len(case_results)]])
        if fisher_p > 0.05:
            fisher_directionality = 'same'
        else:
            if fisher_OR > 1:
                fisher_directionality = 'controls_higher'
            else:
                fisher_directionality = 'cases_higher'
    except ValueError:
        fisher_p, fisher_OR, fisher_directionality = 0, 0, 0

    return {
        'index': gene_or_interval,
        'case_mean': case_mean,
        'control_mean': control_mean,
        'case_median': case_median,
        'control_median': control_median,
        'medians_fisher_OR': fisher_OR,
        'medians_fisher_pvalue': fisher_p,
        'medians_fisher_directionality': fisher_directionality,
        'means_fisher_OR': mean_fisher_OR,
        'means_fisher_pvalue': mean_fisher_p,
        'means_fisher_directionality': mean_fisher_directionality,
        'quantile_20_cases': np.percentile(case_results, 20),
        'quantile_20_controls': np.percentile(control_results, 20),
        'quantile_40_cases': np.percentile(case_results, 40),
        'quantile_40_controls': np.percentile(control_results, 40),
        'quantile_60_cases': np.percentile(case_results, 60),
        'quantile_60_controls': np.percentile(control_results, 60),
        'quantile_80_cases': np.percentile(case_results, 80),
        'quantile_80_controls': np.percentile(control_results, 80),
        'quantile_100_cases': np.percentile(case_results, 100),
        'quantile_100_controls': np.percentile(control_results, 100)
    }


def build_pool_arguments(fractions_case_df, fractions_control_df):
    """Build up the array with arguments that will be fed to the thread pool"""
    pool_args = []

    # For each interval/gene, do a differential coverage analysis between the samples in the cohort vs in the control
    total_intervals = len(fractions_case_df.index)
    for index, interval in enumerate(fractions_case_df.index):
        sys.stderr.write('Adding interval to thread pool: {} ({}/{})\n'.format(interval, index, total_intervals))
        case_results = fractions_case_df.loc[interval]
        control_results = fractions_control_df.loc[interval]

        args = ({'gene_or_interval': interval,
                 'case_results': case_results,
                 'control_results': control_results
                 })
        pool_args.append(args)
    return pool_args


def get_output_df_from_pool_results(results):
    """Process the results derived from running the multithreaded statistical tests"""
    output_column_lists = defaultdict(lambda: [])
    for r in results:
        for k, v in r.items():
            output_column_lists[k].append(v)

    final_df = pd.DataFrame(output_column_lists)
    final_df.index = final_df['index']
    final_df = final_df.drop('index', axis=1)
    return final_df


def exclude_samples(fractions_df, samples_to_exclude_list):
    """Remove all the columns representing the samples in the given list of samples to exclude"""
    sys.stdout.write("Excluding samples: {}\n".format(samples_to_exclude_list))
    all_columns = fractions_df.columns
    columns_to_keep = [c for c in all_columns if c not in samples_to_exclude_list]
    return fractions_df[columns_to_keep]


def test_genes_or_intervals_coverage(
        fractions_case_file,
        fractions_control_file,
        sample_to_exclude_list,
        threads=8):
    """Set up a multithreaded approach to running the Mann Whitney on each interval/gene between control and case"""
    fractions_case = pd.read_csv(fractions_case_file, sep='\t', index_col=0, header='infer')
    fractions_control = pd.read_csv(fractions_control_file, sep='\t', index_col=0, header='infer')

    # Check that the indices of the cases and controls are exactly the same
    case_indices = list(fractions_case.index)
    control_indices = list(fractions_control.index)
    for i, case_i in enumerate(case_indices):
        if not case_i == control_indices[i]:
            sys.exit('Mismatched indices: {} is not equal to {}\n'.format(case_i, control_indices[i]))

    sys.stdout.write('Number of samples in cases: {}\n'.format(len(fractions_case.columns)))
    fractions_case = exclude_samples(fractions_case, sample_to_exclude_list)
    sys.stdout.write('Number of samples in cases after removing samples: {}\n'.format(len(fractions_case.columns)))

    sys.stdout.write('Number of samples in controls: {}\n'.format(len(fractions_control.columns)))
    fractions_control = exclude_samples(fractions_control, sample_to_exclude_list)
    sys.stdout.write('Number of samples in controls after removing samples: {}\n'.format(len(fractions_control.columns)))

    pool = ThreadPool(threads)
    pool_args = build_pool_arguments(fractions_case, fractions_control)

    sys.stdout.write("Running tests\n")
    num_finished = 0
    results = []
    # Execute the mann whitney comparison on each of the pool arguments, appending the outputs
    # to the results array as the results return.
    for result in pool.imap(mann_whitney_u, pool_args):
        results.append(result)
        num_finished += 1
        if num_finished % 250 == 0:
            sys.stderr.write('{}/{} finished\n'.format(num_finished, len(pool_args)))
        if num_finished == len(pool_args):
            sys.stdout.write('{}/{} finished. Done\n'.format(num_finished, len(pool_args)))

    pool.close()
    pool.join()

    final_df = get_output_df_from_pool_results(results)
    return final_df


def load_samples_to_exclude(samples_to_exclude_file):
    """Load the list of samples that should be excluded from this differential coverage analysis (usually because of
    relatedness issues or low sample mean coverage level).

    Format: One column where each row is a different sample ID.
    """
    samples = pd.read_csv(samples_to_exclude_file)
    return samples


def main():
    parser = argparse.ArgumentParser(description='Compare coverage between cohorts at each gene and interval')

    parser.add_argument('--fractions_case', metavar='fractions_case', type=str)
    parser.add_argument('--fractions_control', metavar='fractions_control', type=str)
    parser.add_argument('--samples_to_exclude', metavar='samples_to_exclude', type=str)

    # Gene or interval argument should be 'gene' or 'interval' (duh)
    parser.add_argument('--gene_or_interval', metavar='gene_or_interval', type=str, default="")
    parser.add_argument('--output_folder', metavar='output_folder', type=str)
    args = parser.parse_args()

    # Process arguments
    fractions_case = args.fractions_case
    fractions_control = args.fractions_control

    samples_to_exclude = args.samples_to_exclude
    if samples_to_exclude is None:
        samples_to_exclude = []
    else:
        samples_to_exclude = list(pd.read_csv(samples_to_exclude, header=None, index_col=0).index)

    gene_or_interval = args.gene_or_interval
    output_folder = args.output_folder

    differential_coverage_report = test_genes_or_intervals_coverage(
        fractions_case,
        fractions_control,
        samples_to_exclude,
        threads=8
    )

    differential_coverage_report.to_csv('{}/{}_dc_report.tsv'.format(output_folder,
                                                                     gene_or_interval), sep='\t')


if __name__ == '__main__':
    main()