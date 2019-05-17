import numpy as np
import argparse
import pandas as pd
from multiprocessing.dummy import Pool as ThreadPool
from collections import defaultdict
import sys
from scipy.stats import mannwhitneyu


def mann_whitney_u(pool_args):
    """Conduct the Mann Whitney statistical test, returning the significance and directionality of the comparison"""
    gene_or_interval = pool_args.get('gene_or_interval')
    case_results = np.array(pool_args.get('case_results'))
    control_results = np.array(pool_args.get('control_results'))
    case_median = np.median(case_results)
    control_median = np.median(control_results)

    try:
        statistic, pvalue_one_sided = mannwhitneyu(case_results, control_results)

        directionality = 'same'
        if pvalue_one_sided < 0.05:
            directionality = 'cases_higher' if case_median > control_median else 'controls_higher'
    except ValueError:
        if case_median == case_results.min() == case_results.max() == control_median == control_results.min() == control_results.max():
            directionality = 'values_all_identical'
        else:
            directionality = "ValueError"
        pvalue_one_sided = None

    return {
        'index': gene_or_interval,
        'directionality': directionality,
        'pvalue': pvalue_one_sided,
        'case_median': case_median,
        'control_median': control_median
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
                 'control_results': control_results})
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

    sys.stdout.write("Running Mann-Whitney U test\n")
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