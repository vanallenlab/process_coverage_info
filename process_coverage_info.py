import matplotlib
matplotlib.use('Agg') # If plots don't display, try commenting out this line

from multiprocessing.dummy import Pool as ThreadPool

import argparse
import subprocess
import sys
import glob
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


"""
Generate graphs and compilations of the data generated by a run of GATK's DepthOfCoverage tool.
Example usage:
python process_coverage_info.py
--label
test
--sample
example_sample_summary_files
--cutoff
25
--output
outputs
--interval
example_sample_interval_summary_files
--gene
example_sample_gene_summary_files

"""


def generate_interval_or_gene_coverage_data(interval_or_gene_folder_path, cutoff, cohort_label, output_folder,
                                            gene_or_interval):
    assert(gene_or_interval in ['gene', 'interval'])
    index_col_name = 'Target' if gene_or_interval == 'interval' else 'Gene'

    all_combined = pd.DataFrame()

    sys.stdout.write('Looking at {} samples\n'.format(len(glob.glob('{}/*'.format(interval_or_gene_folder_path)))))
    for filename in glob.glob('{}/*'.format(interval_or_gene_folder_path)):
        try:
            interval_summary = pd.read_csv(filename, sep='\t', header='infer')
            columns = interval_summary.columns

            for column in columns:
                if '%_above_15' in column:
                    fifteen_column_name = column

            slimmed = interval_summary[[index_col_name, fifteen_column_name]]
            sample_name = fifteen_column_name.split('_%')[0]

            slimmed[sample_name] = slimmed[fifteen_column_name]
            slimmed = slimmed.drop(fifteen_column_name, axis=1)

            if all_combined.empty:
                all_combined = slimmed
            else:
                all_combined = pd.concat([all_combined, slimmed[[sample_name]]], axis=1)
        except:
            sys.stdout.write("Ran into issues on file {}\n".format(filename))

    all_combined.index = all_combined[index_col_name]
    means = all_combined.drop(index_col_name, axis=1).mean(axis=1)
    means.to_csv('{}/{}_average_{}_fraction_above_15_coverage.tsv'.format(output_folder, cohort_label, gene_or_interval),
                 sep='\t')

    filename = '{}/{}_{}_fraction_above_15_coverage_per_sample.tsv'.format(output_folder, cohort_label, gene_or_interval)
    all_combined.to_csv(filename,
                        sep='\t', index=False)
    subprocess.call('gzip {}'.format(filename).split())

    plt.figure()
    plt.hist(np.array(means[~np.isnan(means)]))
    plt.title('Mean Fraction Covered at >= 15, Distribution at {} Level:\n{} Cohort'.format(
        gene_or_interval.capitalize(),
        cohort_label))
    plt.ylabel('{}s'.format(gene_or_interval.capitalize()))
    plt.xlabel('{} Mean Coverage'.format(gene_or_interval.capitalize()))
    plt.savefig('{}/{}_{}_mean_coverage'.format(output_folder, cohort_label.replace(' ', '_'), gene_or_interval))


def get_sample_id_and_mean_coverage(pool_args):
    filename = pool_args.get('filename')
    try:
        info = pd.read_csv(filename, sep='\t', header='infer')
        sample_id = info['sample_id'].loc[0]
        sample_mean_coverage = info['mean'].loc[1]
    except:
        sys.stdout.write("Couldn't read sample summary info from {}\n".format(f))
        sample_id = None
        sample_mean_coverage = None
    return {'sample_id': sample_id,
            'sample_mean_coverage': sample_mean_coverage}


def generate_sample_mean_coverage_data(sample_folder_path, cutoff, cohort_label, output_folder):
    sys.stdout.write('Generating sample mean coverage graph\n')
    mean_coverages = []
    sample_id_to_mean_coverage = {}

    sys.stdout.write('Looking at {} samples\n'.format(len(glob.glob('{}/*'.format(sample_folder_path)))))
    pool_args = []
    pool = ThreadPool(10)

    for f in glob.glob('{}/*'.format(sample_folder_path)):
        pool_args.append({'filename': f})

    for result in pool.imap(get_sample_id_and_mean_coverage, pool_args):
        sample_mean_coverage = result.get('sample_mean_coverage')
        sample_id = result.get('sample_id')
        sys.stdout.write('Got info for sample {}\n'.format(sample_id))
        if sample_id:
            mean_coverages.append(sample_mean_coverage)
            sample_id_to_mean_coverage[sample_id] = sample_mean_coverage
            sys.stdout.write('{}\t{}\n'.format(sample_id, sample_mean_coverage))

    sys.stdout.write("Writing sample level mean coverage output\n")
    # Output list of all samples with corresponding mean sample coverage and whether they passed the coverage cutoff
    with open('{}/{}_mean_coverage_by_sample_id.tsv'.format(output_folder, cohort_label), 'w') as f:
        f.write('sample_id\tmean_sample_coverage\tpassed_threshold\n')
        for sample_id, mean_coverage in sample_id_to_mean_coverage.items():
            f.write('{}\t{}\t{}\n'.format(sample_id, mean_coverage, mean_coverage >= cutoff))


def main():
    parser = argparse.ArgumentParser(description='Generate cohort coverage statistics and graphs')

    parser.add_argument('--label', metavar='label', type=str)
    # Folder where sample summary files are
    parser.add_argument('--sample', metavar='bam_file', type=str)
    parser.add_argument('--cutoff', metavar='cutoff', type=int)
    # Folder where sample interval summary files are
    parser.add_argument('--interval', metavar='chromosome', type=str)
    # Folder where sample gene summary files are
    parser.add_argument('--gene', metavar='position', type=str)

    parser.add_argument('--output', metavar='metavar', type=str)
    args = parser.parse_args()

    label = args.label
    cutoff = args.cutoff
    sample_folder = args.sample
    interval_folder = args.interval
    gene_folder = args.gene
    output = args.output

    generate_sample_mean_coverage_data(sample_folder, cutoff, label, output)
    generate_interval_or_gene_coverage_data(interval_folder, cutoff, label, output, 'interval')
    generate_interval_or_gene_coverage_data(gene_folder, cutoff, label, output, 'gene')


if __name__ == '__main__':
    main()