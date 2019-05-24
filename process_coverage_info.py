import matplotlib
matplotlib.use('Agg') # If plots don't display, try commenting out this line

from multiprocessing.dummy import Pool as ThreadPool

import argparse
import subprocess
import sys
import glob
import pandas as pd
import os

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


def get_interval_or_gene_coverage_for_file(pool_args):
    filename = pool_args.get('filename')
    gene_or_interval = pool_args.get('gene_or_interval')
    index_col_name = 'Target' if gene_or_interval == 'interval' else 'Gene'
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
        return {'slimmed': slimmed, 'sample_name': sample_name}
    except:
        sys.stdout.write("Ran into issues on file {}\n".format(filename))
        return {'slimmed': None, 'sample_name': None}


def generate_interval_or_gene_coverage_data(interval_or_gene_folder_path, cutoff, cohort_label, output_folder,
                                            gene_or_interval):
    assert(gene_or_interval in ['gene', 'interval'])
    index_col_name = 'Target' if gene_or_interval == 'interval' else 'Gene'

    sys.stdout.write('Gene/Interval: Looking at {} samples\n'.format(len(glob.glob('{}/*'.format(interval_or_gene_folder_path)))))

    pool_args = []
    pool = ThreadPool(16)

    total_samples = 0
    for filename in glob.glob('{}/*'.format(interval_or_gene_folder_path)):
        total_samples += 1
        sys.stdout.write('Gene/Interval: {}\n'.format(filename))
        pool_args.append({'filename': filename,
                          'gene_or_interval': gene_or_interval})

    sys.stdout.write('Total samples: {}\n'.format(total_samples))
    sys.stdout.write('Gathering gene/interval results\n')

    results_gathered = 0
    for i,  result in enumerate(pool.imap(get_interval_or_gene_coverage_for_file, pool_args)):
        results_gathered += 1
        sys.stdout.write('Processed {}/{} samples for {} info\n'.format(results_gathered,
                                                                        total_samples,
                                                                        gene_or_interval))
        sample_name = result.get('sample_name')
        slimmed = result.get('slimmed')
        if sample_name:
            output = slimmed[[sample_name]]
            if i == 0:
                # Output the index for the first one, but not any subsequent ones so we can efficiently paste
                # them together after
                output.index = slimmed[index_col_name]
                output.to_csv('{}_{}_{}_tmp.tsv'.format(i, sample_name, gene_or_interval), sep='\t', index=True)
            elif i > 0:
                output.to_csv('{}_{}_{}_tmp.tsv'.format(i, sample_name, gene_or_interval), sep='\t', index=False)

    sys.stdout.write('Outputting the results from each sample for {}...\n'.format(gene_or_interval))

    final_filename = '{}/{}_{}_fraction_above_15_coverage_per_sample.tsv'.format(output_folder, cohort_label,
                                                                           gene_or_interval)
    sys.stdout.write('Concatenating using paste for {}\n'.format(gene_or_interval))

    command = 'paste *{}_tmp.tsv > {}; rm *{}_tmp.tsv'.format(gene_or_interval, final_filename, gene_or_interval)

    # This command can be used if the ulimit (number of simulataneously open file descriptors) cannot be increased
    # command = 'for f in *{}_tmp.tsv; do cat {} | paste - $f > temp; cp temp {}; done; rm temp; rm *{}_tmp.tsv'\
    #    .format(gene_or_interval, final_filename, final_filename, gene_or_interval)

    sys.stdout.write('Using command:\n{}\n'.format(command))
    os.system(command)

    sys.stdout.write('Gzipping {}\n'.format(final_filename))
    subprocess.call('gzip {}'.format(final_filename).split())


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

    sys.stdout.write('Loading each sample info...\n')
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

    if sample_folder:
        generate_sample_mean_coverage_data(sample_folder, cutoff, label, output)
    if interval_folder:
        generate_interval_or_gene_coverage_data(interval_folder, cutoff, label, output, 'interval')
    if gene_folder:
        generate_interval_or_gene_coverage_data(gene_folder, cutoff, label, output, 'gene')


if __name__ == '__main__':
    main()