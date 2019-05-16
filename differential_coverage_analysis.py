import matplotlib

import argparse
import subprocess
import sys
import glob
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import argparse
import os
import pandas as pd
from multiprocessing.dummy import Pool as ThreadPool
from collections import defaultdict
import sys


def wilcoxon_test_genes_or_intervals_coverage(gene_or_interval,
                                              fractions_case_file,
                                              fractions_control_file,
                                              sample_to_exclude_list,
                                              threads=8):

    fractions_case = pd.read_csv(fractions_case_file, sep='\t', index_col=0, header='infer')
    fractions_controls = pd.read_csv(fractions_control_file, sep='\t', index_col=0, header='infer')

    pool_args = []
    pool = ThreadPool(threads)


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
    parser.add_argument('--samples_to_exclude', metavar='samples_to_exclude', type=int)

    # Gene or interval argument should be 'gene' or 'interval' (duh)
    parser.add_argument('--gene_or_interval', metavar='gene_or_interval', type=str)
    parser.add_argument('--output', metavar='output', type=str)
    args = parser.parse_args()

    fractions_case = args.fractions_case
    fractions_control = args.fractions_control
    samples_to_exclude = args.samples_to_exclude
    gene_or_interval = args.gene_or_interval
    output = args.output

    samples_to_exclude = [] #load_samples_to_exclude(samples_to_exclude)

    wilcoxon_test_genes_or_intervals_coverage(gene_or_interval,
                                              fractions_case,
                                              fractions_control,
                                              samples_to_exclude,
                                              threads=8)


if __name__ == '__main__':
    main()