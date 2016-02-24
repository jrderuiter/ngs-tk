import argparse

import numpy as np
import pandas as pd
import pysam

from tqdm import tqdm


def region_coverage(bam_file, region, agg_func=np.mean):
    contig, start, end = region

    # Get coverage within region.
    pileups = bam_file.pileup(contig, start, end, stepper='all', truncate=True)
    coverage = [p.n for p in pileups]

    # Handle empty case.
    if len(coverage) == 0:
        coverage = [0]

    return agg_func(coverage)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('bam')
    parser.add_argument('intervals')
    parser.add_argument('output')

    parser.add_argument('--median', default=False, action='store_true')
    parser.add_argument('--verbose', default=False, action='store_true')

    return parser.parse_args()


def main(args):
    # Get bam file.
    bam_file = pysam.AlignmentFile(args.bam)

    # Get intervals from bed file.
    bed = pd.read_csv(args.intervals, sep='\t', header=None)
    bed.columns = ['chrom', 'chromStart', 'chromEnd']

    # Calculate coverage for each region.
    regions = zip(bed['chrom'], bed['chromStart'], bed['chromEnd'])

    if args.verbose:
        regions = tqdm(regions, total=len(bed))

    agg_func = np.median if args.median else np.mean
    coverage = [region_coverage(bam_file, r, agg_func=agg_func)
                for r in regions]

    # Write output.
    bed['coverage'] = coverage
    bed.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    args = parse_args()
    main(args)
