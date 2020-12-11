#!/usr/bin/env python

import os
import sys
import tempfile
import gzip
import argparse
import subprocess as sb
import pandas as pd

DEFAULT_NREADS = 5e4
DEFAULT_THREADS = 2

def get_args():
    parser = argparse.ArgumentParser(description='Check RNA-Seq stranding')
    parser.add_argument('--index', '-i', required=True, type=str, nargs=1, help='kallisto index file', 
                        dest='index')
    parser.add_argument('--num-reads', '-n', type=int, nargs=1, help='number of reads to check (default: 5e4)',
                        dest='num-reads')
    parser.add_argument('--threads', '-t', type=int, nargs=1, help='number of threads for kallisto (default: 2)',
                        dest='threads')
    parser.add_argument('samples', type=str, nargs=1, help='CSV file with fqs to check')
    args = parser.parse_args()
    return args

def get_reads(inputfp, outfp, num_reads):
    if inputfp.endswith('gz'):
        with gzip.open(inputfp, 'rb') as f:
            lines = [next(f).decode("utf-8") for i in range(int(num_reads*4))]
    else:
        with open(inputfp, 'r') as f:
            lines = [next(f) for i in range(int(num_reads*4))]
    with open(outfp, 'w') as f:
        for line in lines:
            f.write(line)

def strand_scoring(uncts, rfcts, frcts):
    if rfcts / frcts < .3 and uncts / rfcts > 3:
        return 'stranded'
    elif rfcts / frcts > 3 and uncts / frcts > 3:
        return 'reverse'
    else:
        return 'unstranded'

def infer_strand(fq1, fq2, idx_fp, nreads=5e4, threads=2):
    print(f"Inferring strand from {fq1}, {fq2}")
    with tempfile.TemporaryDirectory() as tmpdirname:
        fq1_fp = os.path.join(tmpdirname, 'test_1.fq')
        fq2_fp = os.path.join(tmpdirname, 'test_2.fq')

        get_reads(fq1, fq1_fp, nreads)
        get_reads(fq2, fq2_fp, nreads)

        sb.run(['kallisto', 'quant', '-i', idx_fp, '-t', str(threads),'-o', os.path.join(tmpdirname, 'un'), 
                fq1_fp, fq2_fp], stdout=sb.DEVNULL, stderr=sb.DEVNULL)
        sb.run(['kallisto', 'quant', '-i', idx_fp, '-t', str(threads),'-o', os.path.join(tmpdirname, 'rf'), 
                '--rf-stranded', fq1_fp, fq2_fp], stdout=sb.DEVNULL, stderr=sb.DEVNULL)
        sb.run(['kallisto', 'quant', '-i', idx_fp, '-t', str(threads),'-o', os.path.join(tmpdirname, 'fr'), 
                '--fr-stranded', fq1_fp, fq2_fp], stdout=sb.DEVNULL, stderr=sb.DEVNULL)

        res_fps = [os.path.join(tmpdirname, d, 'abundance.tsv') for d in ['un', 'rf', 'fr']]
        dfs = [pd.read_csv(fp, sep='\t') for fp in res_fps]
        ct_sums = [df['est_counts'].sum() for df in dfs]
        uncts, rfcts, frcts = ct_sums
        strand = strand_scoring(uncts, rfcts, frcts)
        return strand

def main():
    args = get_args()
    idx_fp = args['index']
    samples_fp = args['samples']
    
    samples = pd.read_csv(samples_fp)
    
    num_reads = args['num-reads']
    num_threads = args['threads']
    
    if num_reads is None:
        num_reads = DEFAULT_NREADS
    if num_threads is None:
        num_threads = DEFAULT_THREADS

    res = []
    for index, row in samples.iterrows():
        res.append(infer_strand(row['fq1'], row['fq2'], idx_fp, num_reads, num_threads))
    samples['strand'] = res
    fn = samples_fp[0:samples_fp.find('.csv')]
    outfn = f"{fn}.stranded.csv"
    samples.to_csv(outfn, index=False, header=True)

if __name__ == '__main__':
    main()
