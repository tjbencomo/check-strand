import os
import sys
import tempfile
import subprocess as sb
import gzip
import pandas as pd

def get_reads(inputfp, outfp, num_reads=5e4):
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

def infer_strand(fq1, fq2, idx_fp):
    with tempfile.TemporaryDirectory() as tmpdirname:
        fq1_fp = os.path.join(tmpdirname, 'test_1.fq')
        fq2_fp = os.path.join(tmpdirname, 'test_2.fq')

        get_reads(fq1, fq1_fp)
        get_reads(fq2, fq2_fp)

        sb.run(['kallisto', 'quant', '-i', idx_fp, '-t', '2','-o', os.path.join(tmpdirname, 'un'), 
                fq1_fp, fq2_fp], stdout=subprocess.DEVNULL)
        sb.run(['kallisto', 'quant', '-i', idx_fp, '-t', '2','-o', os.path.join(tmpdirname, 'rf'), 
                '--rf-stranded', fq1_fp, fq2_fp], stdout=subprocess.DEVNULL)
        sb.run(['kallisto', 'quant', '-i', idx_fp, '-t', '2','-o', os.path.join(tmpdirname, 'fr'), 
                '--fr-stranded', fq1_fp, fq2_fp], stdout=subprocess.DEVNULL)

        res_fps = [os.path.join(tmpdirname, d, 'abundance.tsv') for d in ['un', 'rf', 'fr']]
        dfs = [pd.read_csv(fp, sep='\t') for fp in res_fps]
        ct_sums = [df['est_counts'].sum() for df in dfs]
        uncts, rfcts, frcts = ct_sums
        strand = strand_scoring(uncts, rfcts, frcts)
        return strand

def main():
    idx_fp = sys.argv[1]
    samples_fp = sys.argv[2]
    samples = pd.read_csv(samples_fp)
    res = []
    for index, row in samples.iterrows():
        res.append(infer_strand(row['fq1'], row['fq2'], idx_fp))
    samples['strand'] = res
    samples.to_csv('strands.csv', index=False, header=True)

if __name__ == '__main__':
    main()
