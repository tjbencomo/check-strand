# check-strand
Command line tool to infer library strandedness for RNA-Seq samples from FASTQs. 
Logic to parse `kallisto` files from Hong Zheng's [blogpost](https://fishycat.netlify.app/en/2017/08/strandness_in_rnaseq/).

## Requirements
* `kallisto`
* `python3`
* `pandas`

## Installation
1. Ensure that the above requirements are installed. `kallisto` should be callable from the command line.
2. Clone the repository
3. (Optional) `chmod +x check-strand.py` to run without calling `python`
4. (Optional) Add the directory to your `$PATH` to use anywhere

## Usage
```
check-strand.py --index INDEX [--num-reads NUM-READS]
                       [--threads THREADS]
                       samples
```

The sample file must be a CSV file with at least `fq1` and `fq2` as columns.

The program saves a file, `[samples].stranded.csv` to the current directory. `[samples]` is the input filename.
This file is a copy of the sample file with an extra column, `strand`, denoting the strand type for each sample.
