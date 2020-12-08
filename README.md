# check-strand
Command line tool to infer library strandedness for RNA-Seq samples from FASTQs. 
Logic to parse `kallisto` files from Hong Zheng's [blogpost](https://fishycat.netlify.app/en/2017/08/strandness_in_rnaseq/).

## Requirements
* `kallisto`
* `python3`
* `pandas`

## Usage
```
check.py [kallisto index] [samples]
```
The sample file must be a CSV file with at least `fq1` and `fq2` as columns.

The program saves a file, `strands.csv` to the current directory. This file is a copy of the sample
file with an extra column, `strand`, denoting the strand type for each sample.
