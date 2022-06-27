# BinaRena-manuscript
These are scripts used in the analysis of the datasets used for validation of [BinaRena](https://github.com/qiyunlab/binarena), from the paper XXX.

### MAQ and TD metagenomes

 - MAQ_tsv_build.R
 - TD_tsv_build.R
 
These two scripts are very similar and essentially merge contig descriptors form various programs (contig-to-bin assingmetns, kraken output, dimensionally reduced kmer frequency, etc.). The td folder also provides a walkthrough for procesing sample #76

### CAMI 2 Marine metagenomes

  - marine_tsv_build.R
 
 This script merges contig-to-bin assingments, coverage, and tetranucleotide frequencies (4kmer). It then runs a principal component analysis on a matirx of 4kmer.  
