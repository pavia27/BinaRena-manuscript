Maquia (MAQ) peatland metagenomes
------

## Data sources

- JGI: Ga0314862-Ga0314867

## Data files

- MAQ_BinaRena_Input.tsv.xz: contig properties (input file for BinaRena)
- MAQ_contigs.fna.xz.00-05: contig sequences
  - To extract: `cat MAQ_contigs.fna.xz.0* | xz -d > MAQ_contigs.fna`

## Script files

- MAQ_tsv_build.R: R script for constructing BinaRena input file.
