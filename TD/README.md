Travelers' Diarrhea (TD) metagenomes
------

## Data source

- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA382010

## Citation

- Zhu Q, Dupont CL, Jones MB, Pham KM, Jiang ZD, DuPont HL, Highlander SK.
  Visualization-assisted binning of metagenome assemblies reveals potential
  new pathogenic profiles in idiopathic travelers’ diarrhea. Microbiome.
  2018 Dec;6(1):1-20.

## Data files

### Sample #76

- 76.fna.xz: contig sequences
- 76.tsv.xz: contig properties
- 76.bins.tsv: contig lists of investigated bins
- 76.Enterobacteriaceae.map: mapping of contigs to _Enterobacteriaceae_
  specific marker genes
- Enterobacteriaceae.lst: CheckM marker gene list for family
  _Enterobacteriaceae_

### Sample #50076

- 50076.fna.xz: contig sequences
- 50076.tsv.xz: contig properties
- 50076.bins.tsv: contig lists of investigated bins
- 50076.Escherichia_coli.map mapping of contigs to _Escherichia coli_ specific
  marker genes
- Escherichia_coli.lst: CheckM marker gene list for species _Escherichia coli_

## Script files

- TD_tsv_build.R: R script for constructing BinaRena input files

## Instruction

(#76 for example)

1. Unzip 76.tsv.xz, drag & drop into BinaRena.
2. Select "tSNE1_k6" and "tSNE2_k6" for x- and y-axes, respectively.
3. Drag & drop 76.Enterobacteriaceae.map into BinaRena, select "feature set"
   data type.
4. Drag & drop Enterobacteriaceae.lst into BinaRena, associate it with "76.
   Enterobacteriaceae".
5. Select a group of contigs. In the "selected" panel, click dropdown menu,
   then click "Feature group" - "Enterobacteriaceae".
