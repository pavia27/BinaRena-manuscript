Travelers' Diarrhea (TD) metagenomes
------

Data source:

- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA382010

Citation:

- Zhu Q, Dupont CL, Jones MB, Pham KM, Jiang ZD, DuPont HL, Highlander SK.
  Visualization-assisted binning of metagenome assemblies reveals potential
  new pathogenic profiles in idiopathic travelers’ diarrhea. Microbiome.
  2018 Dec;6(1):1-20.

Data files:

- 76.tsv.xz: Sample #76, contig properties
- 76.Enterobacteriaceae.map: Sample #76, mapping of contigs to Enterobacteria-
  ceae specific marker genes
- 50076.tsv.xz: Sample #50076, contig properties
- 50076.Escherichia_coli.map Sample #50076, mapping of contigs to Escherichia
  coli specific marker genes
- Enterobacteriaceae.lst: CheckM marker gene list for family Enterobacteriaceae
- Escherichia_coli.lst: CheckM marker gene list for species Escherichia coli

Instruction (#76 for example):

1. Launch BinaRena.
2. Unzip 76.tsv.xz, drag & drop into BinaRena.
3. Select "tSNE1_k6" and "tSNE2_k6" for x- and y-axes, respectively.
4. Drag & drop 76.Enterobacteriaceae.map into BinaRena, select "feature set"
   data type.
5. Drag & drop Enterobacteriaceae.lst into BinaRena, associate it with "76.
   Enterobacteriaceae".
6. Select a group of contigs. In the "selected" panel, click dropdown menu,
   then click "Feature group" - "Enterobacteriaceae".
