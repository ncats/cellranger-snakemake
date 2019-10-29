# run_cellranger.pl

Create and run a snakemake pipeline that launches cellranger commands to biowulf.

View the entire readme with `run_cellranger.pl -h`

## Parts of the cellranger workflow

This script produces a Snakemake file that will invoke a series of Cellranger
commands to analyze single-cell RNA-seq data.  The steps invoked are:

 1. mkfastq
 2. count
 3. aggr [Optional]
 4. mat2csv

Additionally, Seurat objects will be made for each sample.  [Not yet implemented].

## Preparing for the run

Input consists of a tab-delimited file with the following column order, with each
row representing a sample:

 1. Library ID
 2. Sample ID [must be unique]
 3. Raw Data Directory
 4. Lane
 5. Sample Barcode Index
 6. 'Batch' [Optional; Not yet implemented]

## TODO:

[ ] Add a step to create Seurat objects per sample
[ ] Add a step to merge Seurat objects per batch
[ ] Cleanup intermediate files (could be done using native snakemake nomenclature)

## Author

  Jason Inman
  inmanjm@nih.gov

