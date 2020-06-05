# run_cellranger.pl

Create and run a snakemake pipeline that launches cellranger commands to biowulf.

View the entire readme with `run_cellranger.pl -h`

## Parts of the cellranger workflow

This script produces a Snakemake file that will invoke a series of Cellranger
commands to analyze single-cell RNA-seq data.  The Cellranger steps invoked are:

 1. mkfastq
 2. count
 3. mat2csv

Additionally, Seurat objects will be made for each sample.

## Preparing for the run

Input consists of a tab-delimited file with the following column order, with each
row representing a sample:

 1. Library ID
 2. Sample ID [must be unique]
 3. Raw Data Directory
 4. Lane
 5. Sample Barcode Index
 6. Group 

'Sample ID' must be unique.
'Group' determines what samples will be merged together during the Seurat merge 
step.  All samples with a given Group ID will be merged for that group.  Samples 
can be in more than one group (separate Group ID values with commas), or can be
ignored during merging by setting this field to 'NA'.

## TODO:

	1. [ ] Default to using the biomart file found in ./bin, recreate in workdir if asked.
	2. [ ] Better documentation describing setup, etc.
	3. [ ] Create tSNE, UMAP plots as pdfs
	4. [ ] Add FindAllMarkers and create various heatmaps
	5. [ ] Add EnrichR after FindAllMarkers
	6. [ ] Add filtration of various 'noisy' genes
	7. [ ] Cleanup intermediate files (could be done using native snakemake nomenclature)
	8. [ ] Restore %dopar% functionality for main processing loop of seurat objects.

## Author

  Jason Inman
  inmanjm@nih.gov

