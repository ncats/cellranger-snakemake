#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))
library(data.table)
library(biomaRt)
library(stringr)
library(iterators)
library(optparse)
suppressPackageStartupMessages(library(doParallel))


# NOTE: This script assumes:
# 1. That the files within the supplied directory are all part of the same run
# (This assumption means we can convert a single file's Ensemble Gene IDs into
# Gene Symbols via Ensembl Mart and it will map to all other input files)
# 2. That all .csv files within the supplied directory are intended to be used.
# (This assumption means this script should be carefully used when not invoked
# as part of a pipeline where this assumption can be trusted.)

workdir.path <- '.'
map_file_path <- './mapping_file'
num_cores <- 7

# Handle args for command-line runs
if ( !interactive()) {
  
  option_list = list(
      
    make_option(c("-w", "--workdir"), type="character", default=".", help="path to working dir", metavar="character"),
    
    make_option(c("-m", "--mapfile"), type="character", default="mapping_file", help="mapping file name [default= %default]", metavar="character"),
    
    make_option(c("-c", "--numcores"), type="integer", default=7, help="Number of available cores for dopar tasks", metavar="character")
    
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)

  if (!is.null(opt$workdir)) {
    workdir.path <- opt$workdir
  }
  
  if (!is.null(opt$mapfile)) {
    map_file_path <- opt$mapfile
  }

  if (!is.null(opt$numcores)) {
    num_cores <- opt$numcores
  }
  
}

registerDoParallel(cores=num_cores)

setwd(workdir.path)
files <- Sys.glob('*.csv')

countfiles.path <- file.path(workdir.path, 'Countfiles_gene_symbol')
if (!dir.exists(countfiles.path) ) {
  dir.create( countfiles.path )
}
seurat_files.path <- file.path(workdir.path, 'seurat_files')
if (!dir.exists( seurat_files.path) ) {
  dir.create( seurat_files.path )
}


# Map the first file's Ensemble IDs to gene symbols:
ensg.genes.path <- file.path(workdir.path,'ensg.genes')
if ( file_test( '-f', ensg.genes.path ) ) {
  
  # If it already exists, use the pre-existing map.
  message(paste0("Found ",ensg.genes.path,', and will use it instead of making it again.'))
  ensg.genes = fread(ensg.genes.path)
  
} else {
  
  # biomart is notoriously fickle.  Might have to do this a few times
  tries <- 0
  max_tries <- 3
  
  while( !file_test( '-f', ensg.genes.path) && tries < max_tries ) {

    tries <- tries + 1
    message(paste0('Creating ensg2symbol mapping file, attempt ',tries))    
    dt <- as.data.frame(fread(files[1]))
    names(dt)[1] <- 'ENSG'
    dt <- as.data.table(dt)
    ensg.genes <- data.table('ENSG' = dt$ENSG)
    genes <- ensg.genes$ENSG
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    G_list <- getBM(filters='ensembl_gene_id', attributes=c('ensembl_gene_id','external_gene_name'),values=genes,mart=mart)
    G_list <- as.data.table(G_list)
    ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG", by.y="ensembl_gene_id")
    ensg.genes <- na.omit(ensg.genes)
    # Now can use G_list to merge the list from each file together by symbol.  So why not save it for this run?
    fwrite(ensg.genes, ensg.genes.path, col.names=T, row.names=F, sep=',')
    
  }
  
  if (!file_test( '-f', ensg.genes.path) ) {
    stop("Couldn't create the ensg2symbol mapping.")
  }
  
}

# TODO: Filter out 'blacklisted' genes

# Create group.map based on map_file.  Bring in library but call it 'project' ala Seurat
map_file <- as.data.frame(read.table(map_file_path))
col_names <- c('library','sample_id','file_path','lanes','library_id','group')
colnames(map_file) <- col_names
group.map <- data.frame(matrix(ncol=2,nrow=0))
for (row in 1:nrow(map_file)) {
  
  # store group-to-sample map, with multiple rows for samples mappin to multiple groups
  if ( grepl( ',', map_file[row,'group'])) {

    for (group in as.numeric(unlist(strsplit(as.character(map_file$group[row]),',')))) {

        group.map <- rbind( group.map, 
                            list( map_file$sample_id[row],
                                  group
                                ),
                            stringsAsFactors=FALSE
                          )
    }
    
  } else {
    
    if ( !is.na(map_file[row,'group'])) {
      
      group.map <- rbind( group.map,
                          list( map_file$sample_id[row],
                                map_file$group[row]
                              ),
                          stringsAsFactors=FALSE
                        )
      
    }
    
  }
  
}
colnames(group.map) <- c('sample_id','group')

# Just kinda pulling this out to avoid wall o' code
reformat_for_seurat <- function(x, samplename){
  x <- x[external_gene_name != "NA",]
  x <- as.data.frame(x)
  
  #message("sample: ",samplename)
  
  row.names(x) <- x$external_gene_name
  x <- x[-1]
  barcodes <- names(x)
  barcodes <- paste0(barcodes, paste0(".", samplename))
  names(x) <- barcodes
  x <- CreateSeuratObject( counts = x, 
                           project = as.character(map_file$library[map_file$sample_id == samplename])
                          )
  x@meta.data$sample <- samplename
  return(x)
}


seurat_files <- vector()
foreach( file=iter(files) ) %dopar% {

  sample_id <- str_split_fixed( file, ".csv", 2)[1]
  message(paste("Working on",sample_id,sep=" "))
  
  dt <- as.data.frame(fread(file))
  names(dt)[1] <- 'ENSG'
  dt <- as.data.table(dt)
  dt <- subset(dt, dt$ENSG %in% ensg.genes$ENSG )
  dt <- merge(dt, ensg.genes, by = 'ENSG')
  dt <- dt[,ENSG:=NULL]
  dt <- dt[!duplicated(external_gene_name),]
  cols <- names(dt)[c(1: (length(names(dt))-1) )]
  dt<- subset(dt, select=c("external_gene_name",cols))
  fwrite(dt, paste0(countfiles.path,'/',sample_id,'_gene_symbol.csv'), col.names=T, row.names=F, quote=F, sep=",")
  
  # Make the seurat object (so)
  so <- reformat_for_seurat(dt, sample_id)
  seurat_filename <- paste0(seurat_files.path,'/',sample_id,'.seurat.Rdata')
  #message("seurat filename: ", seurat_filename)
  save(so, file=seurat_filename )
  
  message(paste("Finished with",sample_id,sep=" "))
  
}

# merge each group as neccessary:
message("Merging seurat objects")
setwd(seurat_files.path)
foreach (group_x =iter(sort(unique(group.map$group)))) %dopar% {
  
  message("merging group ",group_x)
  
  group.subset <- subset(group.map,group==group_x)
  sample_one.file <- paste0(group.subset$sample_id[1],'.seurat.Rdata')
  # loading created the object names 'so' because it was named that way initially 
  load(sample_one.file)
  group.merged <- so
    if (nrow(group.subset) > 1) {
    for (sample_id in group.subset[2:nrow(group.subset),]$sample_id) {
      sample_n.file <- paste0(sample_id,'.seurat.Rdata')
      load(sample_n.file)
      group.merged <- merge(group.merged,so)
    }
  }
  group.merged.file <- paste0('group_',group_x,'.merged.seurat.Rdata')
  save(group.merged,file=file.path(group.merged.file))
  
}

# TODO: Add a %dopar% loop over the merged files to run seurat analysis steps.
# files <- 
# for (file in files) {
# 
#   x = as.data.table(fread(file))
#   
#   message( 'Creating Seurat object.')
#   x <- CreateSeuratObject(counts = x, project = "IS022")
#   
#   message('Normalizing Data')
#   x <- NormalizeData(x)
#   
#   message('Findin Variable Features')
#   x <- FindVariableFeatures(x)
#   
#   message('Scaling Data')
#   x <- ScaleData(x, features = row.names(as.data.frame(x@assays$RNA@data)))
#   Idents(x) <- 'sample'
#   
#   message( 'Running PCA')
#   x <- RunPCA(x)
#   
#   message( 'Finding Neighbor')
#   x <- FindNeighbors(x)   
#   
#   message( 'Finding Clusters')
#   x <- FindClusters(x)
#   
#   message( 'Running tSNE')
#   x <- RunTSNE(x)
#   
#   message( 'Running UMAP')
#   x <- RunUMAP(x)
#   
#   message( 'Saving')
#   save( x, paste( args[1], '.Seurat', sep = '' ))
#   message( 'Sike')
#   
# }

# TODO: THen continue with some downstream tools like enrichR, etc.
