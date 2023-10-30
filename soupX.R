" Goal: Cleanup source data for etx scRNAseq 
Date:230129
Author:Carsten Knutsen
"
#install.packages('SoupX')
#install.packages('installr')
#library(installr)
#install.Rtools()
#install.packages('Matrix')
#BiocManager::install("DropletUtils")
library(DropletUtils)
library(Matrix)
library(SoupX)
library(Seurat)
data_dir <- '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/CellRanger_output'
output_dir <-'/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/soupx'
sub_fols <- list.dirs(path = data_dir, full.names = FALSE, recursive = FALSE)
for (fol in sub_fols)
{
  print(fol)
  subdir <- sprintf('%s/%s/outs',data_dir,fol)
  subdir_out <- sprintf('%s/%s/outs',output_dir,fol)
  dir.create(subdir_out, recursive = TRUE,showWarnings = FALSE)
  cellnames <- read.csv(sprintf('%s/filtered_feature_bc_matrix/barcodes.tsv.gz', subdir),header =FALSE)
  filt.matrix <- Read10X_h5(sprintf("%s/filtered_feature_bc_matrix.h5",subdir),use.names = F)
  raw.matrix <- Read10X_h5(sprintf("%s/raw_feature_bc_matrix.h5",subdir),use.names = F)
  soup.channel = SoupChannel(raw.matrix, filt.matrix)
  srat <- CreateSeuratObject(counts = filt.matrix)
  srat <-RenameCells(srat, new.names = cellnames$V1)
  srat    <- SCTransform(srat, verbose = F)
  srat    <- RunPCA(srat, verbose = F)
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat    <- FindClusters(srat, verbose = T)
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel  <- setDR(soup.channel, umap)
  soup.channel  <- autoEstCont(soup.channel)
  adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
  DropletUtils:::write10xCounts(sprintf("%s/soupx_filt",subdir_out), adj.matrix)
  
}

