
###########################
## load required packages
###########################
library(Seurat)
library(ggplot2)

## colors will be used 
colors<-c("#4d648d", "#0450fb", "#11aac4", "#42e8f3", "#AEC7E8", "#2CA02C", "#98DF8A", "#9eccaf", "#daf400", "#983b59", "#e81f3f", "#ff8b94", "#ffd3b6", "#f9ae34", "#ffdb00", "#723584", "#9264eb", "#ff00ff", "#E377C2", "#de94e4", "#F7B6D2", "#C5B0D5", "#8C564B", "#C49C94", "#BCBD22", "#DBDB8D", "#7F7F7F", "#C7C7C7", "#a7a7b7")

###########################
## Data setup
###########################
## This is an intrgrated comparative analysis of two single cell datasets using Seurat development version 3.0.0.9000 [@satija2015spatial; @butler2018integrating; @hafemeister2019normalization].
## First read in the scRNA-seq data of the datasets and set up Seurat objects (this includes filtering cells/genes in each dataset). 

## read in the 10X data
## first read in the sample information sheet
## need to specify the count matrices location generated from cellRanger in the sample information file
sampleInfo<-read.table("Sample_info.txt", header=T, sep="\t") ## or set path to Sample_info.txt
  
## read in data
scRNA.list<-list()
for (i in 1:nrow(sampleInfo)) {
  ## read in 10X data
  scData<-Read10X(data.dir = as.character(sampleInfo[i,3]))
  ## change cell names
  scData@Dimnames[2][[1]]<-paste0(as.character(sampleInfo[i,2]), "_", scData@Dimnames[2][[1]])
  ## create Seurat object
  scRNA <- CreateSeuratObject(scData, min.cells = 5, min.features = 200, project = as.character(sampleInfo[i,1]))
  scRNA[['sample']] <- as.character(sampleInfo[i,2]) ## or scRNA@meta.data$sample<-sample
  ## mitochondrial genes
  ref<-as.character(sampleInfo[i,4])
  percent.mito <- Matrix::colSums(GetAssayData(scRNA, slot = 'counts')[grep(pattern = "^MT-", rownames(scRNA), value = TRUE), ]) / Matrix::colSums(GetAssayData(scRNA, slot = 'counts'))
  scRNA[['percent.mito']] <- percent.mito
  
  ## filter data
  scRNA <- subset(scRNA, subset = nFeature_RNA > sampleInfo[i,5] & nFeature_RNA < sampleInfo[i,6] & percent.mito < sampleInfo[i,7])
  ## assign data to data list
  scRNA.list[[i]]<-scRNA
  rm(scRNA, scData)
}
## add name to the data list
names(scRNA.list)<-as.character(sampleInfo[,2])

## Data normalization and variable features identification
for (i in 1:length(scRNA.list)) {
    scRNA.list[[i]] <- NormalizeData(scRNA.list[[i]], verbose = FALSE)
    scRNA.list[[i]] <- FindVariableFeatures(scRNA.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

###########################
# Data Integration
###########################
## Find integration anchors and integrate data
scRNA.anchors <- FindIntegrationAnchors(scRNA.list, dims = 1:30)
scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)

## switch to integrated assay
## automatically set during IntegrateData
DefaultAssay(object = scRNA.integrated) <- "integrated"

## Run the standard workflow of data scaling, dimension reduction
scRNA.integrated <- ScaleData(object = scRNA.integrated, verbose = FALSE)
scRNA.integrated <- RunPCA(object = scRNA.integrated, npcs = 50, verbose = FALSE)
scRNA.integrated <- RunUMAP(object = scRNA.integrated, reduction = "pca", dims = 1:12)
scRNA.integrated <- RunTSNE(object = scRNA.integrated, reduction = "pca", dims = 1:12)

## Graph-based clustering
scRNA.integrated <- FindNeighbors(scRNA.integrated, dims = 1:12, verbose=FALSE)
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.8)

## T-SNE plots
DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.8", label = TRUE, repel = TRUE, cols=colors[1:length(unique(scRNA.integrated@active.ident))]) + NoLegend()
DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.8", label = FALSE, cols=colors[1:length(unique(scRNA.integrated@active.ident))])

## UMAP plots
DimPlot(scRNA.integrated, reduction = "umap", group.by = "sample", cols=c("grey", "firebrick1"))
DimPlot(scRNA.integrated, reduction = "umap", group.by = "integrated_snn_res.0.8", label = TRUE, cols=colors[1:length(unique(scRNA.integrated@active.ident))]) + NoLegend()
DimPlot(scRNA.integrated, reduction = "umap", group.by = "integrated_snn_res.0.8", label = FALSE, cols=colors[1:length(unique(scRNA.integrated@active.ident))])

## split by sample
DimPlot(scRNA.integrated, reduction = "umap", split.by = "sample", label = TRUE, combine=FALSE, cols=colors[1:length(unique(scRNA.integrated@active.ident))]) + NoLegend()

## Identify conserved cell type markers

DefaultAssay(scRNA.integrated) <- "RNA"

## get the cell clusters
uniqueC<-unique(scRNA.integrated@active.ident)[order(unique(scRNA.integrated@active.ident))]
## if no clusters present in only one condition
ct.markers<-lapply(uniqueC, function(i) FindConservedMarkers(scRNA.integrated, ident.1 = i, assay.type = "RNA", grouping.var = "sample"))

## violin plots for the top cluster marker genes or specified gene list
## geneList
vlnP<-lapply(seq_along(uniqueC), function(i) VlnPlot(scRNA.integrated, features= head(row.names(ct.markers[[i]]), 9), ncol=3, cols=colors[1:length(unique(scRNA.integrated@active.ident))]) + theme(legend.position="none"))

## save the integrated Seurat object
saveRDS(scRNA.integrated, file = "InformativeFileName.rds")

## DE genes within a cluster between conditions

