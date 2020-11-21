library(BUSpaRse)
library(DropletUtils)
library(Matrix)
library(tidyverse)
library(Seurat)
library(ggpointdensity)
library(scico)
library(scales)
library(useful)
library(dplyr)
library(ggplot2)
library(stringr)
library(parallel)
library(RColorBrewer)
library(ape)
library(mltools)
theme_set(theme_bw())
set.seed(1)

YourMatrixObject_initial <- read.csv(file = "/home/ubuntu/habenula_scrnaseq_vtp/data/countdata/hab161103.1.csv", header = T, row.names = 1)

intermediate_matrix_1 <- data.table(YourMatrixObject_initial)

rownames(intermediate_matrix_1)<-rownames(YourMatrixObject_initial)
colnames(intermediate_matrix_1)<-colnames(YourMatrixObject_initial)
intermediate_matrix_2<-sparsify(intermediate_matrix_1, sparsifyNAs = FALSE, naCols = "none")
rownames(intermediate_matrix_2)<-rownames(intermediate_matrix_1)

YourMatrixObject <- intermediate_matrix_2

tot_counts <- colSums(YourMatrixObject)
lib_sat <- tibble(nCount = tot_counts,
                  nGene = colSums(YourMatrixObject > 0))

options(repr.plot.width=9, repr.plot.height=6)
ggplot(lib_sat, aes(nCount, nGene)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_x_log10() + scale_y_log10() + annotation_logticks()

ggplot(lib_sat, aes(nCount, nGene)) +
  geom_bin2d(bins = 50) +
  scale_fill_scico(palette = "devon", direction = -1, end = 0.95) +
  scale_x_log10() + scale_y_log10() + annotation_logticks()

bc_rank <- barcodeRanks(YourMatrixObject, lower = 10)

knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>%
    distinct() %>%
    dplyr::filter(total > 0) # nothing has total > 0 here
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}

options(repr.plot.width=9, repr.plot.height=6)

knee_plot(bc_rank)

dim(YourMatrixObject)


YourMatrixObject <- YourMatrixObject[, tot_counts > metadata(bc_rank)$inflection]
YourMatrixObject <- YourMatrixObject[Matrix::rowSums(YourMatrixObject) > 0, ]

dim(YourMatrixObject) 

Your_Seurat_Object <- CreateSeuratObject(counts = YourMatrixObject, project = "habenula", min.cells = 3, min.features = 200) 

Your_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Your_Seurat_Object, pattern = "^mt.") 

head(Your_Seurat_Object) 

VlnPlot(Your_Seurat_Object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Your_Seurat_Object, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(Your_Seurat_Object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

CombinePlots(plots = list(plot1, plot2))

#cp3 starts

Your_Seurat_Object <- subset(Your_Seurat_Object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 500 & nCount_RNA < 18000)  

Your_Seurat_Object <- NormalizeData(Your_Seurat_Object, normalization.method = "LogNormalize", scale.factor = 10000)

options(repr.plot.width=9, repr.plot.height=6)
Your_Seurat_Object <- FindVariableFeatures(Your_Seurat_Object, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(Your_Seurat_Object), 10)

# Plot variable features with and without labels

plot1 <- VariableFeaturePlot(Your_Seurat_Object, log = FALSE)

LabelPoints(plot = plot1, points = top10, repel = TRUE)

#CP4 STARTS

Your_Seurat_Object <- ScaleData(Your_Seurat_Object)

all.genes <- rownames(Your_Seurat_Object)

Your_Seurat_Object <- ScaleData(Your_Seurat_Object, features = all.genes)

Your_Seurat_Object <- RunPCA(Your_Seurat_Object, features = VariableFeatures(object = Your_Seurat_Object))
ElbowPlot(Your_Seurat_Object)

print(Your_Seurat_Object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Your_Seurat_Object, dims = 1:2, reduction = "pca")

DimPlot(Your_Seurat_Object, reduction = "pca")

DimHeatmap(Your_Seurat_Object, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(Your_Seurat_Object, dims = 1:15, cells = 500)

options(repr.plot.width=6, repr.plot.height=8)
FeaturePlot(Your_Seurat_Object, reduction = "pca", features = "Clu")

#CP5 Starts

Your_Seurat_Object <- JackStraw(Your_Seurat_Object, num.replicate = 100)

Your_Seurat_Object <- ScoreJackStraw(Your_Seurat_Object, dims = 1:20)

JackStrawPlot(Your_Seurat_Object, dims = 1:20)  # plotting the JackStraw

Your_Seurat_Object <- FindNeighbors(Your_Seurat_Object, dims = 1:20)

Your_Seurat_Object <- FindClusters(Your_Seurat_Object, resolution = 0.15) 

head(Idents(Your_Seurat_Object), 5)

Your_Seurat_Object<-RunTSNE(Your_Seurat_Object, dims = 1:20)

TSNEPlot(Your_Seurat_Object, label = TRUE) + xlim(-50,50) + ylim(-50,50)

#CP6 Starts

cluster1.markers <- FindMarkers(Your_Seurat_Object, ident.1 = 1, min.pct = 0.1, logfc.threshold = 0.25)

head(Your_Seurat_Object, n = 5)

cluster5.markers <- FindMarkers(Your_Seurat_Object, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.1, logfc.threshold = 0.25)

head(cluster5.markers, n = 5)

Your_Seurat_Object.markers <- FindAllMarkers(Your_Seurat_Object, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

Your_Seurat_Object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

Your_Seurat_Object.markers_table<-data.frame(Your_Seurat_Object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC))

cluster1.markers <- FindMarkers(Your_Seurat_Object, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#You can plot the the violing plot distribution of any gene by cluster. N.B. these genes are placeholders, feel free investigate your genes of choice in addiation to these.  
VlnPlot(Your_Seurat_Object, features = c("Tac2", "Gap43"))

VlnPlot(Your_Seurat_Object, features = c("F2r", "Kcnma1"), slot = "counts", log = TRUE)

#CP7 Starts

FeaturePlot(Your_Seurat_Object, features = c("Gabbr1"))

FeaturePlot(Your_Seurat_Object, features = c("Cxx1a", "Fos", "Cck", "Htr2c", "Cldn11", "Hpca", "Necab1", "Epas1", "C1qc", "Apoe"))

top10_markers <- Your_Seurat_Object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(Your_Seurat_Object, features = top10_markers$gene) + NoLegend()

#CP8

new.cluster.ids <- c("Cxx1a", "Fos", "Cck", "Htr2c", "Plp1", "Hpca", "Necab1")

names(new.cluster.ids) <- levels(Your_Seurat_Object)

Your_Seurat_Object <- RenameIdents(Your_Seurat_Object, new.cluster.ids) 
DimPlot(Your_Seurat_Object, label = TRUE, pt.size = 0.5, repel = T) + NoLegend() 

Idents(Your_Seurat_Object) <-  Your_Seurat_Object@meta.data$RNA_snn_res.0.15 #your res value can be changed

#DotPlot for candidate genes by cluster:
DotPlot(Your_Seurat_Object, features = rev(c("Cxx1a", "Fos", "Cck", "Htr2c", "Cldn11", "Hpca", "Necab1", "Epas1", "Apoe"))) + RotatedAxis()

DotPlot(Your_Seurat_Object, features = rev(new.cluster.ids)) + RotatedAxis()