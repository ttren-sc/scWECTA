######------- All function of data preprocess -------#######
# Import libraries
library(Seurat)
library(scran)
library(M3Drop)
library(SingleCellExperiment)
library(xfun)


# Preprocess train dataset
trainProcess <- function(counts, meta, path = "/home/rtt/Music/Projects/pyProjects/exp2/data", species = "human", tissue = "pancreas", name) {
  ## Create seurat object
  obj <- pipeline(counts, meta, num = 1000)
  obj@meta.data$ct <- as.factor(meta$label)
  obj@active.ident <- obj@meta.data$ct
  markers <- FindAllMarkers(obj, group.by = "ct", logfc.threshold = 0.5, only.pos = TRUE)

  # Gene selection
  markers_filter <- markers[which(markers$avg_log2FC>1), "gene"]
  obj_hvgs <- obj@assays$RNA@var.features # 2000
  obj_m3drop <- rownames(fsM3drop(obj@assays$RNA@counts))
  obj_scran <- fsScran(obj@assays$RNA@counts, num=1000)
  # genes <- Reduce(union, list(baron_markers_filter, baron_hvgs, baron_m3drop, baron_scran))

  dir_path <- paste0(path, "/", species, "_", tissue, "/train_", species, "_", name)
  if (!dir_exists(dir_path))
    dir.create(dir_path, recursive = TRUE)
  
  write.csv(t(obj@assays$RNA@data), paste0(dir_path, "/train.csv"))
  write.csv(meta, paste0(dir_path, "/train_label.csv"))
  write.csv(data.frame(gene = markers_filter), paste0(dir_path, "/seurat_markers.csv"))
  write.csv(data.frame(gene = obj_hvgs), paste0(dir_path, "/seurat_hvgs.csv"))
  write.csv(data.frame(gene = obj_m3drop), paste0(dir_path, "/m3drop_hvgs.csv"))
  write.csv(data.frame(gene = obj_scran), paste0(dir_path, "/scran_hvgs.csv"))

  # assign(paste0(name, ".obj"), obj)
  obj
}

# Preprocess test dataset
testProcess <- function(counts, meta, path = "/home/rtt/Music/Projects/pyProjects/exp2/data", species = "human", tissue = "pancreas", name = "baron") {
  ## seurat objectname = "xin") {
  ## seurat object
  obj <- pipeline(counts, meta, num = 1000)

  obj@meta.data$ct <- as.factor(meta$label)
  obj@active.ident <- obj@meta.data$ct

  dir_path <- paste0(path, "/", species, "_", tissue, "/test_", species, "_", name)
  if (!dir_exists(dir_path))
    dir.create(dir_path, recursive = TRUE)
  
  write.csv(t(obj@assays$RNA@data), paste0(dir_path, "/test.csv"))
  write.csv(meta, paste0(dir_path, "/test_label.csv"))

  obj
}


train_test_process <- function(train, test, train_meta, test_meta, path = "/home/rtt/Music/Projects/pyProjects/exp2/data/common_genes", species = "human", tissue = "pancreas", train_name = "baron", test_name = "muraro") {
  ## Get common genes of train and test data
  train_genes <- rownames(train)
  test_genes <- rownames(test)
  common_genes <- intersect(train_genes, test_genes)
  train <- train[common_genes, ]
  test <- test[common_genes, ]
  
  ## Create seurat object of train
  obj <- pipeline(train, train_meta, num = 1000)
  obj@meta.data$ct <- as.factor(train_meta$label)
  obj@active.ident <- obj@meta.data$ct
  markers <- FindAllMarkers(obj, group.by = "ct", logfc.threshold = 0.5, only.pos = TRUE)
  # Gene selection
  if (dim(markers)[1] > 2000) {
    markers <- markers[which(markers$avg_log2FC>1), "gene"]
  }
  obj_hvgs <- obj@assays$RNA@var.features # 2000
  obj_m3drop <- rownames(fsM3drop(obj@assays$RNA@counts))
  obj_scran <- fsScran(obj@assays$RNA@counts, num=2000)
  # genes <- Reduce(union, list(baron_markers_filter, baron_hvgs, baron_m3drop, baron_scran))
  
  dir_path <- paste0(path, "/", train, "_pred_", test, "/", species, "_", tissue, "/train_", species, "_", train_name)
  if (!dir_exists(dir_path))
    dir.create(dir_path, recursive = TRUE)
  
  write.csv(t(obj@assays$RNA@data), paste0(dir_path, "/train.csv"))
  write.csv(train_meta, paste0(dir_path, "/train_label.csv"))
  write.csv(data.frame(gene = markers), paste0(dir_path, "/seurat_markers.csv"))
  write.csv(data.frame(gene = obj_hvgs), paste0(dir_path, "/seurat_hvgs.csv"))
  write.csv(data.frame(gene = obj_m3drop), paste0(dir_path, "/m3drop_hvgs.csv"))
  write.csv(data.frame(gene = obj_scran), paste0(dir_path, "/scran_hvgs.csv"))
  
  ## Create seurat object of test
  obj <- pipeline(test, test_meta)
  
  obj@meta.data$ct <- as.factor(test_meta$label)
  obj@active.ident <- obj@meta.data$ct
  
  dir_path <- paste0(path, "/", train, "_pred_", test, "/", species, "_", tissue, "/test_", species, "_", test_name)
  if (!dir_exists(dir_path))
    dir.create(dir_path, recursive = TRUE)
  
  write.csv(t(obj@assays$RNA@data), paste0(dir_path, "/test.csv"))
  write.csv(test_meta, paste0(dir_path, "/test_label.csv"))
  obj
}

# Filter genes and cells for scRNA-seq
filter_GC <- function(counts, meta) {
  if (class(meta) != "data.frame") {
    if (is.null(names(meta)))
      names(meta) <- colnames(counts)
    meta <- data.frame(label = meta)
  }
  
  filter <- counts[which(rowSums(counts) > 0), which(colSums(counts) >0)]
  
  pheno_filter <- data.frame(label = meta[colnames(filter), "label"], row.names = colnames(filter))

  return(c(list(filter), list(pheno_filter)))
}

# Gene selection by M3Drop
fsM3drop <- function(mat) {
  library(M3Drop)
  # M3Drop_DE <- BrenneckeGetVariableGenes(mat)
  norm <- M3DropConvertData(mat, is.counts = TRUE)
  M3Drop_DE <- M3DropFeatureSelection(norm, mt_method = 'fdr', mt_threshold = 0.05, suppress.plot = FALSE)
  return(M3Drop_DE)
}

# Gene selection by scran
fsScran <- function(count, num=1000) {
  library(scran)
  # fit <- fitTrendVar(sce, parametric = TRUE)
  # dec <- decomposeVar(sce, fit)
  library(SingleCellExperiment)
  library(scuttle)
  sce <- SingleCellExperiment(assays = list(counts = count))
  sce <- logNormCounts(sce)
  dec <- modelGeneVar(sce)
  hvgs <- dec[order(dec$bio, decreasing = TRUE), ]
  return(rownames(hvgs)[1:num])
}

# Data pre-process and gene selection by seurat
pipeline <- function(count, meta, num = 1000) {
  require(Seurat)
  sce <- CreateSeuratObject(count, meta.data = meta)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize")
  sce <- FindVariableFeatures(sce, nfeatures = num) # 2000
  sce <- ScaleData(sce)
  sce <- RunPCA(sce)
  sce <- FindNeighbors(sce)
  sce <- FindClusters(sce)
  sce <- RunUMAP(sce, dims = 1:50)
  sce <- RunTSNE(sce, dims = 1:50)
  sce
}

