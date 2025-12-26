# Código R - UMAPs individuales por día

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})
set.seed(1234)

# Carpeta con los archivos .RDS de scRNA-seq
data_dir <- "C:/Users/Usuario/Desktop/GSE TIPO DE LESION/iSNAT aplastamiento crush jovenes 0 1 3 7d"
setwd(data_dir)

# Función para procesar un objeto Seurat (lectura, filtrado, normalización, reducción, clustering, UMAP)
run_umap_simple <- function(file) {
  seu <- readRDS(file)
  if (!inherits(seu, "Seurat")) {
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
      seu <- tryCatch(as.Seurat(seu, counts = "counts", data = "logcounts"),
                      error = function(e) stop("No puedo convertir a Seurat: ", basename(file)))
    } else stop("Objeto no-Seurat y no puedo convertir: ", basename(file))
  }
  # Añadir porcentaje de genes mitocondriales si no existe
  if (!"percent.mt" %in% colnames(seu@meta.data))
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, "^mt-")
  # Filtrado de células por número de genes detectados y porcentaje mitocondrial
  seu <- subset(seu, nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

  # Filtrar genes poco expresados (presentes en <3 células)
  counts <- tryCatch(GetAssayData(seu, assay="RNA", layer="counts"),
                     error = function(e) NULL)
  if (is.null(counts)) {
    counts <- tryCatch(GetAssayData(seu, assay="RNA", slot="counts"),
                       error = function(e) NULL)
  }
  if (!is.null(counts)) {
    keep <- Matrix::rowSums(counts > 0) >= 3
    seu <- subset(seu, features = rownames(counts)[keep])
  }

  # Normalización y reducción de dimensionalidad
  seu <- SCTransform(seu, verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE, npcs = 50)
  seu <- FindNeighbors(seu, dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  return(seu)
}

# Detecta y procesa solo los días de interés (0, 3, 7)
files <- list.files(pattern = "\\.RDS$", full.names = TRUE, ignore.case = TRUE)
get_tp <- function(f) regmatches(basename(f), regexpr("day_\\d+", basename(f)))
tps   <- sub("day_", "", get_tp(files))
sel   <- !is.na(tps) & tps %in% c("0","3","7")
files <- files[sel]; tps <- tps[sel]

objs <- lapply(files, run_umap_simple)
names(objs) <- paste0("day", tps)

# Guardado de objetos procesados y figuras UMAP por día
for (nm in names(objs)) {
  p <- DimPlot(objs[[nm]], label = TRUE, pt.size = 0.6) +
    ggtitle(paste("UMAP", nm)) +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16, face="bold"))
  saveRDS(objs[[nm]], file = file.path(data_dir, paste0("iSNAT_", nm, "_processed_simple.rds")))
  ggsave(file.path(data_dir, paste0("UMAP_", nm, "_clusters_simple.png")),
         p, width = 7, height = 6, dpi = 300)
}
