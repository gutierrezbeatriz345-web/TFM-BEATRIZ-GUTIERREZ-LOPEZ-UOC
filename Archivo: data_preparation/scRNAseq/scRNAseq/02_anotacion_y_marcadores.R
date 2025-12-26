# Código R - Anotación SingleR y marcadores por día

suppressPackageStartupMessages({
  library(Seurat); library(SingleR); library(celldex)
  library(dplyr); library(ggplot2); library(readr)
})
set.seed(1234)

data_dir <- "C:/Users/Usuario/Desktop/GSE TIPO DE LESION/iSNAT aplastamiento crush jovenes 0 1 3 7d"
out_dir  <- file.path(data_dir, "outputs_annotation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Usa los objetos de la sesión si existen; si no, carga los .rds procesados
seu_list <- list()
if (exists("objs")) {
  seu_list <- objs
} else {
  cand <- list.files(data_dir, pattern = "iSNAT_day.*_processed_simple\\.rds$", 
                     full.names = TRUE, ignore.case = TRUE)
  stopifnot(length(cand) > 0)
  get_nm <- function(p) sub("^.*(day\\d+).*", "\\1", basename(p))
  nms <- vapply(cand, get_nm, character(1))
  seu_list <- setNames(lapply(cand, readRDS), nms)
}

# Cargar referencia para SingleR (dataset MouseRNAseqData del paquete celldex)
ref <- celldex::MouseRNAseqData()

# Función para anotar un objeto Seurat y obtener marcadores
annotate_and_markers <- function(seu, name_tag) {
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", 
                       scale.factor = 10000, verbose = FALSE)
  test_mat <- GetAssayData(seu, assay = "RNA", layer = "data")
  ann <- SingleR(test = test_mat, ref = ref, labels = ref$label.main)
  seu$celltype <- ann$labels

  # Guardar UMAP con anotación SingleR
  p_anno <- DimPlot(seu, group.by = "celltype", label = TRUE, repel = TRUE) +
    ggtitle(paste0("UMAP ", name_tag, " · Tipos celulares (SingleR)"))
  ggsave(file.path(out_dir, paste0("UMAP_", name_tag, "_SingleR.png")),
         p_anno, width = 9, height = 6, dpi = 300)

  # Encontrar marcadores de todos los clusters
  DefaultAssay(seu) <- if ("SCT" %in% Assays(seu)) "SCT" else "RNA"
  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, 
                            logfc.threshold = 0.25, verbose = FALSE)
  readr::write_csv(markers, file.path(out_dir, paste0("markers_", name_tag, "_clusters.csv")))

  # Top 5 marcadores por cluster
  top_markers <- markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
  readr::write_csv(top_markers, file.path(out_dir, paste0("markers_", name_tag, "_top5_por_cluster.csv")))

  # Comparación de expresión entre cluster 0 vs 1 (si existen)
  Idents(seu) <- "seurat_clusters"
  clust_levels <- levels(Idents(seu))
  if (all(c("0","1") %in% clust_levels)) {
    de01 <- FindMarkers(seu, ident.1 = 0, ident.2 = 1, 
                        min.pct = 0.25, logfc.threshold = 0.25)
    de01$gene <- rownames(de01)
    de01$significant <- ifelse(de01$p_val_adj < 0.05 & abs(de01$avg_log2FC) > 0.5, "Sí", "No")
    readr::write_csv(de01, file.path(out_dir, paste0("DE_cluster0_vs1_", name_tag, ".csv")))

    p_volc <- ggplot(de01, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
      geom_point(alpha = 0.8) + scale_color_manual(values = c("grey", "red")) +
      theme_minimal() + labs(title = paste0("Volcán: clúster 0 vs 1 (", name_tag, ")"),
                             x = "Log2 cambio de expresión", y = "-log10 p-valor ajustado")
    ggsave(file.path(out_dir, paste0("Volcan_0vs1_", name_tag, ".png")),
           p_volc, width = 7, height = 5, dpi = 300)
  }
  saveRDS(seu, file.path(out_dir, paste0("iSNAT_", name_tag, "_annotated.rds")))
  return(invisible(seu))
}

# Anotar todos los objetos (días) y buscar marcadores
for (nm in names(seu_list)) {
  message("Anotando y buscando marcadores: ", nm)
  seu_list[[nm]] <- annotate_and_markers(seu_list[[nm]], nm)
}

message("Listo. Revisa: ", out_dir)
