# Código R - Renombrar clusters y generar UMAPs nombrados

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)

# Día 0: renombrar clusters del objeto objs$day0
Idents(objs$day0) <- factor(objs$day0$seurat_clusters)
new.cluster.ids.day0 <- c(
  "Fibro_eMES", "Neurona_DRG", "Endotelio_BNB", "Mieloide_resid", "SC_mieliniz",
  "Epitelio_contam", "APC_mixta", "Fibro_CXCL", "Eritro", "Mesotelio",
  "Prog_mesen", "vSMC", "Pericito", "Neurona_noci", "vSMC_alt", "cDC_conv",
  "EC_linfa", "Masto", "Linf_T", "EC_fenestr", "Epitelio_simp", "Prolif_S"
)
names(new.cluster.ids.day0) <- levels(objs$day0)
objs$day0 <- RenameIdents(objs$day0, new.cluster.ids.day0)

# Día 3: renombrar clusters del objeto objs$day3
Idents(objs$day3) <- factor(objs$day3$seurat_clusters)
new.cluster.ids.day3 <- c(
  "Mono_inflam", "Fibro_dMES", "Mac_GPNMB", "Fibro_CXCL", "SC_satGFAP",
  "cDC_MoDC", "Mac_inflam", "Epitelio_contam", "EC_angio", "Fibro_LRRC15",
  "SC_repar", "vSMC", "Mesotelio_contam", "cDC_madura", "NK",
  "Neutro", "Prolif_hist", "cDC1", "Masto"
)
names(new.cluster.ids.day3) <- levels(objs$day3)
objs$day3 <- RenameIdents(objs$day3, new.cluster.ids.day3)

# Día 7: renombrar clusters del objeto objs$day7
Idents(objs$day7) <- factor(objs$day7$seurat_clusters)
new.cluster.ids.day7 <- c(
  "Fibro_ECM", "Mac_GPNMB_B", "Mac_FOLR2", "SC_repar", "cDC_CD209a", "Linf_T_act",
  "SC_Sox10", "Perivas_ECM", "Epitelio_contam", "NK", "Axon_TdT", "Mac_infla",
  "EC_ROBO4", "Fibro_COL2", "Prog_DLk1", "cDC2_MoDC", "cDC1_prolif",
  "Pericito_Notch3", "Prolif_mito", "Mesot_MSLN", "vSMC_OLFR", "Cel_B",
  "Mac_alt", "Masto", "Mieloide_IFN"
)
names(new.cluster.ids.day7) <- levels(objs$day7)
objs$day7 <- RenameIdents(objs$day7, new.cluster.ids.day7)

# Generar UMAPs con clusters renombrados para cada día
p0 <- DimPlot(objs$day0, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP day0") +
  theme(plot.title = element_text(size = 16, face = "bold"))
p3 <- DimPlot(objs$day3, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP day3") +
  theme(plot.title = element_text(size = 16, face = "bold"))
p7 <- DimPlot(objs$day7, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP day7") +
  theme(plot.title = element_text(size = 16, face = "bold"))

# Guardar figuras individuales
dir.create("outputs_figs", showWarnings = FALSE)
ggsave("outputs_figs/UMAP_day0_named.png", p0, width = 9, height = 7, dpi = 300)
ggsave("outputs_figs/UMAP_day3_named.png", p3, width = 9, height = 7, dpi = 300)
ggsave("outputs_figs/UMAP_day7_named.png", p7, width = 9, height = 7, dpi = 300)

# Panel combinado de UMAPs (día 0 y 3 arriba, día 7 abajo)
panel <- (p0 | p3) / p7
ggsave("outputs_figs/UMAP_panel_named.png", panel, width = 15, height = 10, dpi = 300)
