# Código R - Análisis de expresión génica (bulk RNA-seq: jóvenes vs viejos, cut vs UI)

# =========================================================
# Bulk RNA-seq  (Mouse Ensembl 112; Salmon quantification)
# Análisis integral: Infiltración inmune + Tregs + Senescencia/Resolución + Vías
# QC + DESeq2 + HK automáticos + Heatmaps/Log2FC + Scores + Excel
# Orden de grupos fijo: Young_UI, Young_Cut3dpi, Old_UI, Old_Cut3dpi
# =========================================================

## 0) Cargar paquetes necesarios (e instalar si falta alguno)
cran_pkgs <- c(
  "data.table","dplyr","tidyr","ggplot2","readr","stringr","purrr",
  "forcats","openxlsx","jsonlite","ggrepel","pheatmap","writexl",
  "matrixStats","tibble"
)
to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bio_pkgs <- c("DESeq2","tximport","ComplexHeatmap","org.Mm.eg.db","AnnotationDbi")
bio_missing <- bio_pkgs[!sapply(bio_pkgs, requireNamespace, quietly=TRUE)]
if (length(bio_missing)) BiocManager::install(bio_missing, ask=FALSE, update=FALSE)

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(stringr); library(purrr); library(forcats); library(openxlsx); library(jsonlite)
  library(ggrepel); library(pheatmap); library(DESeq2); library(tximport); library(matrixStats)
  library(ComplexHeatmap); library(org.Mm.eg.db); library(AnnotationDbi); library(tibble)
})

## Helpers (funciones auxiliares)
normalize_names <- function(v) { 
  v <- as.character(v); 
  v[is.na(v)] <- ""; 
  gsub("\\.\\d+$","", v) %>% toupper() 
}
pick_first_existing <- function(paths, is_dir=TRUE) {
  ok <- if (is_dir) vapply(paths, dir.exists, logical(1)) else vapply(paths, file.exists, logical(1))
  if (!any(ok)) stop("Ninguna ruta existe:\n", paste(paths, collapse="\n"))
  paths[which(ok)[1]]
}

## 1) Rutas de datos (ajustar si cambias carpeta de trabajo)
root <- pick_first_existing(c(
  "C:/Users/Usuario/Desktop/JOSE-PR0022-41335294/quant_all_merged",
  "C:/Users/Usuario/Desktop/INFLAMOSOMA/JOSE-PR0022-41335294/quant_all_merged"
), is_dir = TRUE)
gtf  <- pick_first_existing(c(
  "C:/Users/Usuario/Desktop/JOSE-PR0022-41335294/mouse_Ensembl112/Mus_musculus.GRCm39.112.gtf",
  "C:/Users/Usuario/Desktop/INFLAMOSOMA/JOSE-PR0022-41335294/mouse_Ensembl112/Mus_musculus.GRCm39.112.gtf"
), is_dir = FALSE)

fig_dir <- file.path(root, "figs_inmunidad");    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
res_dir <- file.path(root, "results_inmunidad"); dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
outdir  <- file.path(root, "_tximport_outputs"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
message("root: ", root); message("GTF : ", gtf)

## 2) Detectar muestras con resultados de Salmon (archivos quant.sf)
dirs_all <- list.dirs(root, full.names = TRUE, recursive = FALSE)
dirs     <- dirs_all[file.exists(file.path(dirs_all, "quant.sf"))]
if (!length(dirs)) stop("No hay quant.sf en: ", root)
samples <- basename(dirs)
files   <- setNames(file.path(dirs, "quant.sf"), samples)
message("Muestras: ", paste(samples, collapse = ", "))

## 3) Crear tabla de metadatos (automática + manual)
sample_table <- tibble::tibble(sample = samples) %>%
  mutate(Age = case_when(grepl("^Old|^UNOLD", sample, ignore.case = TRUE) ~ "Old", TRUE ~ "Young"),
         Injury = ifelse(grepl("^UN", sample, ignore.case = TRUE), "UI", "Cut3dpi"),
         Group = paste(Age, Injury, sep = "_"))

# Información manual adicional de muestras (sexo, etiquetas originales, etc.)
manual_meta <- tibble::tribble(
 ~sample,          ~NM,     ~WT,  ~Gender,   ~Injury_label,
 "OldWT2_S10",     "NM0177","WT2","Male",    "Cut 3dpi",
 "OldWT3_S11",     "NM0178","WT3","Female",  "Cut 3dpi",
 "OldWT4_S12",     "NM0179","WT4","Male",    "Cut 3dpi",
 "UNOLDWT2_S16",   "NM0183","WT2","Male",    "UI",
 "UNOLDWT3_S17",   "NM0184","WT3","Male",    "UI",
 "UNOLDWT4_S18",   "NM0185","WT4","Mix",     "UI",
 "UNYoungWT1_S19", "NM0186","WT1","Male",    "UI",
 "UNYoungWT2_S20", "NM0187","WT2","Female",  "UI",
 "UNYoungWT3_S21", "NM0188","WT3","Male",    "UI",
 "YoungWT1_S1",    "NM0168","WT1","Male",    "Cut 3dpi",
 "YoungWT2_S2",    "NM0169","WT2","Male",    "Cut 3dpi",
 "YoungWT3_S3",    "NM0170","WT3","Female",  "Cut 3dpi",
 "YoungWT4_S4",    "NM0171","WT4","Female",  "Cut 3dpi"
)

# Integrar metadatos manuales en la tabla principal
sample_table <- sample_table %>%
  left_join(manual_meta, by = "sample") %>%
  mutate(Gender = if_else(is.na(Gender), "NA", Gender),
         WT = if_else(is.na(WT), "NA", WT),
         NM = if_else(is.na(NM), "NA", NM),
         Group = factor(Group, levels = c("Young_UI","Young_Cut3dpi","Old_UI","Old_Cut3dpi")),
         Age = factor(Age, levels = c("Young","Old")),
         Injury = factor(Injury, levels = c("UI","Cut3dpi")),
         quant_path = files[sample],
         meta_path = file.path(root, sample, "aux_info", "meta_info.json"))

## 4) Generar mapa de transcriptos a genes (tx2gene) desde el GTF
gtf_dt <- data.table::fread(
  gtf, sep="\t", header=FALSE, data.table=TRUE, showProgress=TRUE,
  col.names = c("seqname","source","feature","start","end","score","strand","frame","attr")
)
get_attr <- function(x, key) {
  m <- regexpr(paste0(key, ' "([^"]+)"'), x)
  o <- rep(NA_character_, length(x))
  ok <- m > 0
  o[ok] <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x[ok])
  return(o)
}
tx2gene_full <- gtf_dt[feature == "transcript", .(
  transcript_id = get_attr(attr,"transcript_id"),
  gene_id       = get_attr(attr,"gene_id"),
  gene_name     = get_attr(attr,"gene_name")
)]
tx2gene_full <- unique(tx2gene_full[!is.na(transcript_id) & !is.na(gene_id)])
tx2gene_tx_gid <- as.data.frame(tx2gene_full[, .(transcript_id, gene_id)])
gene_map <- unique(as.data.frame(tx2gene_full[, .(gene_id, gene_name)]))

## 5) Importar conteos con tximport (agregando al nivel génico)
txi <- tximport(files = files, type = "salmon", tx2gene = tx2gene_tx_gid,
                countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

## 6) Calcular % de lecturas mapeadas (Salmon)
map_stats <- sample_table %>%
  mutate(percent_mapped = purrr::map_dbl(meta_path, ~ {
    if (file.exists(.x)) {
      mi <- jsonlite::read_json(.x, simplifyVector = TRUE)
      if (!is.null(mi$percent_mapped)) mi$percent_mapped else NA_real_
    } else {
      NA_real_
    }
  })) %>%
  select(sample, Age, Injury, Group, percent_mapped)

## 7) Crear objeto DESeq2 y ejecutar análisis diferencial
coldata <- sample_table %>%
  select(sample, Age, Injury, Group) %>%
  tibble::column_to_rownames("sample")
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ Group)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

## 8) Transformación VST y preparación de matriz de expresión
vst_mat <- assay(vst(dds, blind = TRUE))
gene_symbol <- gene_map$gene_name[match(rownames(vst_mat), gene_map$gene_id)]
gene_symbol[is.na(gene_symbol)] <- rownames(vst_mat)[is.na(gene_symbol)]
rownames(vst_mat) <- make.unique(gene_symbol)
ord_groups <- c("Young_UI","Young_Cut3dpi","Old_UI","Old_Cut3dpi")
ord <- order(match(coldata$Group, ord_groups))
vst_mat <- vst_mat[, ord]; coldata <- coldata[ord, , drop = FALSE]

## 9) Gráficos de control de calidad (opcionales)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Genes detectados por muestra
genes_gt0 <- colSums(counts(dds) > 0) %>% 
  tibble::enframe(name = "sample", value = "genes_gt0") %>% 
  left_join(sample_table, by = "sample")
ggsave(file.path(fig_dir, "genes_detectados_por_muestra.png"),
  ggplot(genes_gt0, aes(x = reorder(sample, genes_gt0), y = genes_gt0, fill = Age)) +
    geom_col() + coord_flip() + facet_wrap(~Injury, nrow = 1, scales = "free_x") +
    labs(x = "Muestra", y = "Genes con conteo > 0", title = "Genes detectados por muestra") + theme_bw(),
  width = 12, height = 6, dpi = 300)

# Tamaño de librería por muestra
lib_size <- colSums(round(txi$counts)) %>% 
  tibble::enframe(name = "sample", value = "assigned_reads") %>% 
  left_join(sample_table, by = "sample")
ggsave(file.path(fig_dir, "tamano_libreria_por_muestra.png"),
  ggplot(lib_size, aes(x = reorder(sample, assigned_reads), y = assigned_reads, fill = Age)) +
    geom_col() + coord_flip() + facet_wrap(~Injury, nrow = 1, scales = "free_x") +
    labs(x = "Muestra", y = "Lecturas asignadas (genes)", title = "Tamaño de librería por muestra") + theme_bw(),
  width = 12, height = 6, dpi = 300)

# Porcentaje de lecturas mapeadas por muestra
ggsave(file.path(fig_dir, "porcentaje_lecturas_mapeadas.png"),
  ggplot(map_stats, aes(x = reorder(sample, percent_mapped), y = percent_mapped, fill = Age)) +
    geom_col() + coord_flip() + facet_wrap(~Injury, nrow = 1, scales = "free_x") +
    labs(x = "Muestra", y = "% mapeado (Salmon)", title = "Porcentaje de lecturas mapeadas") + theme_bw(),
  width = 12, height = 6, dpi = 300)

# PCA de las muestras (usando matriz VST)
pca <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
pca_df <- as_tibble(pca$x[,1:2], rownames = "sample") %>% 
  left_join(as_tibble(coldata, rownames = "sample"), by = "sample")
var_expl <- round(100 * summary(pca)$importance[2,1:2], 1)
ggsave(file.path(fig_dir, "PCA_VST.png"),
  ggplot(pca_df, aes(PC1, PC2, color = Injury, shape = Age, label = sample)) +
    geom_point(size = 3) + ggrepel::geom_text_repel(show.legend = FALSE) +
    labs(title = "PCA  VST", x = paste0("PC1 (", var_expl[1], "%)"), y = paste0("PC2 (", var_expl[2], "%)")) + theme_bw(),
  width = 7, height = 5, dpi = 300)

## 10) Definir módulos de genes (nuevos, según preguntas de investigación)
# Q1: Infiltración inmunitaria
Q1 <- list(
  Myeloid_core    = c("Lyz2","Itgam","Adgre1","Cx3cr1","Cd68","Ms4a7","Tyrobp","Csf1r"),
  Macrophage_M1   = c("Nos2","Ptgs2","Tnf","Il1b","Cxcl9","Cxcl10","Stat1","Irf5","Il12b"),
  Macrophage_M2   = c("Mrc1","Arg1","Mrc2","Retnla","Chil3","Cd163","Maf","Il10","Tgfb1"),
  Neutrophil      = c("S100a8","S100a9","Ly6g","Mpo","Elane","Lcn2","Cxcr2","Csf3r"),
  Monocyte        = c("Ly6c2","Ccr2","Sell","Spi1","Irf8","Itgal"),
  Dendritic       = c("Itgax","Xcr1","Batf3","Zbtb46","Flt3","Cd74","Ciita"),
  Tcell_CD4       = c("Cd3e","Cd4","Trac","Icos","Il7r"),
  Tcell_CD8       = c("Cd3e","Cd8a","Gzmb","Prf1","Ifng"),
  NK              = c("Ncr1","Klrb1c","Eomes","Gzma","Gzmb","Prf1"),
  Bcell           = c("Cd19","Ms4a1","Cd79a","Cd79b","Pax5","Ighm")
)

# Q2: Senescencia y resolución inflamatoria
Q2 <- list(
  Senescence_core = c("Cdkn2a","Cdkn1a","Trp53","Trp63","Trp73","Gadd45a","Gadd45b","Rbl1","Rbl2"),
  SASP            = c("Il6","Il1a","Il1b","Tnf","Cxcl1","Cxcl2","Ccl2","Ccl7","Serpine1","Mmp3","Mmp12","Mmp13","Igfbp3","Igfbp5","Fn1","Thbs1","Lox","Vegfa","Timp1"),
  Resolution      = c("Il10","Il10ra","Il10rb","Il1rn","Nfe2l2","Hmox1","Gpx1","Gpx4","Prdx1","Prdx6","Cat","Sod1","Sod2","Pparg","Klf4","Mafb",
                      "Nr4a1","Nr4a2","Nr4a3","Stat6","Tgfb1","Smad2","Smad3","Smad7","Ppard","Ppargc1a","Rora","Atg5","Atg7","Becn1","Tfeb","Map1lc3b",
                      "Mertk","Axl","Gas6"),
  Efferocytosis   = c("Mertk","Axl","Tyro3","Gas6","Pros1","Itgav","Itgb3","Mfge8")
)

# Q3: Genes asociados a células Treg (comparativa joven vs vieja, lesionado vs UI)
Q3 <- list(
  Treg_core                  = c("Foxp3","Il2ra","Ctla4","Ikzf2","Ikzf4","Nrp1","Tnfrsf18","Tnfrsf4","Tnfrsf9","Il7r"),
  Treg_suppressive           = c("Il10","Tgfb1","Lag3","Pdcd1","Tigit","Havcr2","Entpd1","Nt5e","Ebi3","Il12a"),
  Treg_repair                = c("Areg","Il1rl1","Pparg","Bmp7","Vegfa","Spp1"),
  Treg_trafficking           = c("Ccr7","Cxcr3","Sell","Itgae","Cxcr5","S1pr1"),
  Treg_activation_exhaustion = c("Prdm1","Batf","Irf4","Gata3","Bcl6","Nr4a1","Klrg1","Gzmb","Ms4a4b"),
  Treg_senescent             = c("Cdkn2a","Cdkn1a","Trp53","Klrg1","Tnfrsf1a","Tnfrsf1b")
)

# Q4: Rutas de señalización (ej. IL-2, TGF-β, dianas FoxP3, metabolismo)
Q4 <- list(
  IL2_signaling   = c("Il2ra","Il2rb","Il2rg","Jak1","Jak3","Stat5a","Stat5b","Ptpn2","Socs1","Socs3"),
  TGFb_signaling  = c("Tgfb1","Tgfbr1","Tgfbr2","Smad2","Smad3","Smad4","Smad7","Skil"),
  Foxp3_targets   = c("Foxp3","Ctla4","Il2ra","Ikzf2","Ikzf4","Nrp1","Tnfrsf18","Tnfrsf4","Tnfrsf9"),
  OXPHOS          = c("Ppargc1a","Ppargc1b","Cox4i1","Cox7a2","Ndufa9","Ndufs3","Sdha","Sdhb","Atp5f1a","Atp5f1b","Uqcrc1","Uqcrc2"),
  Glycolysis      = c("Slc2a1","Hk2","Pfkfb3","Pkm","Ldha")
)

# Q5: Subpoblaciones de células Treg
Q5 <- list(
  Treg_suppressor = Q3$Treg_suppressive,
  Treg_repair     = Q3$Treg_repair,
  Treg_effector   = c("Prdm1","Batf","Irf4","Gata3","Tnfrsf9","Tnfrsf4","Klrg1"),
  Treg_central    = c("Sell","Ccr7","Lef1","Tcf7","Il7r"),
  Treg_senescent  = Q3$Treg_senescent
)

MODULES <- list(
  Q1_Infiltracion      = Q1,
  Q2_SenescenciaResol  = Q2,
  Q3_TregEdad          = Q3,
  Q4_Rutas             = Q4,
  Q5_SubpoblacionesTreg = Q5
)

## 11) Comparaciones de expresión diferencial (DESeq2) para contrastes de interés
res_youngCut_vs_youngUI <- results(dds, contrast = c("Group","Young_Cut3dpi","Young_UI"))
res_oldCut_vs_oldUI     <- results(dds, contrast = c("Group","Old_Cut3dpi","Old_UI"))
res_oldCut_vs_youngCut  <- results(dds, contrast = c("Group","Old_Cut3dpi","Young_Cut3dpi"))

## 12) Identificar genes housekeeping (HK) automáticamente
tpm_mat <- as.data.frame(txi$abundance)
tpm_mat$gene_id   <- rownames(tpm_mat)
tpm_mat$gene_name <- gene_map$gene_name[match(tpm_mat$gene_id, gene_map$gene_id)]
tpm_sym <- tpm_mat %>%
  mutate(gene = coalesce(gene_name, gene_id)) %>%
  select(gene, all_of(colnames(vst_mat)))

tpm_median <- tpm_sym %>%
  mutate(medTPM = matrixStats::rowMedians(as.matrix(select(., -gene)), na.rm = TRUE)) %>%
  select(gene, medTPM)
vst_var <- tibble(
  gene  = rownames(vst_mat),
  varVST = matrixStats::rowVars(as.matrix(vst_mat), na.rm = TRUE)
)

add_names <- function(res) {
  tb <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")
  gn <- gene_map$gene_name[match(tb$gene_id, gene_map$gene_id)]
  tb$gene <- coalesce(gn, tb$gene_id)
  tb %>% select(gene, log2FoldChange, padj)
}

de_all <- dplyr::bind_rows(
  add_names(res_oldCut_vs_youngCut)  %>% mutate(comp = "OldCut_vs_YoungCut"),
  add_names(res_oldCut_vs_oldUI)     %>% mutate(comp = "OldCut_vs_OldUI"),
  add_names(res_youngCut_vs_youngUI) %>% mutate(comp = "YoungCut_vs_YoungUI")
)
de_summary <- de_all %>%
  group_by(gene) %>%
  summarise(max_absLFC = max(abs(log2FoldChange), na.rm = TRUE),
            min_padj    = suppressWarnings(min(padj, na.rm = TRUE)),
            .groups = "drop")
expr_var <- tpm_median %>% inner_join(vst_var, by = "gene") %>% inner_join(de_summary, by = "gene")

q25_var <- quantile(expr_var$varVST, 0.25, na.rm = TRUE)
hk_auto <- expr_var %>%
  filter(medTPM >= 10, varVST <= q25_var, min_padj >= 0.5, max_absLFC <= 0.2) %>%
  arrange(varVST, desc(medTPM)) %>%
  slice_head(n = 8) %>%
  pull(gene)
if (length(hk_auto) == 0) {
  hk_auto <- expr_var %>%
    filter(medTPM >= 5, varVST <= quantile(varVST, 0.35, na.rm = TRUE),
           min_padj >= 0.3, max_absLFC <= 0.3) %>%
    arrange(varVST, desc(medTPM)) %>%
    slice_head(n = 8) %>%
    pull(gene)
}
hk_genes <- unique(hk_auto)
message("HK auto seleccionados: ", paste(hk_genes, collapse = ", "))

## 13) Heatmaps de expresión (matriz VST real y matriz de log2FC entre grupos, con indicación de HK)
anno <- as.data.frame(coldata[, c("Age","Injury")])
ann_colors <- list(Age = c(Young = "#4DAF4A", Old = "#E41A1C"),
                   Injury = c(UI = "#377EB8", Cut3dpi = "#FF7F00"))

plot_module_heatmap_vst <- function(genes, module_name, mat, anno, cap_quant = c(0.02, 0.98)) {
  rn_norm <- normalize_names(rownames(mat))
  want <- normalize_names(genes)
  sel <- rownames(mat)[rn_norm %in% want]
  if (length(sel) < 2) {
    message("Módulo ", module_name, ": pocos genes presentes")
    return(invisible(NULL))
  }
  M <- mat[sel, , drop = FALSE]
  is_hk <- normalize_names(rownames(M)) %in% normalize_names(hk_genes)
  row_anno <- data.frame(Category = ifelse(is_hk, "Housekeeping", "Module"), row.names = rownames(M))
  ord <- order(!is_hk, rownames(M))
  M <- M[ord, , drop = FALSE]; row_anno <- row_anno[rownames(M), , drop = FALSE]
  gaps_row <- if (any(is_hk)) sum(is_hk[ord]) else NULL
  breaks <- if (!is.null(cap_quant)) {
    lims <- quantile(as.numeric(M), cap_quant, na.rm = TRUE)
    seq(lims[1], lims[2], length.out = 101)
  } else seq(min(M), max(M), length.out = 101)
  png_h <- max(8, 0.20 * nrow(M))
  pheatmap::pheatmap(M, scale = "none",
           annotation_col = anno, annotation_row = row_anno,
           annotation_colors = c(ann_colors, list(Category = c(Housekeeping = "grey30", Module = "white"))),
           cluster_cols = FALSE, cluster_rows = TRUE, gaps_row = gaps_row,
           show_rownames = TRUE, fontsize_row = 7,
           main = paste0("Heatmap ", module_name, "  VST (orden fijo)"),
           color = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101),
           breaks = breaks,
           filename = file.path(fig_dir, paste0("heatmap_VST_", gsub("[ /]", "_", module_name), ".png")),
           width = 8, height = png_h)
}

# Función auxiliar para obtener p-valor ajustado de resultados en formato tabla
.add_names <- function(res, comp_label) {
  tb <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")
  gn <- gene_map$gene_name[match(tb$gene_id, gene_map$gene_id)]
  tb$gene <- dplyr::coalesce(gn, tb$gene_id)
  tb %>% select(gene, padj) %>% rename(!!comp_label := padj)
}

# Versión robusta de heatmap group log2FC (resumen por grupos con estrellas significativas)
heatmap_group_log2fc <- function(genes, module_name,
                                 ref_group = "Young_UI",
                                 order_groups = c("Young_Cut3dpi","Old_UI","Old_Cut3dpi"),
                                 range = c(-2, 2), show_stars = TRUE) {
  sel <- intersect(rownames(vst_mat), genes)
  if (length(sel) < 2L) {
    message("Heatmap log2FC: pocos genes (", module_name, ")")
    return(invisible(NULL))
  }
  # Calcular matriz de log2FC por grupo (cada gen centrado respecto al grupo ref)
  vst_means <- as.data.frame(vst_mat[sel, , drop = FALSE]) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "vst") %>%
    left_join(as_tibble(coldata, rownames = "sample"), by = "sample") %>%
    group_by(gene, Group) %>%
    summarise(vst_mean = mean(vst), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Group, values_from = vst_mean)
  stopifnot(ref_group %in% colnames(vst_means))
  for (g in setdiff(colnames(vst_means), c("gene", ref_group))) {
    vst_means[[g]] <- vst_means[[g]] - vst_means[[ref_group]]
  }
  mat_fc <- vst_means %>%
    select(gene, all_of(order_groups)) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  # Añadir estrellas indicando significancia estadística
  numbers <- NULL
  if (show_stars) {
    min_or_na <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else min(x) }
    padj_tab <- dplyr::bind_rows(
      .add_names(res_youngCut_vs_youngUI, "Young_Cut3dpi"),
      .add_names(res_oldCut_vs_oldUI,     "Old_UI"),
      .add_names(res_oldCut_vs_youngCut,  "Old_Cut3dpi")
    ) %>%
      group_by(gene) %>%
      summarise(
        Young_Cut3dpi = min_or_na(Young_Cut3dpi),
        Old_UI        = min_or_na(Old_UI),
        Old_Cut3dpi   = min_or_na(Old_Cut3dpi),
        .groups = "drop"
      ) %>%
      tibble::column_to_rownames("gene") %>%
      as.matrix()
    # Alinear y convertir p-vals a símbolos de estrellas
    common_rows <- intersect(rownames(mat_fc), rownames(padj_tab))
    mat_fc   <- mat_fc[common_rows, , drop = FALSE]
    padj_tab <- padj_tab[common_rows, colnames(mat_fc), drop = FALSE]
    star_fun <- function(p) ifelse(is.na(p), "", ifelse(p <= 0.001, "***",
                                             ifelse(p <= 0.01, "**",
                                             ifelse(p <= 0.05, "*", ""))))
    numbers <- apply(padj_tab, c(1, 2), star_fun)
  }
  # Ordenar filas (HK al final) y preparar anotaciones
  is_hk <- normalize_names(rownames(mat_fc)) %in% normalize_names(hk_genes)
  row_anno <- data.frame(Category = ifelse(is_hk, "Housekeeping", "Module"), row.names = rownames(mat_fc))
  ord <- order(!is_hk, rownames(mat_fc))
  mat_fc <- mat_fc[ord, , drop = FALSE]
  row_anno <- row_anno[rownames(mat_fc), , drop = FALSE]
  if (!is.null(numbers)) numbers <- numbers[rownames(mat_fc), , drop = FALSE]
  gaps_row <- if (any(is_hk)) sum(is_hk[ord]) else NULL

  png_h <- max(8, 0.20 * nrow(mat_fc))
  pheatmap::pheatmap(
    mat_fc, scale = "none",
    cluster_cols = FALSE, cluster_rows = TRUE, gaps_row = gaps_row,
    annotation_row = row_anno,
    annotation_colors = list(Category = c(Housekeeping = "grey30", Module = "white")),
    fontsize_row = 7,
    main = paste0("Resumen por grupo  ", module_name, " (log2FC vs ", ref_group, ")"),
    color = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101),
    breaks = seq(range[1], range[2], length.out = 101),
    display_numbers = numbers,
    number_color = "black",
    fontsize_number = 9,
    filename = file.path(fig_dir, paste0("heatmap_log2FC_", gsub("[ /]", "_", module_name), "_vs_", ref_group, ".png")),
    width = 6.5, height = png_h
  )
}

## 14) Calcular scores de expresión por módulo (utilizando método de ranking) y graficar
rank_based_score <- function(expr_mat, genes_up) {
  rk <- apply(expr_mat, 2, rank, ties.method = "average")
  g_up <- intersect(rownames(expr_mat), genes_up)
  if (length(g_up) < 3) warning("Conjunto pequeño presente: ", length(g_up), " genes.")
  avg_rk <- colMeans(rk[g_up, , drop = FALSE])
  G <- nrow(expr_mat)
  score <- (avg_rk - 1) / (G - 1)
  return(score)
}
compute_module_scores <- function(vst_mat, modules_named_list) {
  scores <- lapply(names(modules_named_list), function(nm) {
    g <- unique(modules_named_list[[nm]])
    s <- rank_based_score(vst_mat, g)
    setNames(as.numeric(s), colnames(vst_mat))
  })
  names(scores) <- names(modules_named_list)
  as.data.frame(scores) %>% tibble::rownames_to_column("sample")
}
plot_score_box <- function(df, title, file_out) {
  g <- ggplot(df, aes(x = Group, y = score, fill = Age)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.8) +
    scale_x_discrete(limits = c("Young_UI","Young_Cut3dpi","Old_UI","Old_Cut3dpi")) +
    labs(title = title, y = "Score (0–1)", x = "Grupo") + theme_bw(base_size = 12)
  ggsave(file_out, g, width = 7.5, height = 5, dpi = 300)
}

## 15) Ejecutar: generar heatmaps y boxplots de scores para cada super-módulo y submódulo
# Heatmaps para todos los genes de cada super-módulo
purrr::iwalk(MODULES, function(mod_list, super_name) {
  genes_all <- unique(unlist(mod_list))
  plot_module_heatmap_vst(genes_all, paste0(super_name, "_All"), vst_mat, anno)
  heatmap_group_log2fc(genes_all, paste0(super_name, "_All"))
})
# Heatmaps para cada sub-módulo individual
purrr::iwalk(MODULES, function(mod_list, super_name) {
  purrr::iwalk(mod_list, function(gset, sub_name) {
    plot_module_heatmap_vst(gset, paste0(super_name, "__", sub_name), vst_mat, anno)
    heatmap_group_log2fc(gset, paste0(super_name, "__", sub_name))
  })
})

# Scores por sub-módulo (cálculo y boxplots)
scores_long <- list()
for (super_name in names(MODULES)) {
  mod_list <- MODULES[[super_name]]
  sc <- compute_module_scores(vst_mat, mod_list)
  sc_long <- sc %>%
    tidyr::pivot_longer(-sample, names_to = "Module", values_to = "score") %>%
    left_join(as_tibble(coldata, rownames = "sample"), by = "sample") %>%
    mutate(Group = factor(Group, levels = c("Young_UI","Young_Cut3dpi","Old_UI","Old_Cut3dpi")))
  scores_long[[super_name]] <- sc_long
  plot_score_box(sc_long, paste0("Scores  ", super_name),
                 file.path(fig_dir, paste0("Scores_", gsub("[ /]", "_", super_name), ".png")))
}
scores_all <- dplyr::bind_rows(scores_long, .id = "Question")

## 16) Exportar resultados a tablas (incluyendo Excel con múltiples hojas)
add_gene_names <- function(res) {
  tb <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id") %>% as_tibble()
  tb$gene_name <- gene_map$gene_name[match(tb$gene_id, gene_map$gene_id)]
  tb %>% relocate(gene_id, gene_name) %>% arrange(padj, desc(abs(log2FoldChange)))
}
tab1 <- add_gene_names(res_oldCut_vs_youngCut)  %>% mutate(comparison = "Old_Cut3dpi_vs_Young_Cut3dpi")
tab2 <- add_gene_names(res_oldCut_vs_oldUI)     %>% mutate(comparison = "Old_Cut3dpi_vs_Old_UI")
tab3 <- add_gene_names(res_youngCut_vs_youngUI) %>% mutate(comparison = "Young_Cut3dpi_vs_Young_UI")

# Scores por muestra (formato ancho)
scores_wide <- scores_all %>% select(sample, Module, score) %>%
  tidyr::pivot_wider(names_from = Module, values_from = score) %>%
  left_join(as_tibble(coldata, rownames = "sample"), by = "sample") %>%
  relocate(sample, Age, Injury, Group)

# Scores por grupo (promedios por grupo experimental)
scores_group <- scores_all %>%
  group_by(Question, Module, Group) %>% 
  summarise(score_mean = mean(score), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Group, values_from = score_mean)

# Resumen de expresión (VST) de genes Treg por submódulo (promedio por grupo)
get_group_means <- function(genes) {
  sel <- intersect(rownames(vst_mat), genes)
  if (!length(sel)) return(NULL)
  vst_tbl <- as_tibble(vst_mat[sel, , drop = FALSE], rownames = "gene")
  vst_long <- vst_tbl %>% tidyr::pivot_longer(-gene, names_to = "sample", values_to = "vst") %>%
    left_join(as_tibble(coldata, rownames = "sample"), by = "sample")
  vst_long %>% group_by(gene, Group) %>% summarise(vst = mean(vst), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Group, values_from = vst)
}
Treg_summaries <- purrr::imap(Q3, ~ { gm <- get_group_means(.x); if (!is.null(gm)) mutate(gm, Submodule = .y) }) %>% bind_rows()

# Crear libro Excel con todos los resultados relevantes
wb <- createWorkbook()
addWorksheet(wb, "DE_OldCut_vs_YoungCut");  writeData(wb, 1, tab1)
addWorksheet(wb, "DE_OldCut_vs_OldUI");     writeData(wb, 2, tab2)
addWorksheet(wb, "DE_YoungCut_vs_YoungUI"); writeData(wb, 3, tab3)
addWorksheet(wb, "ModuleScores_perSample"); writeData(wb, 4, scores_wide)
addWorksheet(wb, "ModuleScores_GroupMeans"); writeData(wb, 5, scores_group)
addWorksheet(wb, "TregMarkers_Summary");     writeData(wb, 6, Treg_summaries)
# (Agregar hojas extra de control, si se desea:)
vst_tbl <- as_tibble(vst_mat, rownames = "gene")
addWorksheet(wb, "VST_matrix");              writeData(wb, 7, vst_tbl)
addWorksheet(wb, "HK_selected");             writeData(wb, 8, tibble(HK = hk_genes))

saveWorkbook(wb, file.path(res_dir, "INMUNO_TREG_SENESCENCIA_results.xlsx"), overwrite = TRUE)

message(" Fin. Figuras: ", fig_dir, " | Excel: ", file.path(res_dir, "INMUNO_TREG_SENESCENCIA_results.xlsx"))
