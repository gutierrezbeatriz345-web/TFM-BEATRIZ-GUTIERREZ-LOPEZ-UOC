# Código R - Enriquecimiento GO/KEGG por bloques (módulos Q1–Q5) en distintos contrastes

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in c("clusterProfiler","org.Mm.eg.db","enrichplot")) {
    if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
  }
  library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
  library(dplyr); library(purrr); library(readr); library(ggplot2)
})

# Helper: obtener genes DE significativos para un bloque (lista de genes) en un contraste dado
# Usa resultados DESeq2 ya calculados: res_youngCut_vs_youngUI, res_oldCut_vs_oldUI, res_oldCut_vs_youngCut
get_de_genes_for_block <- function(res_obj, block_genes, padj_thr = 0.05, lfc_thr = 0.5) {
  tb <- as.data.frame(res_obj) %>% tibble::rownames_to_column("gene_id")
  # Mapear gene_id -> símbolo (usando gene_map existente de bulk_RNAseq)
  gn <- gene_map$gene_name[match(tb$gene_id, gene_map$gene_id)]
  tb$gene <- coalesce(gn, tb$gene_id)
  tb <- tb %>%
    filter(!is.na(padj), padj <= padj_thr, abs(log2FoldChange) >= lfc_thr) %>%
    filter(toupper(gsub("[^A-Za-z0-9]", "", gene)) %in% toupper(gsub("[^A-Za-z0-9]", "", block_genes))) %>%
    pull(gene) %>% unique()
  return(tb)
}

# Helper: ejecutar enriquecimiento GO BP y KEGG dado un vector de símbolos génicos
do_enrichments <- function(genes_symbol, universe_symbol = NULL, out_prefix = "Qx_BlockX") {
  if (length(genes_symbol) < 5) return(NULL)
  # Mapear símbolos a IDs de Entrez
  eg <- bitr(genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  if (!is.null(universe_symbol)) {
    bg <- bitr(universe_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    bg_ids <- unique(bg$ENTREZID)
  } else {
    bg_ids <- NULL
  }
  if (nrow(eg) < 5) return(NULL)
  gene_ids <- unique(eg$ENTREZID)
  # Enriquecimiento GO Biological Process (BP)
  ego <- enrichGO(gene         = gene_ids,
                  universe     = bg_ids,
                  OrgDb        = org.Mm.eg.db,
                  keyType      = "ENTREZID",
                  ont          = "BP",
                  pAdjustMethod= "BH",
                  qvalueCutoff = 0.10,
                  readable     = TRUE)
  # Enriquecimiento KEGG
  ekk <- enrichKEGG(gene         = gene_ids,
                    organism     = "mmu",
                    pAdjustMethod= "BH",
                    qvalueCutoff = 0.10)
  # Guardar resultados y graficar top términos
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write_csv(as.data.frame(ego), file.path(res_dir, paste0(out_prefix, "_GO_BP.csv")))
    p1 <- barplot(ego, showCategory = 15)
    ggsave(file.path(fig_dir, paste0(out_prefix, "_GO_BP_bar.png")), p1, width = 8, height = 6, dpi = 300)
  }
  if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0) {
    write_csv(as.data.frame(ekk), file.path(res_dir, paste0(out_prefix, "_KEGG.csv")))
    p2 <- barplot(ekk, showCategory = 15)
    ggsave(file.path(fig_dir, paste0(out_prefix, "_KEGG_bar.png")), p2, width = 8, height = 6, dpi = 300)
  }
  return(invisible(list(GO = ego, KEGG = ekk)))
}

# Preparar bloques (super-módulos) y lista de contrastes
BLOCKS <- lapply(MODULES, function(sublist) unique(unlist(sublist)))
contrastes <- list(
  YoungCut_vs_YoungUI = res_youngCut_vs_youngUI,
  OldCut_vs_OldUI     = res_oldCut_vs_oldUI,
  OldCut_vs_YoungCut  = res_oldCut_vs_youngCut
)

# Recorrer cada bloque y cada contraste para realizar enriquecimiento
for (blk in names(BLOCKS)) {
  block_genes <- BLOCKS[[blk]]
  # Universo de fondo opcional: todos los genes detectados (símbolos)
  universe_symbol <- rownames(vst_mat)
  for (cname in names(contrastes)) {
    de_syms <- get_de_genes_for_block(contrastes[[cname]], block_genes,
                                      padj_thr = 0.05, lfc_thr = 0.5)
    if (length(de_syms) >= 5) {
      do_enrichments(de_syms, universe_symbol, out_prefix = paste0(blk, "__", cname))
    }
  }
}
