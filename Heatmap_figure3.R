library(dplyr)
library(tidyr)
library(pheatmap)
library(gridExtra)
library(grid)

# -----------------------------
# Order of the 9 panels
# -----------------------------
files9 <- c(
  "4hr_60S","4hr_Monosome","4hr_Disome",
  "Bulked_60S","Bulked_Monosome","Bulked_Disome",
  "16hr_60S","16hr_Monosome","16hr_Disome"
)

# -----------------------------
# Sanitize matrix (NA/Inf safe)
# -----------------------------
sanitize_mat <- function(m) {
  m <- apply(m, 2, as.numeric)
  m[!is.finite(m)] <- NA
  m[is.na(m)] <- 0
  m
}

# -----------------------------
# Build ONE heatmap grob
# -----------------------------
make_heatmap_grob <- function(f) {
  
  df <- master_df_sig %>% filter(file == f)
  
  sig_genes_f <- df %>%
    filter(significant_in_wt %in% c("up","down")) %>%
    distinct(gene_id)
  
  if (nrow(sig_genes_f) == 0) {
    return(grid::nullGrob())
  }
  
  hm <- df %>%
    semi_join(sig_genes_f, by = "gene_id") %>%
    select(gene_id, log2fc_wt, log2fc_ufm1) %>%
    group_by(gene_id) %>%
    summarise(
      WT   = mean(log2fc_wt, na.rm = TRUE),
      UFM1 = mean(log2fc_ufm1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()
  
  rownames(hm) <- hm$gene_id
  hm$gene_id <- NULL
  
  hm_mat <- sanitize_mat(as.matrix(hm))
  hm_breaks <- seq(-3, 3, length.out = 101)
  
  # return grob, do NOT plot yet
  pheatmap(
    hm_mat,gaps_col=1,
    cutree_rows=3,
    border_color=NA,
    cellwidth = 15,
    cluster_rows = T,
    cluster_cols = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    show_rownames = FALSE,
    silent = TRUE,
    main = gsub("_", " ", f),
    breaks = hm_breaks,
    color = colorRampPalette(c("#29558b", "white", "#d86d83"))(100)
  )$gtable
}

# -----------------------------
# Build all 9 heatmaps
# -----------------------------
heatmap_grobs <- lapply(files9, make_heatmap_grob)

# -----------------------------
# Draw as ONE 3×3 panel
# -----------------------------
grid.arrange(
  grobs = heatmap_grobs,
  ncol = 3
)
