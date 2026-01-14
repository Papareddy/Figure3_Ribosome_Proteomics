library(dplyr)

# ---------- Setup & Helpers ----------
adjust_col <- function(hex, alpha = 0.75) adjustcolor(hex, alpha.f = alpha)

# Hyperbolic threshold function
threshold_func <- function(x_vals, c_val, x0_val) {
  ifelse(abs(x_vals) > x0_val, c_val / (abs(x_vals) - x0_val), NA_real_)
}

# Threshold parameters
c_val <- 1.5
x0_1 <- 0
x0_2 <- 1

# Desired plot order
fractions <- c(
  "4hr_60S", "4hr_Monosome", "4hr_Disome",
  "Bulked_60S", "Bulked_Monosome", "Bulked_Disome",
  "16hr_60S", "16hr_Monosome", "16hr_Disome"
)

# Genes to highlight (Aesthetic focus)
highlight_genes <- c("AT4G27120", "AT5G06830", "AT3G46220")
highlight_colors <- c(
  "AT4G27120" = "#129990",
  "AT5G06830" = "#27548A",
  "AT3G46220" = "#DDA853"
)

# ---------- 1) Load and Prep Data ----------
master_df <- read.delim("data/Ribosome_Associated_proteome.tsv", header=TRUE, sep="\t")

# Standardize column names if necessary (ensure they match your header)
# Using 'file', 'log2fc_wt', 'pval_wt', 'gene_id' from your snippet
master_df_sig <- master_df %>%
  group_by(file) %>%
  mutate(
    neglog10p_raw = -log10(pval_wt),
    thr_1 = threshold_func(log2fc_wt, c_val, x0_1),
    thr_2 = threshold_func(log2fc_wt, c_val, x0_2),
    sig_mask_wt = pval_wt <= 0.05 & (
      (abs(log2fc_wt) > x0_2 & neglog10p_raw > thr_2) |
      (abs(log2fc_wt) > x0_1 & neglog10p_raw > thr_1)
    ),
    significant_in_wt = case_when(
      sig_mask_wt & log2fc_wt > 0 ~ "up",
      sig_mask_wt & log2fc_wt < 0 ~ "down",
      TRUE ~ "NS"
    )
  ) %>%
  ungroup()

# ---------- 2) Plot WT Volcanoes (Base R Aesthetic) ----------
pdf("results/Figure3_Proteome_Volcanoes.pdf", width = 12, height = 12)
par(mfrow = c(3, 3), mar = c(4.5, 4.5, 3, 1), las = 1, tcl = -0.3, bty = "l")

for (fraction in fractions) {
  
  rb <- master_df_sig %>% filter(file == fraction)
  
  # Skip empty panels cleanly
  if (nrow(rb) == 0) {
    plot.new()
    title(main = gsub("_", " ", fraction))
    next
  }
  
  log2fc_raw <- rb$log2fc_wt
  neglog10p_raw <- -log10(rb$pval_wt)
  
  # Capping for visual consistency
  log2fc <- pmin(pmax(log2fc_raw, -4), 4)
  neglog10p <- pmin(neglog10p_raw, 8)
  
  sig_mask <- rb$significant_in_wt != "NS"
  
  # Aesthetic Colors
  point_color <- ifelse(sig_mask,
                        adjust_col("#d86d83", 0.75), # Muted Rose for Sig
                        adjust_col("#d3d3d8", 0.50)) # Light Grey for NS
  
  point_cex <- rep(0.8, length(log2fc))
  
  # Highlight specific genes
  for (gene in highlight_genes) {
    idx <- which(rb$gene_id == gene)
    if (length(idx) > 0) {
      point_color[idx] <- adjust_col(highlight_colors[gene], 1.0)
      point_cex[idx] <- 1.8
    }
  }
  
  # Base Plot
  plot(
    log2fc, neglog10p,
    pch = 16,
    col = point_color,
    cex = point_cex,
    xlab = expression(log[2]~Fold~Change~(WT)),
    ylab = expression(-log[10]~P-value~(WT)),
    main = gsub("_", " ", fraction),
    xlim = c(-4, 4),
    ylim = c(0, 8),
    axes = FALSE
  )
  
  # Custom Axes for "l" box type aesthetic
  axis(1, at = seq(-4, 4, 2))
  axis(2, at = seq(0, 8, 2))
  
  # Reference Lines
  abline(h = -log10(0.05), col = "grey70", lty = 3)
  
  # Effect-size threshold curve (Curved boundary)
  # Plotting the curve for both positive and negative FC
  curve(threshold_func(x, c_val, x0_1), from = 0.1, to = 4, col = "grey50", lty = 2, add = TRUE, lwd = 1)
  curve(threshold_func(x, c_val, x0_1), from = -4, to = -0.1, col = "grey50", lty = 2, add = TRUE, lwd = 1)
  
  # Stats Annotation
  n_up <- sum(rb$significant_in_wt == "up", na.rm = TRUE)
  n_down <- sum(rb$significant_in_wt == "down", na.rm = TRUE)
  
  legend(
    "topleft",
    legend = c(paste0("Up: ", n_up), paste0("Down: ", n_down)),
    bty = "n", text.col = "black", cex = 0.85, inset = c(0, 0)
  )
}

dev.off()
message(">> Success: Aesthetic Curved Volcanoes generated.")