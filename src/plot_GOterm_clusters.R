library(ggplot2)
library(tidyr)
library(forcats)
library(ggpubr)
library(dplyr)
library(readr)

# --- CATCH STUN OPTIONS FROM TERMINAL ---
args <- commandArgs(trailingOnly = TRUE)

results_dir  <- if(length(args) >= 1) args[1] else "GO_Results_UP"
output_file  <- if(length(args) >= 2) args[2] else "Figure3_GO_UP_Dotplot.pdf"
plot_title   <- if(length(args) >= 3) args[3] else "Ribosome Associated Proteome GO Terms"

# These match the order passed by yolo_run.sh
opt_top_n    <- if(length(args) >= 4) as.numeric(args[4]) else 4
opt_term_max <- if(length(args) >= 5) as.numeric(args[5]) else 300
opt_pval     <- if(length(args) >= 6) as.numeric(args[6]) else 0.01
opt_fold     <- if(length(args) >= 7) as.numeric(args[7]) else 5

message(paste("Processing directory:", results_dir))
message(paste("Stun Filters: TopN =", opt_top_n, "| MaxTerm =", opt_term_max, "| P <", opt_pval, "| Fold >", opt_fold))

# Unified GO data loader utilizing the opt_ parameters
load_go_file <- function(file, ontology) {
  if (!file.exists(file)) return(NULL)
  
  df <- read_tsv(file, show_col_types = FALSE)
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>%
    mutate(
      Fold_Enrichment = round((intersection_size / query_size) / (term_size / effective_domain_size), 2),
      ontology = ontology
    ) %>%
    filter(
      term_size <= opt_term_max,
      p_value <= opt_pval,
      Fold_Enrichment >= opt_fold
    ) %>%
    group_by(Cluster) %>%
    arrange(p_value) %>%
    slice_head(n = opt_top_n) %>%
    ungroup() %>%
    mutate(
      short_name = paste0(native, ": ", sapply(strsplit(name, " "), function(words) paste(head(words, 6), collapse = " ")))
    )
  
  return(df)
}

# Load data
bp_df <- load_go_file(file.path(results_dir, "go_terms_bp.tsv"), "Biological Process")
mf_df <- load_go_file(file.path(results_dir, "go_terms_mf.tsv"), "Molecular Function")
cc_df <- load_go_file(file.path(results_dir, "go_terms_cc.tsv"), "Cellular Component")

all_go_df <- bind_rows(bp_df, mf_df, cc_df)
if (is.null(all_go_df) || nrow(all_go_df) == 0) stop("No GO terms passed filters.")

# --- DYNAMIC CLUSTER ORDERING ---
fractions <- c("4hr_60S", "4hr_Monosome", "4hr_Disome", 
               "Bulked_60S", "Bulked_Monosome", "Bulked_Disome", 
               "16hr_60S", "16hr_Monosome", "16hr_Disome")

is_up <- grepl("UP", results_dir, ignore.case = TRUE)
suffix <- if(is_up) "_up" else "_down"
current_order <- paste0(fractions, suffix)
available_clusters <- intersect(current_order, unique(all_go_df$Cluster))

# Define plotting function
plot_go_terms <- function(df, title_text) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df_plot <- df %>%
    filter(Cluster %in% available_clusters) %>%
    mutate(
      Cluster = factor(Cluster, levels = available_clusters),
      short_name = factor(short_name, levels = rev(unique(short_name))),
      ontology = factor(ontology, levels = c("Biological Process", "Molecular Function", "Cellular Component"))
    )

  ggplot(df_plot, aes(x = Cluster, y = short_name, fill = Fold_Enrichment, size = -log10(p_value))) +
    geom_point(shape = 21, color = "black") +
    scale_fill_viridis_c(option = "D", name = "Fold Enrichment") +
    scale_size_continuous(name = "-log10(FDR)") +
    facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
    labs(title = title_text, x = "Ribosome Fraction", y = "GO Term") +
    theme_minimal(9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 8, colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      strip.text.y = element_text(size = 9, face = "bold")
    )
}

# --- SAVE MULTI-PAGE PDF ---
pdf(output_file, width = 10, height = 10)

# Page 1: BP
p_bp <- plot_go_terms(bp_df, paste(plot_title, "- Biological Process"))
if (!is.null(p_bp)) print(p_bp)

# Page 2: MF
p_mf <- plot_go_terms(mf_df, paste(plot_title, "- Molecular Function"))
if (!is.null(p_mf)) print(p_mf)

# Page 3: CC
p_cc <- plot_go_terms(cc_df, paste(plot_title, "- Cellular Component"))
if (!is.null(p_cc)) print(p_cc)

# Page 4: Combined
p_all <- plot_go_terms(all_go_df, paste(plot_title, "- All Categories"))
if (!is.null(p_all)) print(p_all)

dev.off()
message(">> Success: Multi-page GO Dotplots generated.")