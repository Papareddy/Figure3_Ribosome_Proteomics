library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(biomaRt)

# ---------------------------------------------------------
# 1. Get Ensembl gene info
# ---------------------------------------------------------
get_gene_info <- function() {
  mart <- useMart(
    "plants_mart",
    dataset = "athaliana_eg_gene",
    host = "https://plants.ensembl.org"
  )
  
  getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    mart = mart
  )
}

gene_info <- get_gene_info()

# ---------------------------------------------------------
# 2. Manual naming for tricomplex genes
# ---------------------------------------------------------
manual_map <- c(
  "AT4G27120" = "DDRGK1",
  "AT5G06830" = "C53",
  "AT3G46220" = "UFL1"
)

# ---------------------------------------------------------
# 3. Process each file
# ---------------------------------------------------------
process_file <- function(file_path) {
  
  data <- readxl::read_excel(file_path) %>%
    setNames(tolower(names(.)))
  
  req <- c(
    "gene_id", "peptide_count",
    "ans_vs_dmso_wt.fc",
    "ans_vs_dmso_ufm1.fc",
    "ans_vs_dmso_wt.pval",
    "ans_vs_dmso_ufm1.pval"
  )
  
  if (!all(req %in% names(data))) {
    missing <- setdiff(req, names(data))
    stop(paste0("Missing required columns in ", basename(file_path), ": ",
                paste(missing, collapse = ", ")))
  }
  
  df <- data %>%
    dplyr::select(
      gene_id,
      peptide_count,
      fc_wt     = `ans_vs_dmso_wt.fc`,
      fc_ufm1   = `ans_vs_dmso_ufm1.fc`,
      pval_wt   = `ans_vs_dmso_wt.pval`,
      pval_ufm1 = `ans_vs_dmso_ufm1.pval`
    ) %>%
    dplyr::filter(peptide_count > 0) %>%
    dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
    dplyr::left_join(gene_info, by = c("gene_id" = "ensembl_gene_id")) %>%
    dplyr::mutate(
      external_gene_name = ifelse(
        gene_id %in% names(manual_map),
        manual_map[gene_id],
        external_gene_name
      ),
      fc_wt   = pmax(pmin(fc_wt, 4), -4),
      fc_ufm1 = pmax(pmin(fc_ufm1, 4), -4)
    )
  
  prefix <- tools::file_path_sans_ext(basename(file_path))
  
  df_final <- df %>%
    dplyr::mutate(file = prefix) %>%
    dplyr::select(
      file,
      gene_id,
      external_gene_name,
      description,
      log2fc_wt   = fc_wt,
      log2fc_ufm1 = fc_ufm1,
      pval_wt,
      pval_ufm1
    )
  
  return(df_final)
}

# ---------------------------------------------------------
# 4. Run on all files → ONE dataframe
# ---------------------------------------------------------
excel_dir <- "~/Desktop/Export_MS/Excels/header_change/Excels_backup/"

excel_files <- list.files(
  excel_dir, pattern = "\\.xlsx$", full.names = TRUE
)

master_df <- bind_rows(map(excel_files, process_file))

# ---------------------------------------------------------
# 5. Ordering of fractions/timepoints
# ---------------------------------------------------------
desired_order <- c(
  "4hr_60S",
  "Bulked_60S",
  "16hr_60S",
  "4hr_Monosome",
  "Bulked_Monosome",
  "16hr_Monosome",
  "4hr_Disome",
  "Bulked_Disome",
  "16hr_Disome"
)

master_df$file <- factor(master_df$file, levels = desired_order)

master_df <- master_df %>% arrange(file)

AT5G61140
master_df %>%
  filter(gene_id == "AT4G27650") %>%
  arrange(file)

# ---------------------------------------------------------
# 6. BASE R BARPLOT FOR ANY GENE (WT + UFM1)
# ---------------------------------------------------------
par(mfrow=c(3,4),las=2,tcl=-0.3,bty="o")
gene_name <- "MBF1B"   # <-- set the gene you want

df <- master_df %>%filter(external_gene_name == gene_name) %>%arrange(file)
df

if (nrow(df) == 0) {stop("Gene not found.")}

# Prepare grouped matrix WT vs UFM1
mat <- rbind(WT   = df$log2fc_wt,UFM1 = df$log2fc_ufm1)
mat[is.na(mat)] <- 0
#mat <- pmax(pmin(mat, 1), -1)
# Colors
cols <- c("WT" = "#505193", "UFM1" = "#03A791")

# Barplot (base R)
barplot(
  mat,
  beside = TRUE,
  names.arg = df$file,     # <-- FIXED (actual file names)
  col = cols,
  border=cols,
  las = 2,
  horiz = F,
  ylab = "log2 Fold Change",
  main = paste(gene_name, "Abundance")
)
box()
#abline(h = 0, lwd = 0.5, col="grey20", lty=2)
