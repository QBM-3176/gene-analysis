# ============================================================
# Gene Expression Boxplots from GEO Series Matrix
# Dataset: GSE104954-GPL22945
# Genes: 3026_at, 3426_at, 3075_at, 1305_at, 725_at
# ============================================================

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

# ---- 1. Load the file ----
file_path <- "GSE104954-GPL22945_series_matrix.txt"   # adjust path as needed

lines <- readLines(file_path)

# ---- 2. Extract sample titles and diagnosis labels ----

# Sample titles (first !Sample_title line)
title_line <- grep("^!Sample_title", lines, value = TRUE)[1]
sample_titles <- gsub('"', '', strsplit(title_line, "\t")[[1]][-1])

# Diagnosis: second !Sample_characteristics_ch1 line contains diagnosis
char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
diag_line  <- char_lines[2]                           # second occurrence = diagnosis
diag_vals  <- gsub('"', '', strsplit(diag_line, "\t")[[1]][-1])

# Clean up diagnosis labels
diag_vals <- sub("^diagnosis: ", "", diag_vals)
diag_vals[diag_vals == ""] <- "Living Donor"          # empty = living donor controls

# Abbreviate long labels for readability on plots
diag_vals <- gsub("Diabetic nephropathy",                                         "DN",       diag_vals)
diag_vals <- gsub("Focal segmental glomerulosclerosis and minimal change disease", "FSGS+MCD", diag_vals)
diag_vals <- gsub("Focal segmental glomerulosclerosis",                           "FSGS",     diag_vals)
diag_vals <- gsub("Minimal change disease",                                       "MCD",      diag_vals)
diag_vals <- gsub("ANCA-associated vasculitis",                                   "RPGN",     diag_vals)
diag_vals <- gsub("Tumor nephrectomy",                                            "TN",       diag_vals)
diag_vals <- gsub("Living Donor",                                                 "LD",       diag_vals)

# ---- 3. Parse the expression matrix ----
table_start <- grep("^!series_matrix_table_begin", lines) + 1
table_end   <- grep("^!series_matrix_table_end",   lines) - 1

matrix_lines <- lines[table_start:table_end]
con  <- textConnection(matrix_lines)
expr <- read.table(con, header = TRUE, sep = "\t",
                   row.names = 1, check.names = FALSE,
                   quote = '"')
close(con)

# ---- 4. Select the five genes of interest ----
genes_of_interest <- c("3026_at", "3426_at", "3075_at", "1305_at", "725_at")
missing <- setdiff(genes_of_interest, rownames(expr))
if (length(missing) > 0) stop("Gene(s) not found: ", paste(missing, collapse = ", "))

expr_sub <- expr[genes_of_interest, , drop = FALSE]

# ---- 5. Build a tidy data frame ----
expr_df <- as.data.frame(t(expr_sub))
expr_df$Sample    <- sample_titles
expr_df$Diagnosis <- diag_vals

# Ordered factor so boxes appear in a logical sequence
diag_order <- c("LD", "DN", "FSGS", "FSGS+MCD", "MCD", "RPGN", "TN")
diag_order <- intersect(diag_order, unique(expr_df$Diagnosis))
expr_df$Diagnosis <- factor(expr_df$Diagnosis, levels = diag_order)

long_df <- tidyr::pivot_longer(
  expr_df,
  cols      = all_of(genes_of_interest),
  names_to  = "Gene",
  values_to = "Expression"
)

# ---- 6. Color palette (one colour per diagnosis group) ----
palette <- c(
  LD         = "#4DAF4A",
  DN         = "#E41A1C",
  FSGS       = "#377EB8",
  `FSGS+MCD` = "#984EA3",
  MCD        = "#FF7F00",
  RPGN       = "#A65628",
  TN         = "#F781BF"
)

# ---- 7. Plot function: one PNG per gene ----
plot_gene <- function(gene_id, data, output_dir = ".") {
  
  df <- data[data$Gene == gene_id, ]
  
  p <- ggplot(df, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2.5,
                 outlier.fill  = "white", outlier.stroke = 0.5,
                 width = 0.6, linewidth = 0.55) +
    geom_jitter(width = 0.18, size = 1.2, alpha = 0.55, colour = "grey30") +
    scale_fill_manual(values = palette) +
    labs(
      title    = paste("Gene Expression \u2013", gene_id),
      subtitle = "GSE104954  |  Kidney tubulointerstitium  (GPL22945 platform)",
      x        = "Diagnosis Group",
      y        = "RMA-normalised expression (log2)",
      fill     = "Diagnosis"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title         = element_text(face = "bold", size = 15),
      plot.subtitle      = element_text(colour = "grey40", size = 10),
      axis.text.x        = element_text(angle = 30, hjust = 1, size = 11),
      panel.grid.major.x = element_blank(),
      legend.position    = "right"
    )
  
  out_file <- file.path(output_dir, paste0("boxplot_", gene_id, ".png"))
  ggsave(out_file, plot = p, width = 9, height = 6, dpi = 150)
  message("Saved: ", out_file)
  invisible(p)
}

# ---- 8. Generate one image per gene ----
output_dir <- "."    # change to your preferred output folder

for (gene in genes_of_interest) {
  plot_gene(gene, long_df, output_dir = output_dir)
}

message("\nDone! Five PNG files written to: ", normalizePath(output_dir))
