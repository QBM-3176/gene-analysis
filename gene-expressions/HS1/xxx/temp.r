# ============================================================
# Gene Expression Boxplots from GEO Series Matrix
# Dataset: GSE104954-GPL22945
# Genes: 3026_at, 3426_at, 3075_at, 1305_at, 725_at
# ============================================================

library(GEOquery)
library(ggplot2)
library(tidyr)
library(dplyr)

# ---- 1. Load the file via getGEO ----
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
}

drawBoxplot <- function(GSE_ID_EXT) {

  FILE <- paste0(GSE_ID_EXT, "_series_matrix.txt")

  cat("Reading GEO series matrix ...\n")
  gse   <- getGEO(filename = FILE, GSEMatrix = TRUE)
  expr  <- exprs(gse)
  pdata <- pData(gse)

  cat(sprintf("Expression matrix: %d probes x %d samples\n", nrow(expr), ncol(expr)))

  # ---- 2. Extract sample titles and diagnosis labels ----
  sample_titles <- pdata$title

  # Dynamically find the characteristics column that contains "diagnosis:"
  char_cols <- grep("^characteristics_ch1", colnames(pdata), value = TRUE)
  diag_col  <- char_cols[sapply(char_cols, function(col)
                 any(grepl("^diagnosis:", pdata[[col]])))]
  diag_vals <- as.character(pdata[[diag_col[1]]])

  # Clean up diagnosis labels
  diag_vals <- sub("^diagnosis: ", "", diag_vals)
  diag_vals[diag_vals == ""] <- "LD"          # empty = living donor controls

  # Abbreviate long labels for readability on plots
  diag_vals <- gsub("Diabetic nephropathy",                                         "DN",       diag_vals)
  diag_vals <- gsub("Focal segmental glomerulosclerosis and minimal change disease", "FSGS+MCD", diag_vals)
  diag_vals <- gsub("Focal segmental glomerulosclerosis",                           "FSGS",     diag_vals)
  diag_vals <- gsub("Minimal change disease",                                       "MCD",      diag_vals)
  diag_vals <- gsub("ANCA-associated vasculitis",                                   "RPGN",     diag_vals)
  diag_vals <- gsub("Tumor nephrectomy",                                            "TN",       diag_vals)

  # ---- 3. Select the five genes of interest ----
  genes_of_interest <- c("3026_at", "3426_at", "3075_at", "1305_at", "725_at")
  missing <- setdiff(genes_of_interest, rownames(expr))
  if (length(missing) > 0) stop("Gene(s) not found: ", paste(missing, collapse = ", "))

  expr_sub <- expr[genes_of_interest, , drop = FALSE]

  # ---- 4. Build a tidy data frame ----
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

  # ---- 5. Color palette (one colour per diagnosis group) ----
  palette <- c(
    LD         = "#4DAF4A",
    DN         = "#E41A1C",
    FSGS       = "#377EB8",
    `FSGS+MCD` = "#984EA3",
    MCD        = "#FF7F00",
    RPGN       = "#A65628",
    TN         = "#F781BF"
  )

  # ---- 6. Output folder ----
  out_dir <- "Boxplots"
  if (!dir.exists(out_dir)) dir.create(out_dir)

  # ---- 7. Plot function: one PNG per gene ----
  gene_labels <- c(
    "3026_at" = "Habp2",
    "3426_at" = "Cfi",
    "3075_at" = "Cfh",
    "1305_at" = "Col13a1",
    "725_at"  = "C4bpb"
  )

  plot_gene <- function(gene_id, data) {

    df         <- data[data$Gene == gene_id, ]
    gene_label <- gene_labels[[gene_id]]

    p <- ggplot(df, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
      geom_boxplot(outlier.shape = 21, outlier.size = 2.5,
                   outlier.fill  = "white", outlier.stroke = 0.5,
                   width = 0.6, linewidth = 0.55) +
      geom_jitter(width = 0.18, size = 1.2, alpha = 0.55, colour = "grey30") +
      scale_fill_manual(values = palette) +
      labs(
        title    = paste("Gene Expression \u2013", gene_label),
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

    out_file <- file.path(out_dir, paste0(GSE_ID_EXT, "_boxplot_", gene_label, ".png"))
    ggsave(out_file, plot = p, width = 9, height = 6, dpi = 150)
    message("Saved: ", out_file)
    invisible(p)
  }

  # ---- 8. Generate one image per gene ----
  for (gene in genes_of_interest) {
    plot_gene(gene, long_df)
  }
}

drawBoxplot("GSE104954-GPL22945")
drawBoxplot("GSE104954-GPL24120")
drawBoxplot("GSE104948-GPL22945")
drawBoxplot("GSE104948-GPL24120")
