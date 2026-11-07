# ============================================================
# Gene Expression Boxplots from GEO Series Matrix
### Genes: 3026_at, 3426_at, 3075_at, 1305_at, 725_at
# 3934_at, 26762_at, 5329_at, 7132_at, 7133_at, 337_at, 922_at, 3486_at
# ============================================================

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsignif)

# ---- 1. Resolve script directory ----
script_dir <- if (rstudioapi::isAvailable()) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  dirname(normalizePath(sys.frames()[[1]]$ofile))
}

drawBoxplot <- function(GSE_ID_EXT) {

  FILE <- file.path(script_dir, paste0(GSE_ID_EXT, "_series_matrix.txt"))

  cat("Reading GEO series matrix ...\n")
  lines <- readLines(FILE)

  # Parse !Sample_* metadata lines → data frame (rows = samples, cols = fields)
  meta_lines <- grep("^!Sample_", lines, value = TRUE)
  meta_list  <- lapply(meta_lines, function(ln) {
    parts <- strsplit(ln, "\t")[[1]]
    list(key  = sub("^!Sample_", "", parts[1]),
         vals = gsub('^"|"$', "", parts[-1]))
  })
  meta_keys <- make.unique(sapply(meta_list, `[[`, "key"))
  pdata     <- as.data.frame(
    setNames(lapply(meta_list, `[[`, "vals"), meta_keys),
    stringsAsFactors = FALSE
  )

  # Parse expression table between table_begin / table_end markers
  tbl_start <- which(trimws(lines) == "!series_matrix_table_begin") + 1
  tbl_end   <- which(trimws(lines) == "!series_matrix_table_end")   - 1
  expr      <- as.matrix(read.table(
    text        = paste(lines[tbl_start:tbl_end], collapse = "\n"),
    header      = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, quote = '"'
  ))

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
  # Patterns match both GSE104954 and GSE104948 label variants (case-insensitive)
  diag_vals <- gsub("diabetic nephropathy",                          "DN",       diag_vals, ignore.case = TRUE)
  diag_vals <- gsub("focal segmental.*(minimal change|mcd)",         "FSGS+MCD", diag_vals, ignore.case = TRUE)
  diag_vals <- gsub("focal segmental",                               "FSGS",     diag_vals, ignore.case = TRUE)
  diag_vals <- gsub("minimal change disease",                        "MCD",      diag_vals, ignore.case = TRUE)
  diag_vals <- gsub("anca[- ]associated vasculitis",                 "RPGN",     diag_vals, ignore.case = TRUE)
  diag_vals <- gsub("tumor nephrectomy",                             "TN",       diag_vals, ignore.case = TRUE)

  # ---- 3. Select the five genes of interest ----
  genes_of_interest <- c("3934_at", "26762_at", "5329_at", "7132_at", "7133_at", "337_at", "922_at", "3486_at")
  missing <- setdiff(genes_of_interest, rownames(expr))
  if (length(missing) > 0)
    warning("Gene(s) not found in ", GSE_ID_EXT, ", skipping: ", paste(missing, collapse = ", "))
  genes_of_interest <- intersect(genes_of_interest, rownames(expr))
  if (length(genes_of_interest) == 0) {
    warning("No genes found in ", GSE_ID_EXT, ", skipping dataset.")
    return(invisible(NULL))
  }

  expr_sub <- expr[genes_of_interest, , drop = FALSE]

  # ---- 4. Build a tidy data frame ----
  expr_df <- as.data.frame(t(expr_sub))
  expr_df$Sample    <- sample_titles
  expr_df$Diagnosis <- diag_vals
  # Ordered factor so boxes appear in a logical sequence
  diag_order <- c("LD", "DN", "FSGS", "FSGS+MCD", "MCD", "RPGN", "TN")
  diag_order <- intersect(diag_order, unique(expr_df$Diagnosis))
  expr_df$Diagnosis <- factor(expr_df$Diagnosis, levels = diag_order)
  expr_df <- expr_df[!is.na(expr_df$Diagnosis), ]

  long_df <- tidyr::pivot_longer(
    expr_df,
    cols      = all_of(genes_of_interest),
    names_to  = "Gene",
    values_to = "Expression"
  )

  # ---- 5. Color palette (one colour per diagnosis group) ----
  palette <- c(
    LD         = "#4daf4a",
    DN         = "#E41A1C",
    FSGS       = "#377EB8",
    `FSGS+MCD` = "#984EA3",
    MCD        = "#FF7F00",
    RPGN       = "#A65628",
    TN         = "#F781BF"
  )

  # ---- 6. Output folder ----
  out_dir <- file.path(script_dir, "Boxplots_Other")
  if (!dir.exists(out_dir)) dir.create(out_dir)

  # ---- 7. Plot function: one PNG per gene ----
  gene_labels <- c(
    # "3026_at" = "Habp2",
    # "3426_at" = "Cfi",
    # "3075_at" = "Cfh",
    # "1305_at" = "Col13a1",
    # "725_at"  = "C4bpb",
    "3934_at" = "LCN2",
    "26762_at" = "HAVCR1",
    "5329_at" = "PLAUR",
    "7132_at" = "TNFRSF1A",
    "7133_at" = "TNFRSF1B",
    "337_at" = "APOA4",
    "922_at" = "CD5L",
    "3486_at" = "IGFBP3"
  )

  plot_gene <- function(gene_id, data) {

    df         <- data[data$Gene == gene_id, ]
    gene_label <- gene_labels[[gene_id]]

    # ---- Wilcoxon p-values vs LD reference ----
    diag_levels <- levels(df$Diagnosis)
    ref_vals    <- df$Expression[df$Diagnosis == "LD"]
    y_max       <- max(df$Expression, na.rm = TRUE)
    y_range     <- y_max - min(df$Expression, na.rm = TRUE)

    fmt_pval <- function(p) {
      if (p < 0.001) "p<0.001"
      else if (p < 0.01) paste0("p=", formatC(p, digits = 3, format = "f"))
      else paste0("p=", formatC(p, digits = 2, format = "f"))
    }

    annot_rows <- lapply(setdiff(diag_levels, "LD"), function(grp) {
      vals <- df$Expression[df$Diagnosis == grp]
      if (length(vals) < 2 || length(ref_vals) < 2) return(NULL)
      pv <- wilcox.test(vals, ref_vals, exact = FALSE)$p.value
      data.frame(
        grp   = grp,
        x_grp = as.numeric(factor(grp, levels = diag_levels)),
        label = fmt_pval(pv),
        stringsAsFactors = FALSE
      )
    })
    annot_df <- do.call(rbind, Filter(Negate(is.null), annot_rows))

    # Stack brackets by span: shortest (closest to LD) sits lowest
    y_step <- 0.12 * y_range
    y_base <- y_max + 0.08 * y_range

    if (!is.null(annot_df) && nrow(annot_df) > 0) {
      annot_df <- annot_df[order(annot_df$x_grp), ]
      annot_df$y_brack <- y_base + (seq_len(nrow(annot_df)) - 1) * y_step
    }

    y_top <- if (!is.null(annot_df) && nrow(annot_df) > 0)
      max(annot_df$y_brack) + 0.10 * y_range else y_max + 0.10 * y_range

    comparisons <- lapply(annot_df$grp, function(g) c("LD", g))

    p <- ggplot(df, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
      geom_boxplot(outlier.shape = 16, outlier.size = 1.2,
                   outlier.colour = "grey30", outlier.alpha = 0.55,
                   width = 0.6, linewidth = 0.55) +
      geom_jitter(width = 0.12, size = 2.3, alpha = 0.7) +
      geom_signif(
        comparisons = comparisons,
        annotations = annot_df$label,
        y_position  = annot_df$y_brack,
        tip_length  = 0.01,
        textsize    = 3.5,
        vjust       = -0.2,
        color       = "grey35"
      ) +
      scale_fill_manual(values = palette) +
      coord_cartesian(ylim = c(NA, y_top)) +
      labs(
        title    = paste("Gene Expression \u2013", gene_label),
        subtitle = paste0("Dataset: ", GSE_ID_EXT),
        x        = "Treatment",
        # y        = expression("RMA-normalised expression (log"[2]*")"),
        y        = expression("RMA-normalised expression (log"[2]*")"),
        # fill     = "Diagnosis"
      ) +
      theme_bw(base_size = 13) +
      theme(
      axis.text.x        = element_text(angle = 30, hjust = 1, face = "bold"),
      legend.position    = "none",
      plot.title         = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
      # theme(
      #   plot.title         = element_text(face = "bold", size = 15),
      #   plot.subtitle      = element_text(colour = "grey40", size = 10),
      #   axis.text.x        = element_text(angle = 30, hjust = 1, size = 11),
      #   panel.grid.major.x = element_blank(),
      #   legend.position    = "right"
      # )

    out_file <- file.path(out_dir, paste0(GSE_ID_EXT, "_boxplot_", gene_label, ".png"))
    ggsave(out_file, plot = p, width = 8, height = 6, dpi = 300)
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
