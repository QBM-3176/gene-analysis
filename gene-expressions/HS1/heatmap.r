############################################################
# limma pipeline + Heatmap for GSE104954-GPL22945
# Input: GSE104954-GPL22945_series_matrix.txt
#        (Affymetrix microarray, RMA-normalized log2 — pre-processed)
# Analysis: limma (correct DESeq2-analogue for microarray data)
# Output: Results CSVs + heatmap "GSE104954-GPL22945_heatmap_topXX.png"
############################################################

# Install packages if needed
# install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("GEOquery")
# install.packages("pheatmap")

library(limma)
library(GEOquery)
library(pheatmap)

############################################################
# Load expression matrix
############################################################

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
}

drawHeatmap <- function(GSE_ID_EXT, TOP_GENES_COUNT = 50) {

  FILE <- paste0(GSE_ID_EXT, "_series_matrix.txt")

  cat("Reading series matrix …\n")
  gse   <- getGEO(filename = FILE, GSEMatrix = TRUE)
  expr  <- exprs(gse)          # probes × samples, already log2-normalized
  pdata <- pData(gse)

  cat(sprintf("Expression matrix: %d probes × %d samples\n", nrow(expr), ncol(expr)))

  ############################################################
  # Extract sample metadata (colData equivalent)
  ############################################################

  char_cols <- grep("^characteristics_ch1", colnames(pdata), value = TRUE)
  diag_col  <- char_cols[sapply(char_cols, function(col)
                 any(grepl("^diagnosis:", pdata[[col]])))]
  diag_vals <- sub("^diagnosis: ", "", as.character(pdata[[diag_col[1]]]))
  diag_vals[diag_vals == ""] <- "LD"

  # grepl + assignment replaces the whole string (gsub would only replace the
  # matched substring, leaving e.g. "FSGS glomerulosclerosis" instead of "FSGS")
  # FSGS_MCD must run before FSGS so combined cases are not overwritten
  diag_vals[grepl("diabetic nephropathy",                  diag_vals, ignore.case = TRUE)] <- "DN"
  diag_vals[grepl("focal segmental.*(minimal change|mcd)", diag_vals, ignore.case = TRUE)] <- "FSGS_MCD"
  diag_vals[grepl("focal segmental",                       diag_vals, ignore.case = TRUE)] <- "FSGS"
  diag_vals[grepl("minimal change disease",                diag_vals, ignore.case = TRUE)] <- "MCD"
  diag_vals[grepl("anca[- ]associated vasculitis",         diag_vals, ignore.case = TRUE)] <- "RPGN"
  diag_vals[grepl("tumor nephrectomy",                     diag_vals, ignore.case = TRUE)] <- "TN"
  diag_vals[grepl("hypertensive nephropathy",              diag_vals, ignore.case = TRUE)] <- "HT"
  diag_vals[grepl("iga nephropathy",                       diag_vals, ignore.case = TRUE)] <- "IgA"
  diag_vals[grepl("membranous glomerulonephro",            diag_vals, ignore.case = TRUE)] <- "MGN"
  diag_vals[grepl("systemic lupus erythematosus",          diag_vals, ignore.case = TRUE)] <- "SLE"
  diag_vals[grepl("thin mem",                              diag_vals, ignore.case = TRUE)] <- "TMD"

  level_order <- c("LD", "DN", "FSGS", "FSGS_MCD", "MCD", "RPGN", "TN",
                   "HT", "IgA", "MGN", "SLE", "TMD")
  diagnosis <- factor(diag_vals, levels = intersect(level_order, unique(diag_vals)))

  cat("Sample counts per group:\n")
  print(table(diagnosis))

  ############################################################
  # Filter low-expression probes (analogous to DESeq2 prefilter)
  ############################################################

  keep   <- rowSums(expr > median(expr)) >= 0.10 * ncol(expr)
  expr_f <- expr[keep, ]
  cat(sprintf("After filtering: %d probes retained\n", nrow(expr_f)))

  ############################################################
  # Build design matrix & fit limma model (analogous to DESeqDataSet + DESeq())
  ############################################################

  design <- model.matrix(~ 0 + diagnosis)
  colnames(design) <- levels(diagnosis)
  fit <- lmFit(expr_f, design)

  ############################################################
  # Heatmap: Top TOP_GENES_COUNT most variable genes
  ############################################################

  all_genes_sorted <- names(sort(apply(expr_f, 1, var), decreasing = TRUE))
  top_var_genes    <- all_genes_sorted[1:TOP_GENES_COUNT]

  # Subset expression matrix
  mat <- expr_f[top_var_genes, ]

  # Transpose matrix (genes become columns, samples become rows)
  mat <- t(mat)

  ############################################################
  # Annotation (analogous to annotation_row in RDaseq2.r)
  ############################################################

  annotation_row <- data.frame(
    Diagnosis = diagnosis,
    row.names = rownames(mat)
  )

  # Map sample IDs → diagnosis labels for row display
  row_labels <- setNames(as.character(diagnosis), rownames(mat))

  all_colors <- c(
    LD       = "#00FFFF",
    DN       = "#0000FF",
    FSGS     = "#FF00FF",
    FSGS_MCD = "#00FF00",
    MCD      = "#FFA500",
    RPGN     = "#FF0000",
    TN       = "#FFFF00",
    HT       = "#A0522D",
    IgA      = "#8B008B",
    MGN      = "#2E8B57",
    SLE      = "#FF69B4",
    TMD      = "#708090"
  )
  annotation_colors <- list(Diagnosis = all_colors[levels(diagnosis)])

  ############################################################
  # Draw & save heatmap (same style as RDaseq2.r)
  ############################################################

  out_dir  <- "Heatmaps"
  if (!dir.exists(out_dir)) dir.create(out_dir)
  out_file  <- file.path(out_dir, paste0(GSE_ID_EXT, "_heatmap_top", TOP_GENES_COUNT, ".png"))
  # cell_size <- 5
  cellwidth <- 11
  cellheight <- 5

  # Total width: cells + row dendrogram + row labels + annotation + legend
  # pheatmap honours cellwidth/cellheight and auto-sizes internally;
  # width/height are the canvas upper bound (in inches at 300 dpi).
  # 22 × 300 = 6600 px — ample room for legend text "FSGS_MCD".
  pheatmap(mat,
           scale                    = "column",  # Z-score per gene
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method        = "complete",
           cellwidth  = cellwidth,
           cellheight = cellheight,
           color = colorRampPalette(c("navy", "white", "red"))(TOP_GENES_COUNT),
           breaks = seq(-1.5, 1.5, length.out = TOP_GENES_COUNT + 1),
           annotation_row    = annotation_row,
           annotation_colors = annotation_colors,
           show_rownames = FALSE,
           show_colnames = TRUE,
           fontsize      = 14,
           fontsize_row  = 11,
           fontsize_col  = 11,
           main = paste0(GSE_ID_EXT, "\nTop ", TOP_GENES_COUNT,
                         " Variable Genes (Limma analysis)\n"),
           filename = out_file,
           width = 12, height = 12)

  cat(paste0("Heatmap saved as ", out_file, "\n"))
}

drawHeatmap("GSE104954-GPL22945")
drawHeatmap("GSE104954-GPL24120")
drawHeatmap("GSE104948-GPL22945")
drawHeatmap("GSE104948-GPL24120")
