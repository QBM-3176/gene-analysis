#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — per-dataset RenoRisk boxplots and heatmaps
#
# One boxplot PNG per gene per dataset, following the two Examples/ templates:
#   * microarray -> Examples/"boxplot_rr_ affimetrix.r"
#       series matrix -> panel probes -> boxplot by diagnosis group,
#       Wilcoxon vs the cohort's control/reference arm, geom_signif brackets
#   * RNA-Seq    -> Examples/boxplot_rr_dnaseq.r
#       counts -> DESeqDataSetFromMatrix -> DESeq -> vst -> same boxplot, with
#       DESeq2 Wald padj instead of Wilcoxon
#
# Each dataset is drawn STANDALONE — no batch correction, no pooling. The y axis
# is whatever that submitter deposited, which differs between cohorts; `y_lab`
# in the registry states the scale per dataset so no plot overstates its units.
#
# Plus two heatmap PNGs per dataset, following Examples/heatmap.r (Z-score per
# gene, samples annotated by group): the top `cfg$heatmap_top` most variable
# probes transcriptome-wide, and the 5 RenoRisk panel genes alone.
#
# Output: Results/SingleDatasetPlots/<Raw_filtered subfolder name>/<cohort>_boxplot_<GENE>.png
#         Results/SingleDatasetPlots/<Raw_filtered subfolder name>/<cohort>_heatmap_top<N>.png
#         Results/SingleDatasetPlots/<Raw_filtered subfolder name>/<cohort>_heatmap_5genes.png
#
# NOTE ON THE RNA-SEQ BRANCH
#   GSE142025 (RNA-Seq; user instruction 2026-11-03) is treated like a
#   microarray: its registry entry names a submitter-normalised log2 table
#   (`expr_file`, rows = gene symbols) and it goes through process_microarray().
#   `process_rnaseq()` (raw counts -> DESeq2/VST, Wald padj) fires for an
#   RNA-Seq series with no expression table, or with `counts = TRUE`:
#   GSE141295 and GSE175759 (<id>_counts.txt, by instruction 2026-11-04).
# =============================================================================

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
  min_n       = 2,       # a group needs this many samples to enter a Wilcoxon
  heatmap_top = 50,      # most variable probes per dataset heatmap
  width       = 8,
  height      = 6,
  dpi         = 300
)
cfg$raw_dir <- file.path(cfg$project_dir, "Data", "Raw_filtered")
cfg$results_dir <- file.path(cfg$project_dir, "Results")   # everything the pipeline generates
cfg$out_dir <- file.path(cfg$results_dir, "SingleDatasetPlots")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggsignif); library(AnnotationDbi)
  library(WGCNA); library(writexl)
})
set.seed(1)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

ENTREZ <- c(C4BPA = "722", COL13A1 = "1305", HABP2 = "3026", CFH = "3075", CFI = "3426")

# Source technology shown in each plot title, keyed off the annotation route.
# NOT all microarray cohorts are Affymetrix: GSE116626 is an Illumina BeadChip and
# GSE66494 an Agilent array, so labelling either "Affymetrix" would be false.
# The RNA-Seq path overrides this with "DNA-Seq".
PLATFORM <- c(entrez                      = "Affymetrix",
              hta20transcriptcluster.db   = "Affymetrix",
              hgu133a2.db                 = "Affymetrix",
              illuminaHumanv4.db          = "Illumina",
              hgug4112a.db                = "Agilent",
              symbol                      = "RNA-Seq")

# ---------------------------------------------------------------------------
# 1. Label abbreviation — one shared vocabulary across all 16 microarray cohorts + GSE142025.
#    Order matters: the FSGS+MCD compound must be matched before plain FSGS.
# ---------------------------------------------------------------------------
ABBREV <- list(
  c("glomerulus of control kidney",                      "Ctrl-Glom"),
  c("tubuli of control kidney",                          "Ctrl-Tub"),
  c("glomerulus of dkd kidney",                          "DKD-Glom"),
  c("tubuli of dkd kidney",                              "DKD-Tub"),
  c("focal segmental.*(minimal change|mcd)",             "FSGS+MCD"),
  c("focal (segmental|and segmental).*(sclerosis|scelrosis)", "FSGS"),
  c("minimal_lesion_igan",                               "IgAN-min"),
  c("chronic_lesion_igan",                               "IgAN-chr"),
  c("mixed_lesions_igan",                                "IgAN-mix"),
  c("active_lesion_igan",                                "IgAN-act"),
  c("non.?igan gn",                                      "non-IgAN"),
  c("iga nephropathy|from human igan$|immunoglobulin a nephropathy", "IgAN"),
  c("minimal change disease",                            "MCD"),
  c("anca.?associated vasculitis",                       "RPGN"),
  c("membranous (nephropathy|glomerulonephritis)",       "MN"),
  c("systemic lupus|lupus nephritis",                    "SLE"),
  c("hypertensive nephropathy",                          "HTN"),
  c("thin membrane disease",                             "TMD"),
  c("early diabetic nephropathy",                        "Early_DN"),
  c("non-diabetic control|from human control$",          "Control"),
  c("diabetic nephropathy",                              "DN"),
  c("tumor nephrectom|normal cortex from renal cell carcinoma", "TN"),
  c("living donor|^ld$",                                 "LD"),
  c("glomerular disease, other",                         "Other"),
  c("normal kidney",                                     "Normal"),
  c("renal biopsy specimens",                            "CKD"),
  c("^control$",                                         "Control")
)
abbreviate_groups <- function(x) {
  v <- tolower(trimws(x))
  out <- rep(NA_character_, length(v))
  for (rule in ABBREV) {
    hit <- is.na(out) & grepl(rule[1], v)
    out[hit] <- rule[2]
  }
  ifelse(is.na(out), x, out)   # unmatched labels pass through verbatim
}

PALETTE <- c(
  LD = "#4DAF4A", Control = "#4DAF4A", Normal = "#4DAF4A", TN = "#66C2A5",
  `Ctrl-Glom` = "#4DAF4A", `Ctrl-Tub` = "#1B9E77",
  DN = "#E41A1C", `DKD-Glom` = "#E41A1C", `DKD-Tub` = "#B2182B", CKD = "#E41A1C",
  FSGS = "#377EB8", `FSGS+MCD` = "#984EA3", MCD = "#FF7F00", RPGN = "#A65628",
  MN = "#F781BF", SLE = "#999999", HTN = "#8DA0CB", TMD = "#E7298A",
  IgAN = "#377EB8", `IgAN-min` = "#80B1D3", `IgAN-chr` = "#FB8072",
  `IgAN-mix` = "#BEBADA", `IgAN-act` = "#FDB462", `non-IgAN` = "#B3DE69",
  Other = "#CCCCCC", Early_DN = "#FB8072", Advanced_DN = "#E41A1C"
)

# ---------------------------------------------------------------------------
# 2. Dataset registry.  `dir` is the Raw_filtered subfolder (and the output
#    subfolder name); `id` is the series-matrix basename.
#    `grp` names where the grouping labels live; `ref` is the reference arm for
#    the significance brackets; `pairs` overrides it where a single reference is
#    wrong (GSE30122 carries two compartments, each with its own control).
# ---------------------------------------------------------------------------
COHORTS <- list(
  list(dir = "GSE104948", id = "GSE104948-GPL22945", idtype = "entrez",
       grp = c("char", "diagnosis"), blank = "LD", ref = "LD",
       y_lab = "RMA-normalised expression (log2)"),
  list(dir = "GSE104954", id = "GSE104954-GPL22945", idtype = "entrez",
       grp = c("char", "diagnosis"), blank = "LD", ref = "LD",
       y_lab = "RMA-normalised expression (log2)"),
  list(dir = "GSE108109", id = "GSE108109", idtype = "entrez",
       grp = c("source", NA), ref = "LD",
       y_lab = "Normalised expression (log2)"),
  list(dir = "GSE200828", id = "GSE200828", idtype = "entrez",
       grp = c("source", NA), ref = "LD",
       y_lab = "RMA-normalised expression (log2)"),
  list(dir = "GSE99339", id = "GSE99339-GPL19184", idtype = "entrez",
       grp = c("desc", NA), ref = "TN",
       y_lab = "Normalised expression (log2)"),
  # Added 2026-11-03 (by instruction).
  list(dir = "GSE93798", id = "GSE93798", idtype = "entrez",
       grp = c("source", NA), ref = "Control",
       y_lab = "RMA-normalised expression (log2)"),
  list(dir = "GSE35488", id = "GSE35488", idtype = "entrez",
       grp = c("source", NA), ref = "LD",
       y_lab = "RMA-normalised expression (log2)"),
  list(dir = "GSE111154", id = "GSE111154", idtype = "hta20transcriptcluster.db",
       grp = c("source", NA), ref = "Control",
       y_lab = "Normalised expression (log2, Partek)"),
  list(dir = "GSE121233", id = "GSE121233", idtype = "hta20transcriptcluster.db",
       grp = c("char", "disease state"), ref = "Control",
       y_lab = "RMA-normalised expression (log2)"),
  list(dir = "GSE129973", id = "GSE129973", idtype = "hta20transcriptcluster.db",
       grp = c("source", NA), ref = "TN",
       y_lab = "RMA-normalised expression (log2)"),
  list(dir = "GSE96804", id = "GSE96804", idtype = "hta20transcriptcluster.db",
       grp = c("source", NA), ref = "TN",
       y_lab = "RMA-normalised expression (log2)"),
  # NB: GSE30122 carries BOTH compartments, so reference-vs-all compares
  # Ctrl-Glom against tubular groups too. Kept for consistency with the other 12
  # datasets; restore `pairs` below to confine comparisons within compartment:
  #   pairs = list(c("Ctrl-Glom", "DKD-Glom"), c("Ctrl-Tub", "DKD-Tub")),
  list(dir = "GSE30122", id = "GSE30122", idtype = "hgu133a2.db",
       grp = c("source", NA), ref = "Ctrl-Glom",
       y_lab = "Per-gene centred expression (log2 ratio)"),
  list(dir = "GSE116626", id = "GSE116626", idtype = "illuminaHumanv4.db",
       grp = c("source", NA), ref = "LD",
       y_lab = "Per-gene centred expression (log2 ratio)"),
  list(dir = "GSE66494", id = "GSE66494", idtype = "hgug4112a.db",
       grp = c("char", "sample type"), ref = "Normal", log2 = TRUE,
       y_lab = "log2 expression (Agilent, linear source)"),
  # RNA-Seq: the series matrix is metadata-only; values come from the
  # submitter's table (rows = gene symbols, already log2). C4BPA is absent.
  list(dir = "GSE142025", id = "GSE142025", idtype = "symbol",
       expr_file = "GSE142025_expression.txt",
       grp = c("char", "group"), ref = "Control",
       y_lab = "Submitter-normalised log2 (RNA-Seq)"),
  # RNA-Seq from raw counts (`counts = TRUE`): DESeq2 branch on <id>_counts.txt;
  # its VST equals <id>_expression.txt (both made the same way, rnaseq_counts_vst.R).
  list(dir = "GSE141295", id = "GSE141295", idtype = "symbol",
       expr_file = "GSE141295_expression.txt", counts = TRUE,
       grp = c("char", "tissue1"), ref = "TN",
       y_lab = "DESeq2 VST (RNA-Seq, from counts)"),
  list(dir = "GSE175759", id = "GSE175759", idtype = "symbol",
       expr_file = "GSE175759_expression.txt", counts = TRUE,
       grp = c("char", "diagnosis"), ref = "Control",
       y_lab = "DESeq2 VST (RNA-Seq, from counts)")
)

# ---------------------------------------------------------------------------
# 3. Series-matrix reader (metadata + table; no GPL fetch)
# ---------------------------------------------------------------------------
read_sm <- function(path) {
  lines <- readLines(path, warn = FALSE)
  meta  <- grep("^!Sample_", lines, value = TRUE)
  parts <- lapply(meta, function(ln) {
    p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    list(key = sub("^!Sample_", "", p[1]), vals = gsub('^"|"$', "", p[-1]))
  })
  keys  <- make.unique(vapply(parts, `[[`, "", "key"))
  pdata <- as.data.frame(setNames(lapply(parts, `[[`, "vals"), keys),
                         stringsAsFactors = FALSE)
  rownames(pdata) <- pdata$geo_accession

  b <- which(trimws(lines) == "!series_matrix_table_begin") + 1
  e <- which(trimws(lines) == "!series_matrix_table_end")   - 1
  expr <- if (e >= b)
    as.matrix(read.table(text = paste(lines[b:e], collapse = "\n"),
                         header = TRUE, sep = "\t", row.names = 1,
                         check.names = FALSE, quote = '"', comment.char = ""))
  else matrix(numeric(0), 0, 0)
  list(pheno = pdata, expr = expr)
}

is_rnaseq <- function(pdata) {
  s <- unlist(pdata[, grep("^library_strategy", names(pdata)), drop = FALSE])
  any(grepl("RNA-Seq", s, ignore.case = TRUE))
}

group_labels <- function(pdata, co) {
  kind <- co$grp[1]; fld <- co$grp[2]
  raw <- if (kind == "source") pdata$source_name_ch1
         else if (kind == "desc") pdata$description
         else {
           cols <- grep("^characteristics_ch1", names(pdata), value = TRUE)
           hit  <- cols[vapply(cols, function(c)
                     any(grepl(paste0("^", fld, ":"), pdata[[c]], ignore.case = TRUE)), TRUE)]
           sub(paste0("^", fld, ":\\s*"), "", pdata[[hit[1]]], ignore.case = TRUE)
         }
  raw <- trimws(raw)
  # ERCB encodes living donors as an empty diagnosis value.
  if (!is.null(co$blank)) raw[raw == ""] <- co$blank
  abbreviate_groups(raw)
}

# ---------------------------------------------------------------------------
# 4. Panel extraction (shared with Stage 2/3): probe -> SYMBOL, MaxMean collapse
#    via WGCNA::collapseRows
# ---------------------------------------------------------------------------
panel_probes <- function(idtype, probes) {
  if (idtype == "symbol") {
    hit <- intersect(cfg$genes, probes)
    return(setNames(hit, hit))
  }
  if (idtype == "entrez") {   # "<entrez>_at", or bare "<entrez>" (GSE35488)
    want <- c(paste0(ENTREZ, "_at"), ENTREZ)
    hit  <- intersect(want, probes)
    return(setNames(rep(names(ENTREZ), 2)[match(hit, want)], hit))
  }
  if (!requireNamespace(idtype, quietly = TRUE))
    stop("annotation package not installed: ", idtype)
  db  <- getExportedValue(idtype, idtype)   # ChipDb keeps the .db suffix
  sym <- suppressMessages(AnnotationDbi::mapIds(db, keys = probes, column = "SYMBOL",
                                                keytype = "PROBEID", multiVals = "first"))
  sym <- sym[!is.na(sym) & sym %in% cfg$genes]
  setNames(unname(sym), names(sym))
}

extract_panel <- function(expr, co) {
  map <- panel_probes(co$idtype, rownames(expr))
  if (!length(map)) return(NULL)
  sub <- expr[names(map), , drop = FALSE]
  if (isTRUE(co$log2)) sub <- log2(sub + 1)
  # Multiple probes per gene -> keep the one with the highest mean expression;
  # rows of the result are named by gene symbol.
  WGCNA::collapseRows(sub, rowGroup = unname(map[rownames(sub)]),
                      rowID = rownames(sub), method = "MaxMean")$datETcollapsed
}

# ---------------------------------------------------------------------------
# 5. The plot — one gene, one dataset
# ---------------------------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "p<0.001"
  else if (p < 0.01) paste0("p=", formatC(p, digits = 3, format = "f"))
  else paste0("p=", formatC(p, digits = 2, format = "f"))
}

draw_gene <- function(df, gene, co, out_dir, pfun) {
  lv  <- levels(df$Group)
  cmp <- if (!is.null(co$pairs)) co$pairs else lapply(setdiff(lv, co$ref), function(g) c(co$ref, g))
  cmp <- Filter(function(p) all(p %in% lv) &&
                  all(vapply(p, function(g) sum(df$Group == g) >= cfg$min_n, TRUE)), cmp)

  ann <- vapply(cmp, function(p) fmt_p(pfun(df, p[1], p[2], gene)), character(1))
  rng <- range(df$Expression, na.rm = TRUE); span <- diff(rng)
  if (!is.finite(span) || span == 0) span <- 1
  ypos <- rng[2] + span * 0.10 * seq_along(cmp)

  p <- ggplot(df, aes(Group, Expression, fill = Group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1.2, outlier.colour = "grey30",
                 outlier.alpha = 0.55, width = 0.6, linewidth = 0.55) +
    geom_jitter(width = 0.12, size = 2.3, alpha = 0.7) +
    scale_fill_manual(values = PALETTE, na.value = "grey70") +
    labs(title = paste("Gene Expression –", gene),
         subtitle = sprintf("Dataset: %s (%s)", co$id,
                            if (!is.null(co$src)) co$src else unname(PLATFORM[co$idtype])),
         x = "Group", y = co$y_lab) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          legend.position = "none", plot.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank())

  if (length(cmp))
    p <- p + geom_signif(comparisons = cmp, annotations = ann, y_position = ypos,
                         tip_length = 0.01, textsize = 3.5, vjust = -0.2,
                         colour = "grey35") +
             coord_cartesian(ylim = c(NA, max(ypos) + span * 0.10))

  ggsave(file.path(out_dir, paste0(co$id, "_boxplot_", gene, ".png")), p,
         width = cfg$width, height = cfg$height, dpi = cfg$dpi)
}

# ---------------------------------------------------------------------------
# 5b. The heatmap — one per dataset (Examples/heatmap.r): low-expression
#     prefilter, top-N most variable probes transcriptome-wide, samples as rows
#     annotated by group, Z-score per gene, euclidean/complete clustering.
#     Columns are labelled by gene symbol where the platform maps one (raw
#     probe IDs such as TC01000001.hg.1 are unreadable), else by probe ID.
# ---------------------------------------------------------------------------
probe_symbols <- function(idtype, probes) {
  if (idtype == "symbol") return(make.unique(probes))
  sym <- tryCatch(suppressMessages(
    if (idtype == "entrez")
      AnnotationDbi::mapIds(getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"),
                            keys = sub("_at$", "", probes), column = "SYMBOL",
                            keytype = "ENTREZID", multiVals = "first")
    else
      AnnotationDbi::mapIds(getExportedValue(idtype, idtype), keys = probes,
                            column = "SYMBOL", keytype = "PROBEID", multiVals = "first")),
    error = function(e) rep(NA_character_, length(probes)))
  make.unique(ifelse(is.na(sym), probes, unname(sym)))
}

# Shared renderer: `mat` is samples x genes. Column width and canvas grow with
# the content; the colour scale is fixed at 50 steps (the example's value) so a
# 5-gene heatmap is not drawn with only 5 colours.
render_heatmap <- function(mat, grp, co, out_file, title, cellwidth = 11) {
  ord <- c(co$ref, setdiff(names(sort(table(grp), decreasing = TRUE)), co$ref))
  annotation_row <- data.frame(Group = factor(grp, levels = intersect(ord, unique(grp))),
                               row.names = rownames(mat))
  lv   <- levels(annotation_row$Group)
  cols <- PALETTE[lv]
  cols[is.na(cols)] <- grDevices::hcl.colors(sum(is.na(cols)), "Dark 3")
  names(cols) <- lv

  cellheight <- 5
  # The example's fixed 12 x 12 in canvas clips the larger cohorts (up to 153
  # samples) and leaves the small ones mostly blank; size the height to the
  # rows plus ~7 in for title, column dendrogram and the rotated gene labels.
  height <- nrow(mat) * cellheight / 72 + 7
  # 12 in is the example's width (top-50); a few-gene heatmap gets a canvas
  # sized to its columns plus ~6 in for dendrogram, annotation and legend.
  width  <- if (ncol(mat) >= 20) 12 else ncol(mat) * cellwidth / 72 + 6
  n_col  <- 50
  pheatmap::pheatmap(mat,
           scale                    = "column",  # Z-score per gene
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method        = "complete",
           cellwidth  = cellwidth,
           cellheight = cellheight,
           color  = colorRampPalette(c("navy", "white", "red"))(n_col),
           breaks = seq(-1.5, 1.5, length.out = n_col + 1),
           annotation_row    = annotation_row,
           annotation_colors = list(Group = cols),
           show_rownames = FALSE,
           show_colnames = TRUE,
           fontsize      = 14,
           fontsize_row  = 11,
           fontsize_col  = 11,
           main = title,
           filename = out_file,
           width = width, height = height)
}

# Top-N most variable probes, transcriptome-wide.
draw_heatmap <- function(expr, grp, co, out_dir, n_top = cfg$heatmap_top) {
  expr <- expr[stats::complete.cases(expr), , drop = FALSE]   # hclust cannot take NA
  keep <- rowSums(expr > median(expr)) >= 0.10 * ncol(expr)
  expr <- expr[keep, , drop = FALSE]
  n_top <- min(n_top, nrow(expr))
  top  <- names(sort(apply(expr, 1, var), decreasing = TRUE))[seq_len(n_top)]

  mat <- t(expr[top, , drop = FALSE])                # samples x genes
  colnames(mat) <- probe_symbols(co$idtype, top)
  render_heatmap(mat, grp, co,
                 file.path(out_dir, paste0(co$id, "_heatmap_top", n_top, ".png")),
                 sprintf("%s\nTop %d Variable Genes\n%s", co$id, n_top, co$y_lab))
  n_top
}

# The RenoRisk panel only: `panel` is the collapsed genes x samples matrix from
# extract_panel(). Wider cells — five 11 pt columns would be slivers.
draw_panel_heatmap <- function(panel, grp, co, out_dir) {
  panel <- panel[intersect(cfg$genes, rownames(panel)), , drop = FALSE]
  ok    <- stats::complete.cases(t(panel))
  mat   <- t(panel[, ok, drop = FALSE])               # samples x genes
  render_heatmap(mat, grp[ok], co,
                 file.path(out_dir, paste0(co$id, "_heatmap_", nrow(panel), "genes.png")),
                 sprintf("%s\nRenoRisk panel (%d genes)\n%s", co$id, nrow(panel), co$y_lab),
                 cellwidth = 40)
}

# ---------------------------------------------------------------------------
# 6a. Microarray path  (Examples/"boxplot_rr_ affimetrix.r")
# ---------------------------------------------------------------------------
process_microarray <- function(sm, co, out_dir) {
  mat <- extract_panel(sm$expr, co)
  if (is.null(mat)) stop("none of the panel genes map on this platform")
  grp <- group_labels(sm$pheno, co)[match(colnames(mat), rownames(sm$pheno))]

  ord <- c(co$ref, setdiff(names(sort(table(grp), decreasing = TRUE)), co$ref))
  wilcox_p <- function(d, a, b, ...) {
    x <- d$Expression[d$Group == a]; y <- d$Expression[d$Group == b]
    if (length(x) < cfg$min_n || length(y) < cfg$min_n) return(NA_real_)
    suppressWarnings(wilcox.test(y, x, exact = FALSE)$p.value)
  }
  for (g in rownames(mat)) {
    df <- data.table(Expression = as.numeric(mat[g, ]),
                     Group = factor(grp, levels = intersect(ord, unique(grp))))
    df <- df[!is.na(Group) & !is.na(Expression)]
    draw_gene(df, g, co, out_dir, wilcox_p)
  }
  draw_panel_heatmap(mat, grp, co, out_dir)
  full <- if (isTRUE(co$log2)) log2(sm$expr + 1) else sm$expr
  draw_heatmap(full, group_labels(sm$pheno, co)[match(colnames(full), rownames(sm$pheno))],
               co, out_dir)
  miss <- setdiff(cfg$genes, rownames(mat))
  data.table(genes = nrow(mat), samples = ncol(mat), groups = uniqueN(grp),
             test = "Wilcoxon vs reference",
             missing_genes = if (length(miss)) paste(miss, collapse = ", ") else NA_character_)
}

# ---------------------------------------------------------------------------
# 6b. RNA-Seq path  (Examples/boxplot_rr_dnaseq.r) — DESeq2 -> VST -> Wald padj.
#     Needs a counts matrix in the dataset folder; a GEO series matrix for an
#     RNA-Seq series carries an empty expression table, so one must be supplied.
# ---------------------------------------------------------------------------
process_rnaseq <- function(sm, co, out_dir, folder) {
  co$src <- "DNA-Seq"        # shown in the plot titles instead of an array vendor
  cf <- list.files(folder, "count.*\\.(csv|tsv|txt)(\\.gz)?$", full.names = TRUE,
                   ignore.case = TRUE)
  if (!length(cf))
    stop("RNA-Seq cohort but no counts matrix in ", basename(folder),
         " - the series matrix has no expression table; supply counts or run a ",
         "quantification stage first")
  for (p in c("DESeq2", "SummarizedExperiment"))
    if (!requireNamespace(p, quietly = TRUE)) stop("package not installed: ", p)

  raw <- fread(cf[1])
  gid <- as.character(raw[[grep("gene|geneid|ensembl", names(raw), ignore.case = TRUE)[1]]])
  num <- names(raw)[vapply(raw, is.numeric, TRUE)]
  cnt <- round(as.matrix(raw[, ..num])); rownames(cnt) <- gid

  grp <- group_labels(sm$pheno, co)[match(colnames(cnt), rownames(sm$pheno))]
  cd  <- data.frame(row.names = colnames(cnt), condition = factor(grp))
  dds <- DESeq2::DESeqDataSetFromMatrix(cnt, cd, design = ~ condition)
  dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
  # VST shifts slightly with factor-level names/order; plot the analysis values instead.
  if (nrow(sm$expr)) vsd <- sm$expr[rownames(vsd), colnames(vsd), drop = FALSE]

  # Rows may be Ensembl or symbol; keep whichever carries the panel.
  sym  <- if (any(cfg$genes %in% rownames(vsd))) rownames(vsd) else
          suppressMessages(AnnotationDbi::mapIds(
            getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"), rownames(vsd),
            "SYMBOL", "ENSEMBL", multiVals = "first"))
  hit <- which(sym %in% cfg$genes)
  if (!length(hit)) stop("no panel gene survived DESeq2 prefiltering")

  wald_p <- function(d, a, b, gene) {
    r <- DESeq2::results(dds, contrast = c("condition", b, a))
    rn <- rownames(vsd)[hit][match(gene, sym[hit])]
    if (rn %in% rownames(r)) r[rn, "padj"] else NA_real_
  }
  ord <- c(co$ref, setdiff(names(sort(table(grp), decreasing = TRUE)), co$ref))
  for (i in hit) {
    g  <- sym[i]
    df <- data.table(Expression = as.numeric(vsd[i, ]),
                     Group = factor(grp, levels = intersect(ord, unique(grp))))
    df <- df[!is.na(Group) & !is.na(Expression)]
    draw_gene(df, g, co, out_dir, wald_p)
  }
  draw_heatmap(vsd, grp, co, out_dir)
  pnl <- vsd[hit, , drop = FALSE]; rownames(pnl) <- sym[hit]
  draw_panel_heatmap(pnl, grp, co, out_dir)
  data.table(genes = length(hit), samples = ncol(vsd), groups = uniqueN(grp),
             test = "DESeq2 Wald padj")
}

# ---------------------------------------------------------------------------
# 7. Driver — one Results/SingleDatasetPlots/<subfolder> per Raw_filtered dataset
# ---------------------------------------------------------------------------
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
report <- list()

for (co in COHORTS) {
  folder <- file.path(cfg$raw_dir, co$dir)
  f <- list.files(folder, paste0("^", co$id, "_series_matrix\\.txt$"), full.names = TRUE)
  if (!length(f)) { warning(co$id, ": series matrix not found"); next }

  out_dir <- file.path(cfg$out_dir, co$dir)    # named after the Raw_filtered folder
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  sm   <- read_sm(f[1])
  if (!is.null(co$expr_file)) {
    x <- fread(file.path(folder, co$expr_file))
    sm$expr <- as.matrix(x[, -1]); rownames(sm$expr) <- x[[1]]
    stopifnot(all(colnames(sm$expr) %in% rownames(sm$pheno)))
  }
  type <- if (is_rnaseq(sm$pheno)) "rnaseq" else "microarray"
  r <- tryCatch(
    # Normalised RNA-Seq values take the microarray path; DESeq2 only for counts.
    if (type == "rnaseq" && (!nrow(sm$expr) || isTRUE(co$counts))) process_rnaseq(sm, co, out_dir, folder)
    else                  process_microarray(sm, co, out_dir),
    error = function(e) { warning(co$id, ": ", conditionMessage(e)); conditionMessage(e) })

  report[[co$id]] <- if (is.data.table(r))
    data.table(dataset = co$dir, cohort = co$id, type, status = "done", r)
  else
    data.table(dataset = co$dir, cohort = co$id, type, status = "failed",
               genes = NA_integer_, samples = NA_integer_, groups = NA_integer_, test = r)
  msg(co$dir, " [", type, "] -> ", out_dir)
}

tab <- rbindlist(report, fill = TRUE)
print(tab)
write_xlsx(tab, file.path(cfg$out_dir, "boxplots_summary.xlsx"))
msg("PNGs written: ", length(list.files(cfg$out_dir, "\\.png$", recursive = TRUE)))
writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_dir, "boxplots_sessionInfo.txt"))
