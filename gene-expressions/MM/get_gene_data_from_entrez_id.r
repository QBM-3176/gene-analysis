# Install once if needed:
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi"))

library(AnnotationDbi)
library(org.Hs.eg.db)   # human
library(org.Mm.eg.db)   # mouse

# Generic helper: Entrez ID -> Ensembl gene ID + symbol
entrez_to_ensembl <- function(entrez_ids, species = c("human", "mouse")) {
  species <- match.arg(species)
  db <- if (species == "human") org.Hs.eg.db else org.Mm.eg.db

  entrez_ids <- as.character(entrez_ids)

  AnnotationDbi::select(
    db,
    keys     = entrez_ids,
    keytype  = "ENTREZID",
    columns  = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
  )
}

# --- Examples ---

human_ids <- c("725")
entrez_to_ensembl(human_ids, species = "human")

# mouse_ids <- c("22059", "12189", "17869")
# entrez_to_ensembl(mouse_ids, species = "mouse")