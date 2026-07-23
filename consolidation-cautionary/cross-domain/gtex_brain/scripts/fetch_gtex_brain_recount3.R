#!/usr/bin/env Rscript
# Fetch GTEx BRAIN (n~2931) from recount3, score Hallmark pathways in R, and
# export a compact pathway matrix + brain-region labels for the PSF pipeline.
#
# Large-N (public, no controlled access) demonstration of the brain-region confound at scale:
# GTEx brain is neurotypical, 13 subregions. Scoring pathways in R keeps the
# export tiny (samples x 50 pathways) instead of a ~100MB gene matrix.
#
# Output: results/gtex_brain_pathway_scores.tsv  (samples x pathways + region col)
suppressMessages({library(recount3); library(SummarizedExperiment)})

HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)))
if (length(HERE) == 0 || HERE == "") HERE <- "."
PANEL <- file.path(HERE, "../../../panels/hallmark_200genes.gmt")
OUT   <- file.path(HERE, "../results/gtex_brain_pathway_scores.tsv")

read_gmt <- function(path) {
  lines <- readLines(path)
  sets <- list()
  for (ln in lines) {
    f <- strsplit(ln, "\t")[[1]]
    if (length(f) > 2) sets[[f[1]]] <- f[3:length(f)][nzchar(f[3:length(f)])]
  }
  sets
}

cat("Querying recount3 for GTEx BRAIN ...\n")
ap <- available_projects(organism = "human")
sel <- ap[ap$project == "BRAIN" & ap$file_source == "gtex", ]
rse <- create_rse(sel)
cat("  RSE:", nrow(rse), "genes x", ncol(rse), "samples\n")

# coverage-scale raw counts -> TPM-like, then log2
assays(rse)$counts <- transform_counts(rse)
expr <- log2(assays(rse)$counts + 1)

# gene symbols (recount3 rowData has gene_name)
gn <- rowData(rse)$gene_name
rownames(expr) <- make.unique(as.character(gn))

# brain subregion label (GTEx smtsd = specific tissue detail)
cd <- as.data.frame(colData(rse))
region_col <- grep("smtsd", tolower(names(cd)), value = TRUE)[1]
if (is.na(region_col)) region_col <- grep("smts", tolower(names(cd)), value = TRUE)[1]
region <- cd[[names(cd)[tolower(names(cd)) == region_col][1]]]
cat("  region field:", region_col, "| unique regions:", length(unique(region)), "\n")

# z-score each gene across samples, then Hallmark pathway means
z <- t(scale(t(expr)))                       # gene-wise z across samples
z[is.na(z)] <- 0
pw <- read_gmt(PANEL)
scores <- sapply(names(pw), function(n) {
  g <- intersect(pw[[n]], rownames(z))
  if (length(g) < 2) return(rep(NA, ncol(z)))
  colMeans(z[g, , drop = FALSE])
})
scores <- scores[, colSums(is.na(scores)) == 0, drop = FALSE]
out <- data.frame(sample = colnames(z), region = region, scores, check.names = FALSE)
out <- out[!is.na(out$region), ]
write.table(out, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", OUT, ":", nrow(out), "samples x", ncol(scores), "pathways\n")
cat("Region counts:\n"); print(sort(table(out$region), decreasing = TRUE))
