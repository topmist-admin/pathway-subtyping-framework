#!/usr/bin/env Rscript
# Flag which candidate SRA studies are in recount3 (the uniformly-processed
# Track-A set). Reads results/candidate_studies.tsv (from enumerate_candidates.py),
# checks each SRP against recount3 available_projects, writes the Track-A table.
suppressMessages(library(recount3))

HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(HERE) == 0 || HERE == "") HERE <- "."
IN  <- file.path(HERE, "../results/candidate_studies.tsv")
OUT <- file.path(HERE, "../results/track_a_recount3.tsv")

cand <- read.delim(IN, stringsAsFactors = FALSE)
cand <- cand[cand$keep == "True" & nzchar(cand$srp), ]
cat("candidates with SRP:", nrow(cand), "\n")

ap <- available_projects(organism = "human")
sra <- ap[ap$file_source == "sra", c("project", "n_samples")]
m <- merge(cand, sra, by.x = "srp", by.y = "project", all.x = TRUE)
m$in_recount3 <- !is.na(m$n_samples.y)
m <- m[order(-m$in_recount3, -m$n_samples.y), ]

write.table(m[, c("gse", "srp", "n_samples.x", "n_samples.y", "in_recount3")],
            OUT, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("gse", "srp", "n_geo", "n_recount3", "in_recount3"))

inr <- m[m$in_recount3, ]
cat("IN recount3 (Track-A uniform):", nrow(inr), "studies,",
    sum(inr$n_samples.y, na.rm = TRUE), "recount3 samples\n")
for (i in seq_len(nrow(inr)))
  cat("  ", inr$gse[i], inr$srp[i], "n_recount3 =", inr$n_samples.y[i], "\n")
cat("\nWrote", OUT, "\n")
