cat > $DEX/run_dexseq_4_comparisons.R << 'RSCRIPT'
suppressPackageStartupMessages({
  library(DEXSeq)
})

base <- "/mnt/scratch/hqurashi/12_Nihar_DEXSeq/dexseq"
samples <- read.delim(file.path(base, "sample_table.tsv"), stringsAsFactors = FALSE)
samples$countFile <- file.path(base, "counts", paste0(samples$sample, ".dexseq.txt"))
flattened <- file.path(base, "annotation", "GRCh38.dexseq.gff")
out_base <- file.path(base, "results")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

stopifnot(all(file.exists(samples$countFile)))
stopifnot(file.exists(flattened))

run_one <- function(name, samples_sub, condition_col = "condition") {
  message("=== ", name, " ===")
  outdir <- file.path(out_base, name)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  samples_sub[[condition_col]] <- factor(samples_sub[[condition_col]])
  samples_sub$sample <- factor(samples_sub$sample)

  dxd <- DEXSeqDataSetFromHTSeq(
    countfiles    = samples_sub$countFile,
    sampleData    = samples_sub,
    design        = as.formula(paste0("~ sample + exon + ", condition_col, ":exon")),
    flattenedfile = flattened
  )

  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = condition_col)
  dxr <- DEXSeqResults(dxd)

  # Exon-bin results
  res <- as.data.frame(dxr)
  res$featureID <- dxr$featureID
  res$groupID <- dxr$groupID
  write.table(res, file.path(outdir, "dexseq_exon_results.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  # Gene-level q-values
  gq <- perGeneQValue(dxr)
  gq_df <- data.frame(groupID = names(gq), perGeneQValue = as.numeric(gq))
  gq_df <- gq_df[order(gq_df$perGeneQValue), ]
  write.table(gq_df, file.path(outdir, "dexseq_gene_qvalues.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  # Convenience: significant genes at q<0.05
  sig <- subset(gq_df, perGeneQValue < 0.05)
  write.table(sig, file.path(outdir, "dexseq_sig_genes_qlt0.05.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  # --- Splicing-based clustering (PCA) using exon-usage proportions ---
  cts <- counts(dxd, normalized = TRUE)       # exonBins x samples
  gene <- dxr$groupID
  gene_totals <- rowsum(cts, group = gene)   # genes x samples
  denom <- gene_totals[gene, , drop=FALSE]   # exonBins x samples

  usage <- cts / pmax(denom, 1e-6)           # proportion within gene
  keep <- rowSums(cts) >= 50                 # filter low info bins
  usage <- usage[keep, , drop=FALSE]

  # logit transform for PCA stability
  eps <- 1e-6
  u <- pmin(pmax(usage, eps), 1 - eps)
  u_tr <- log2(u / (1 - u))

  vars <- apply(u_tr, 1, var)
  top <- order(vars, decreasing = TRUE)[seq_len(min(2000, length(vars)))]
  mat <- t(u_tr[top, , drop=FALSE])

  pca <- prcomp(mat, scale. = TRUE)

  pca_df <- data.frame(sample = rownames(pca$x), pca$x[,1:5, drop=FALSE])
  meta <- as.data.frame(colData(dxd))
  meta$sample <- rownames(meta)
  pca_df <- merge(pca_df, meta, by="sample", all.x=TRUE)
  write.table(pca_df, file.path(outdir, "splicing_PCA_coordinates.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  pdf(file.path(outdir, "splicing_PCA_PC1_PC2.pdf"), width=7, height=6)
  plot(pca$x[,1], pca$x[,2], pch=19,
       xlab="PC1 (exon-usage)", ylab="PC2 (exon-usage)")
  text(pca$x[,1], pca$x[,2], labels=rownames(pca$x), pos=3, cex=0.7)
  dev.off()

  saveRDS(dxd, file.path(outdir, "dxd.rds"))
  saveRDS(dxr, file.path(outdir, "dxr.rds"))
}

# 1) C vs E
run_one("C_vs_E", subset(samples, condition %in% c("C","E")))

# 2) C vs NE
run_one("C_vs_NE", subset(samples, condition %in% c("C","NE")))

# 3) E vs NE
run_one("E_vs_NE", subset(samples, condition %in% c("E","NE")))

# 4) C vs (E + NE685)
samples4 <- samples
samples4$condition4 <- NA
samples4$condition4[samples4$condition == "C"] <- "C"
samples4$condition4[samples4$condition == "E"] <- "CASE"
samples4$condition4[grepl("^NE685_", samples4$sample)] <- "CASE"
samples4 <- subset(samples4, !is.na(condition4))
run_one("C_vs_Eplus_NE685", samples4, condition_col="condition4")
RSCRIPT