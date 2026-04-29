suppressPackageStartupMessages({
  library(DEXSeq)
  library(BiocParallel)
  library(ggplot2)
  library(ggrepel)
  library(rlang)
})

DEX <- Sys.getenv("DEX")
stopifnot(nchar(DEX) > 0)

counts_dir <- file.path(DEX, "counts_fmt3")
flat_gff   <- file.path(DEX, "annotation", "GRCh38.dexseq.gff")
sample_tsv <- file.path(DEX, "sample_table.tsv")

stopifnot(dir.exists(counts_dir))
stopifnot(file.exists(flat_gff))
stopifnot(file.exists(sample_tsv))

outroot <- file.path(DEX, "results", "DEXSeq_4comparisons")
dir.create(outroot, showWarnings = FALSE, recursive = TRUE)

padj_cut   <- 0.05
gene_q_cut <- 0.05

ncores <- as.integer(Sys.getenv("NCORES", "8"))
if (is.na(ncores) || ncores < 1) ncores <- 1
BPP <- if (ncores > 1) BiocParallel::MulticoreParam(workers = ncores) else BiocParallel::SerialParam()

force <- as.integer(Sys.getenv("FORCE","0")) == 1

message("DEX = ", DEX)
message("Using ncores = ", ncores)
message("FORCE = ", force)

# -------------------------
# read sample table
# -------------------------
st <- read.delim(sample_tsv, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(st)[1:3] <- c("sample","condition","bam")

# Base grouping (robust): prefer sample-name prefix; fallback to provided condition
group_base <- st$condition
group_base[grepl("^C",  st$sample)] <- "C"
group_base[grepl("^E",  st$sample)] <- "E"
group_base[grepl("^NE", st$sample)] <- "NE"

is_NE685 <- grepl("^NE685", st$sample) | st$condition == "NE685"

st$group_base <- group_base
st$is_NE685   <- is_NE685

message("Condition levels seen (raw): ", paste(sort(unique(st$condition)), collapse=", "))
message("Group_base levels used: ", paste(sort(unique(st$group_base)), collapse=", "))
message("NE685 samples matched: ", sum(st$is_NE685))

# -------------------------
# helper: PCA on exon usage
# -------------------------
make_exon_usage_pca <- function(st_sub, countFiles_sub, flat_gff, outdir, tag) {

  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

  sampleTable <- data.frame(
    row.names = st_sub$sample,
    condition = factor(st_sub$condition)
  )

  dxd <- DEXSeqDataSetFromHTSeq(
    countfiles    = countFiles_sub,
    sampleData    = sampleTable,
    design        = ~ sample + exon + condition:exon,
    flattenedfile = flat_gff
  )

  dxd <- estimateSizeFactors(dxd)

  # Use normalized counts directly; then select only exon=="this"
  cn <- counts(dxd, normalized=TRUE)
  this_cols <- colData(dxd)$exon == "this"
  nc <- cn[, this_cols, drop=FALSE]
  colnames(nc) <- as.character(colData(dxd)$sample[this_cols])

  gene <- as.character(rowRanges(dxd)$groupID)
  gene_sum <- rowsum(nc, group=gene, reorder=FALSE)
  idx <- match(gene, rownames(gene_sum))
  usage <- nc / gene_sum[idx, ]
  usage[!is.finite(usage)] <- NA

  mean_nc <- rowMeans(nc, na.rm=TRUE)
  sd_use  <- apply(usage, 1, sd, na.rm=TRUE)

  keep <- mean_nc >= 10 & sd_use >= 0.02
  usage_f <- usage[keep, , drop=FALSE]

  eps <- 1e-5
  pp <- pmin(pmax(usage_f, eps), 1 - eps)
  logit <- log(pp / (1 - pp))
  logit[!is.finite(logit)] <- NA
  logit <- logit[complete.cases(logit), , drop=FALSE]

  if (nrow(logit) < 10) {
    warning("Too few exon bins after filtering for PCA in ", tag, " (n=", nrow(logit), "). Skipping PCA plots.")
    return(invisible(list(dxd=dxd, pca=NULL, kept_bins=nrow(logit))))
  }

  pca <- prcomp(t(logit), center=TRUE, scale.=TRUE)
  pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100

  pcs <- pca$x[, 1:min(5, ncol(pca$x)), drop=FALSE]
  df <- data.frame(sample = rownames(pcs), pcs, stringsAsFactors=FALSE)
  df <- merge(df, st_sub[,c("sample","condition")], by="sample", all.x=TRUE)

  write.table(df, file.path(outdir, paste0("PCA_coordinates_", tag, ".tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)

  plot_pair <- function(xpc, ypc) {
    ggplot(df, aes(x = .data[[xpc]], y = .data[[ypc]], color=condition, label=sample)) +
      geom_point(size=3) +
      ggrepel::geom_text_repel(max.overlaps=100, size=3) +
      theme_bw() +
      labs(
        title = paste0("PCA exon-usage: ", tag),
        x = sprintf("%s (%.1f%%)", xpc, pct[as.integer(sub("PC","",xpc))]),
        y = sprintf("%s (%.1f%%)", ypc, pct[as.integer(sub("PC","",ypc))])
      )
  }

  p12 <- plot_pair("PC1","PC2")
  ggsave(file.path(outdir, paste0("PCA_", tag, "_PC1_PC2.pdf")), p12, width=9, height=7)

  if (ncol(pca$x) >= 3) {
    p13 <- plot_pair("PC1","PC3")
    p23 <- plot_pair("PC2","PC3")
    ggsave(file.path(outdir, paste0("PCA_", tag, "_PC1_PC3.pdf")), p13, width=9, height=7)
    ggsave(file.path(outdir, paste0("PCA_", tag, "_PC2_PC3.pdf")), p23, width=9, height=7)
  }

  invisible(list(dxd=dxd, pca=pca, kept_bins=nrow(logit)))
}

# -------------------------
# helper: flatten DEXSeqResults to plain df for TSV
# -------------------------
flatten_dxr <- function(dxr) {
  df <- as.data.frame(dxr)
  df$exon_bin_id <- rownames(dxr)

  # transcripts: list -> string
  if ("transcripts" %in% names(df)) {
    tr <- dxr$transcripts
    df$transcripts <- vapply(tr, function(x) paste(x, collapse=";"), character(1))
  }

  # genomicData: GRanges -> columns
  if ("genomicData" %in% names(df)) {
    gd <- dxr$genomicData
    df$genomicData.seqnames <- as.character(GenomicRanges::seqnames(gd))
    df$genomicData.start    <- GenomicRanges::start(gd)
    df$genomicData.end      <- GenomicRanges::end(gd)
    df$genomicData.width    <- GenomicRanges::width(gd)
    df$genomicData.strand   <- as.character(GenomicRanges::strand(gd))
    df$genomicData <- NULL
  }

  # countData: matrix -> many numeric columns
  if ("countData" %in% names(df)) {
    cd <- dxr$countData
    cd_df <- as.data.frame(cd)
    colnames(cd_df) <- paste0("countData.", colnames(cd))
    df <- cbind(df, cd_df)
    df$countData <- NULL
  }

  # Any remaining list columns -> stringify (last-resort)
  is_list_col <- vapply(df, is.list, logical(1))
  if (any(is_list_col)) {
    for (nm in names(df)[is_list_col]) {
      df[[nm]] <- vapply(df[[nm]], function(x) paste(x, collapse=";"), character(1))
    }
  }

  df
}

# -------------------------
# helper: run DEXSeq + write TSVs
# -------------------------
run_dexseq <- function(st_sub, countFiles_sub, flat_gff, outdir, tag, BPPARAM, padj_cut, gene_q_cut) {

  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

  sampleTable <- data.frame(
    row.names = st_sub$sample,
    condition = factor(st_sub$condition)
  )

  dxd <- DEXSeqDataSetFromHTSeq(
    countfiles    = countFiles_sub,
    sampleData    = sampleTable,
    design        = ~ sample + exon + condition:exon,
    flattenedfile = flat_gff
  )

  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM)
  dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

  dxr <- DEXSeqResults(dxd)

  exons_df <- flatten_dxr(dxr)

  exons_all_path <- file.path(outdir, paste0("exons_all_", tag, ".tsv"))
  write.table(exons_df, exons_all_path, sep="\t", quote=FALSE, row.names=FALSE)

  exons_sig <- exons_df[!is.na(exons_df$padj) & exons_df$padj < padj_cut, , drop=FALSE]
  exons_sig_path <- file.path(outdir, paste0("exons_sig_padj", padj_cut, "_", tag, ".tsv"))
  write.table(exons_sig, exons_sig_path, sep="\t", quote=FALSE, row.names=FALSE)

  # gene-level q-values
  q <- perGeneQValue(dxr)
  genes_df <- data.frame(
    groupID = names(q),
    perGeneQValue = as.numeric(q),
    stringsAsFactors = FALSE
  )

  n_tested <- as.integer(table(exons_df$groupID)[genes_df$groupID])
  n_tested[is.na(n_tested)] <- 0
  genes_df$n_exons_tested <- n_tested

  n_sig <- as.integer(table(exons_sig$groupID)[genes_df$groupID])
  n_sig[is.na(n_sig)] <- 0
  genes_df$n_exons_sig_padj <- n_sig

  minpadj <- tapply(exons_df$padj, exons_df$groupID, function(x) suppressWarnings(min(x, na.rm=TRUE)))
  genes_df$min_exon_padj <- as.numeric(minpadj[genes_df$groupID])

  genes_all_path <- file.path(outdir, paste0("genes_qvalues_", tag, ".tsv"))
  write.table(genes_df, genes_all_path, sep="\t", quote=FALSE, row.names=FALSE)

  genes_sig <- genes_df[!is.na(genes_df$perGeneQValue) & genes_df$perGeneQValue < gene_q_cut, , drop=FALSE]
  genes_sig_path <- file.path(outdir, paste0("genes_sig_q", gene_q_cut, "_", tag, ".tsv"))
  write.table(genes_sig, genes_sig_path, sep="\t", quote=FALSE, row.names=FALSE)

  summary_df <- data.frame(
    tag = tag,
    n_samples = nrow(st_sub),
    n_groups = length(unique(st_sub$condition)),
    n_exon_bins = nrow(exons_df),
    n_exons_sig_padj = nrow(exons_sig),
    n_genes = nrow(genes_df),
    n_genes_sig_q = nrow(genes_sig),
    stringsAsFactors = FALSE
  )
  write.table(summary_df, file.path(outdir, paste0("summary_", tag, ".tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)

  pdf(file.path(outdir, paste0("dispersion_plot_", tag, ".pdf")))
  plotDispEsts(dxd)
  dev.off()

  saveRDS(dxd, file.path(outdir, paste0("dxd_", tag, ".rds")))
  saveRDS(dxr, file.path(outdir, paste0("dxr_", tag, ".rds")))

  invisible(list(dxd=dxd, dxr=dxr))
}

# -------------------------
# define 4 comparisons (your exact intent)
# -------------------------
comparisons <- list()

# C vs E
s1 <- st[st$group_base %in% c("C","E"), , drop=FALSE]
s1$condition <- s1$group_base
comparisons[[length(comparisons)+1]] <- list(tag="C_vs_E", st_sub=s1)

# C vs NE
s2 <- st[st$group_base %in% c("C","NE"), , drop=FALSE]
s2$condition <- s2$group_base
comparisons[[length(comparisons)+1]] <- list(tag="C_vs_NE", st_sub=s2)

# E vs NE
s3 <- st[st$group_base %in% c("E","NE"), , drop=FALSE]
s3$condition <- s3$group_base
comparisons[[length(comparisons)+1]] <- list(tag="E_vs_NE", st_sub=s3)

# C vs (all E + the specific NE685 sample)  => label as "EplusNE685"
s4 <- st[(st$group_base %in% c("C","E")) | st$is_NE685, , drop=FALSE]
s4$condition <- ifelse(s4$group_base == "C", "C", "EplusNE685")
s4$condition[s4$is_NE685] <- "EplusNE685"
comparisons[[length(comparisons)+1]] <- list(tag="C_vs_Eplus_NE685", st_sub=s4)

# -------------------------
# run each comparison
# -------------------------
for (cmp in comparisons) {
  tag <- cmp$tag
  st_sub <- cmp$st_sub
  st_sub$condition <- factor(st_sub$condition)

  message("\n=== Running: ", tag, " ===")
  message("Samples: ", nrow(st_sub), " | Groups: ", paste(levels(st_sub$condition), collapse=", "))

  out_cmp <- file.path(outroot, tag)
  out_pca <- file.path(out_cmp, "PCA")
  out_dx  <- file.path(out_cmp, "DEXSeq")

  summary_path <- file.path(out_dx, paste0("summary_", tag, ".tsv"))
  if (!force && file.exists(summary_path)) {
    message("Skipping (already completed): ", tag)
    next
  }

  countFiles_sub <- file.path(counts_dir, paste0(st_sub$sample, ".dexseq.txt"))
  if (!all(file.exists(countFiles_sub))) {
    missing <- countFiles_sub[!file.exists(countFiles_sub)]
    stop("Missing count files (first few):\n", paste(head(missing, 10), collapse="\n"))
  }

  pca_res <- make_exon_usage_pca(st_sub, countFiles_sub, flat_gff, out_pca, tag)
  message("PCA kept exon bins: ", pca_res$kept_bins)

  run_dexseq(st_sub, countFiles_sub, flat_gff, out_dx, tag, BPP, padj_cut, gene_q_cut)

  message("Done: ", tag)
}

message("\nAll comparisons finished. Outputs in:\n", outroot)