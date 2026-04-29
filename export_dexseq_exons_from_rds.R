#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DEXSeq)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 4) {
  stop(
    "Usage: Rscript export_dexseq_exons_from_rds.R <dxr_rds> <outdir> <tag> [padj_cut]\n",
    "Example:\n",
    "  Rscript export_dexseq_exons_from_rds.R dxr_PPCD_vs_Control_Select_withBatch.rds ./ PPCD_vs_Control_Select_withBatch 0.05"
  )
}

dxr_rds  <- args[1]
outdir   <- args[2]
tag      <- args[3]
padj_cut <- if (length(args) >= 4) as.numeric(args[4]) else 0.05

if (!file.exists(dxr_rds)) stop("Missing dxr RDS: ", dxr_rds)
if (is.na(padj_cut)) stop("padj_cut must be numeric")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

flatten_dxr_exact <- function(dxr) {
  df <- as.data.frame(dxr, stringsAsFactors = FALSE, check.names = FALSE)
  df$exon_bin_id <- rownames(dxr)

  if ("transcripts" %in% names(df)) {
    tr <- dxr$transcripts
    df$transcripts <- vapply(tr, function(x) paste(x, collapse = ";"), character(1))
    rm(tr)
  }

  if ("genomicData" %in% names(df)) {
    gd <- dxr$genomicData
    df$genomicData.seqnames <- as.character(GenomicRanges::seqnames(gd))
    df$genomicData.start <- GenomicRanges::start(gd)
    df$genomicData.end <- GenomicRanges::end(gd)
    df$genomicData.width <- GenomicRanges::width(gd)
    df$genomicData.strand <- as.character(GenomicRanges::strand(gd))
    df$genomicData <- NULL
    rm(gd)
  }

  if ("countData" %in% names(df)) {
    cd <- dxr$countData
    cd_df <- as.data.frame(cd, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(cd_df) <- paste0("countData.", colnames(cd_df))
    df <- cbind(df, cd_df)
    df$countData <- NULL
    rm(cd, cd_df)
  }

  is_list_col <- vapply(df, is.list, logical(1))
  if (any(is_list_col)) {
    for (nm in names(df)[is_list_col]) {
      df[[nm]] <- vapply(df[[nm]], function(x) paste(x, collapse = ";"), character(1))
    }
  }

  df
}

message("Reading dxr from: ", dxr_rds)
dxr <- readRDS(dxr_rds)
message("Loaded dxr with ", nrow(dxr), " exon bins")

message("Flattening full exon table...")
exons_df <- flatten_dxr_exact(dxr)
rm(dxr)
gc()

outfile_all <- file.path(outdir, paste0("exons_all_", tag, ".tsv"))
message("Writing full exon table: ", outfile_all)
fwrite(
  as.data.table(exons_df),
  outfile_all,
  sep = "\t",
  quote = FALSE
)

message("Subsetting significant exon bins at padj < ", padj_cut, " ...")
exons_sig <- exons_df[!is.na(exons_df$padj) & exons_df$padj < padj_cut, , drop = FALSE]

outfile_sig <- file.path(outdir, paste0("exons_sig_padj", padj_cut, "_", tag, ".tsv"))
message("Writing significant exon table: ", outfile_sig)
fwrite(
  as.data.table(exons_sig),
  outfile_sig,
  sep = "\t",
  quote = FALSE
)

message("Significant exon bins written: ", nrow(exons_sig))

rm(exons_df, exons_sig)
gc()

message("Done.")