#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# PPCD vs Control SELECT DEXSeq all-in-one workflow
#
# Assumes:
#   - conda env "dexseq_grch38" is already activated
#   - BAMs are in nf-core/rnaseq star_salmon output directory
#   - you want SELECT samples only:
#       C23 C38 PPCD1S1 PPCD1S2 PPCD1S3 PPCD1S4 C18 C56 C59 PPCD1S5P0
#   - no batch correction in the DEXSeq model
#   - dexseq_prepare_annotation.py should use -r no
#
# Outputs:
#   /media/pontikos_nas2/HarisQurashi/projects/12_Nihar_DEXSeq/PPCD_vs_Control_Select/outputs/DEXSeq
###############################################################################

############################
# USER CONFIG
############################
STAR_SALMON_DIR="/media/pontikos_nas2/HarisQurashi/projects/12_Nihar_DEXSeq/PPCD_vs_Control/outputs/star_salmon"
GTF_GZ="/media/pontikos_nas2/HarisQurashi/refs/Homo_sapiens.GRCh38.115.gtf.gz"
OUT_BASE="/mnt/scratch/hqurashi/12_Nihar_DEXSeq_PPCD_vs_Control/dexseq_NEW"

STRAND="reverse"         # reverse / yes / no
COUNT_THREADS=6          # for dexseq_count.py via GNU parallel
FORCE_ANNOTATION=0       # 1 = remake local GTF + flattened GFF
FORCE_COUNTS=0           # 1 = rerun exon-bin counting for selected samples
FORCE_CLEAN=0            # 1 = recreate counts_clean in R step
FORCE_RESULTS=0          # 1 = wipe DEXSeq results dir before R analysis

PADJ_CUT="0.05"
GENE_Q_CUT="0.05"

############################
# CHECKS
############################
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found. Activate dexseq_grch38."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found. Activate dexseq_grch38."; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "GNU parallel not found. Activate dexseq_grch38."; exit 1; }
DEXSEQ_PREPARE="$(Rscript -e 'cat(system.file("python_scripts","dexseq_prepare_annotation.py",package="DEXSeq"))')"
[[ -n "$DEXSEQ_PREPARE" && -f "$DEXSEQ_PREPARE" ]] || {
  echo "Could not find dexseq_prepare_annotation.py via R DEXSeq. Got: '$DEXSEQ_PREPARE'"
  exit 1
}

[[ -d "$STAR_SALMON_DIR" ]] || { echo "Missing STAR_SALMON_DIR: $STAR_SALMON_DIR"; exit 1; }
[[ -f "$GTF_GZ" ]] || { echo "Missing GTF_GZ: $GTF_GZ"; exit 1; }

DEXROOT="$OUT_BASE"
ANN_DIR="$DEXROOT/annotation"
MERGED_DIR="$DEXROOT/merged_bams"
COUNTS_DIR="$DEXROOT/counts"
COUNTS_CLEAN_DIR="$DEXROOT/counts_clean"
LOG_DIR="$DEXROOT/logs"
COUNT_LOG_DIR="$LOG_DIR/count_logs"
RESULTS_ROOT="$DEXROOT/results/DEXSeq_PPCD_vs_Control_Select"
SCRIPT_DIR="$DEXROOT/scripts_generated"

mkdir -p "$ANN_DIR" "$MERGED_DIR" "$COUNTS_DIR" "$COUNTS_CLEAN_DIR" \
         "$LOG_DIR" "$COUNT_LOG_DIR" "$RESULTS_ROOT" "$SCRIPT_DIR"

echo "STAR_SALMON_DIR = $STAR_SALMON_DIR"
echo "GTF_GZ          = $GTF_GZ"
echo "OUT_BASE        = $OUT_BASE"
echo "STRAND          = $STRAND"
echo "COUNT_THREADS   = $COUNT_THREADS"
echo "FORCE_ANNOTATION= $FORCE_ANNOTATION"
echo "FORCE_COUNTS    = $FORCE_COUNTS"
echo "FORCE_CLEAN     = $FORCE_CLEAN"
echo "FORCE_RESULTS   = $FORCE_RESULTS"
echo

############################
# SAMPLE DEFINITIONS (SELECT SET ONLY)
############################
samples=(
  C23 C38 PPCD1S1 PPCD1S2 PPCD1S3 PPCD1S4
  C18 C56 C59 PPCD1S5P0
)

declare -A batch=(
  [C23]=B1 [C38]=B1 [PPCD1S1]=B1 [PPCD1S2]=B1 [PPCD1S3]=B1 [PPCD1S4]=B1
  [C18]=B2 [C56]=B2 [C59]=B2 [PPCD1S5P0]=B2
)

declare -A cond=(
  [C23]=Control [C38]=Control [C18]=Control [C56]=Control [C59]=Control
  [PPCD1S1]=PPCD [PPCD1S2]=PPCD [PPCD1S3]=PPCD [PPCD1S4]=PPCD [PPCD1S5P0]=PPCD
)

############################
# STEP 1: PREPARE LOCAL GTF + DEXSEQ FLATTENED GFF (-r no)
############################
GTF_LOCAL="$ANN_DIR/Homo_sapiens.GRCh38.115.gtf"
FLAT_GFF="$ANN_DIR/GRCh38.115.dexseq.r_no.gff"

if [[ "$FORCE_ANNOTATION" == "1" ]]; then
  rm -f "$GTF_LOCAL" "$FLAT_GFF"
fi

if [[ ! -f "$GTF_LOCAL" ]]; then
  echo "[1/5] Decompressing GTF to local annotation directory..."
  gzip -dc "$GTF_GZ" > "$GTF_LOCAL"
else
  echo "[1/5] Reusing local GTF: $GTF_LOCAL"
fi

echo -n "Local GTF exon count: "
awk '$3=="exon"{c++} END{print c+0}' "$GTF_LOCAL"

if [[ ! -f "$FLAT_GFF" ]]; then
  echo "[1/5] Creating flattened DEXSeq GFF with -r no ..."
  python "$DEXSEQ_PREPARE" -r no "$GTF_LOCAL" "$FLAT_GFF"
else
  echo "[1/5] Reusing flattened GFF: $FLAT_GFF"
fi

############################
# STEP 2: BUILD SAMPLE TABLE + LINK SELECT BAMs
############################
SAMPLE_TSV="$DEXROOT/sample_table.tsv"

echo "[2/5] Building sample table for SELECT samples ..."
{
  printf "sample\tcondition\tbatch\tbam\n"
  for s in "${samples[@]}"; do
    [[ -n "${batch[$s]:-}" ]] || { echo "Missing batch for $s"; exit 1; }
    [[ -n "${cond[$s]:-}"  ]] || { echo "Missing condition for $s"; exit 1; }

    bam="$STAR_SALMON_DIR/${s}.markdup.sorted.bam"
    if [[ ! -f "$bam" ]]; then
      shopt -s nullglob
      matches=("$STAR_SALMON_DIR/${s}"*.markdup.sorted.bam)
      shopt -u nullglob
      if (( ${#matches[@]} == 1 )); then
        bam="${matches[0]}"
      else
        echo "Could not uniquely find BAM for sample $s in $STAR_SALMON_DIR" >&2
        echo "Tried exact and glob match." >&2
        exit 1
      fi
    fi

    link="$MERGED_DIR/${s}.markdup.sorted.bam"
    ln -sf "$bam" "$link"

    if [[ -f "${bam}.bai" ]]; then
      ln -sf "${bam}.bai" "${link}.bai"
    elif [[ -f "${bam%.bam}.bai" ]]; then
      ln -sf "${bam%.bam}.bai" "${link}.bai"
    else
      samtools index -@ 4 "$link"
    fi

    printf "%s\t%s\t%s\t%s\n" "$s" "${cond[$s]}" "${batch[$s]}" "$link"
  done
} > "$SAMPLE_TSV"

echo "Wrote: $SAMPLE_TSV"
if command -v column >/dev/null 2>&1; then
  column -t "$SAMPLE_TSV"
else
  cat "$SAMPLE_TSV"
fi
echo

############################
# STEP 3: GENERATE EXON-BIN COUNT FILES
############################
DEXSEQ_COUNT="$(Rscript -e 'cat(system.file("python_scripts","dexseq_count.py",package="DEXSeq"))')"
[[ -n "$DEXSEQ_COUNT" && -f "$DEXSEQ_COUNT" ]] || {
  echo "Could not find dexseq_count.py via R DEXSeq. Got: '$DEXSEQ_COUNT'"
  exit 1
}

echo "[3/5] Using dexseq_count.py: $DEXSEQ_COUNT"

COUNT_JOBS_TSV="$LOG_DIR/count_jobs.tsv"
awk -v DEXDIR="$DEXROOT" 'BEGIN{OFS="\t"} NR>1 { print $1, $4, DEXDIR"/counts/"$1".dexseq.txt" }' \
  "$SAMPLE_TSV" > "$COUNT_JOBS_TSV"

echo "Wrote count job list: $COUNT_JOBS_TSV"
head -n 5 "$COUNT_JOBS_TSV" || true
echo

if [[ "$FORCE_COUNTS" == "1" ]]; then
  echo "[3/5] FORCE_COUNTS=1, removing existing selected count files ..."
  while IFS=$'\t' read -r sample bam out; do
    rm -f "$out"
  done < "$COUNT_JOBS_TSV"
fi

export DEXSEQ_COUNT FLAT_GFF STRAND DEXROOT COUNT_LOG_DIR

parallel -j "$COUNT_THREADS" --colsep '\t' '
  set -euo pipefail
  sample={1}
  bam={2}
  out={3}

  if [[ -s "$out" ]]; then
    echo "[skip] ${sample} (counts already present)"
    exit 0
  fi

  echo "[count] ${sample}"
  python "$DEXSEQ_COUNT" \
    -f bam \
    -p yes \
    -r pos \
    -s "$STRAND" \
    "$FLAT_GFF" \
    "$bam" \
    "$out" \
    > "$COUNT_LOG_DIR/${sample}.log" 2>&1
' :::: "$COUNT_JOBS_TSV"

echo
echo "[3/5] Counts written to: $COUNTS_DIR"
ls -lh "$COUNTS_DIR" | head
echo

############################
# STEP 4: WRITE R SCRIPT
############################
R_SCRIPT="$SCRIPT_DIR/run_DEXSeq_PPCD_vs_Control_Select_noBatch.R"

cat > "$R_SCRIPT" << 'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DEXSeq)
  library(BiocParallel)
  library(ggplot2)
  library(ggrepel)
  library(rlang)
  library(data.table)
})

# =========================
# CONFIG
# =========================
DEXROOT <- Sys.getenv("DEXROOT")
if (nchar(DEXROOT) == 0) stop("DEXROOT env var not set")

sample_tsv <- file.path(DEXROOT, "sample_table.tsv")
counts_dir_raw <- file.path(DEXROOT, "counts")
counts_dir_clean <- file.path(DEXROOT, "counts_clean")
flat_gff <- file.path(DEXROOT, "annotation", "GRCh38.115.dexseq.r_no.gff")
outroot <- file.path(DEXROOT, "results", "DEXSeq_PPCD_vs_Control_Select")

dir.create(counts_dir_clean, recursive = TRUE, showWarnings = FALSE)
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)

padj_cut   <- as.numeric(Sys.getenv("PADJ_CUT", "0.05"))
gene_q_cut <- as.numeric(Sys.getenv("GENE_Q_CUT", "0.05"))

force_clean   <- as.integer(Sys.getenv("FORCE_CLEAN", "0")) == 1
force_results <- as.integer(Sys.getenv("FORCE_RESULTS", "0")) == 1

message("DEXROOT          = ", DEXROOT)
message("sample_tsv       = ", sample_tsv)
message("counts_dir_raw   = ", counts_dir_raw)
message("counts_dir_clean = ", counts_dir_clean)
message("flat_gff         = ", flat_gff)
message("outroot          = ", outroot)
message("PADJ_CUT         = ", padj_cut)
message("GENE_Q_CUT       = ", gene_q_cut)
message("FORCE_CLEAN      = ", force_clean)
message("FORCE_RESULTS    = ", force_results)

stopifnot(file.exists(sample_tsv))
stopifnot(dir.exists(counts_dir_raw))
stopifnot(file.exists(flat_gff))

st <- read.delim(sample_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("sample", "condition", "batch", "bam") %in% colnames(st)))

st$condition <- factor(st$condition, levels = c("Control", "PPCD"))
st$batch <- factor(st$batch, levels = unique(st$batch))

message("Samples: ", paste(st$sample, collapse = ", "))
message("Condition x batch:\n", paste(capture.output(print(table(st$condition, st$batch))), collapse = "\n"))

# =========================
# Clean DEXSeq count files into STRICT 2-col format
# =========================
clean_one_countfile <- function(infile, outfile) {
  dt <- data.table::fread(
    infile, header = FALSE, sep = "\t", quote = "",
    fill = TRUE, data.table = TRUE, showProgress = FALSE
  )
  if (ncol(dt) < 2) stop("Unexpected <2 columns in: ", infile)

  dt <- dt[, 1:2]
  data.table::setnames(dt, c("feature", "count_raw"))

  dt$feature <- gsub('"', "", as.character(dt$feature), fixed = TRUE)
  dt$count_raw <- as.character(dt$count_raw)

  keep <- nzchar(dt$feature) &
    !(startsWith(dt$feature, "_") | startsWith(dt$feature, "__")) &
    nzchar(dt$count_raw)
  dt <- dt[keep]

  dt <- dt[grepl(":", dt$feature, fixed = TRUE)]

  dt[, count := suppressWarnings(as.integer(count_raw))]
  dt <- dt[!is.na(count)]

  if (nrow(dt) == 0) stop("After cleaning, no usable rows left in: ", infile)

  data.table::fwrite(
    dt[, .(feature, count)],
    file = outfile,
    sep = "\t", col.names = FALSE, quote = FALSE
  )
}

raw_files <- file.path(counts_dir_raw, paste0(st$sample, ".dexseq.txt"))
clean_files <- file.path(counts_dir_clean, paste0(st$sample, ".dexseq.txt"))

missing_raw <- raw_files[!file.exists(raw_files)]
if (length(missing_raw) > 0) {
  stop("Missing raw count files:\n", paste(missing_raw, collapse = "\n"))
}

message("\nCleaning raw counts -> counts_clean ...")
for (i in seq_along(raw_files)) {
  inf <- raw_files[i]
  outf <- clean_files[i]
  if (!force_clean && file.exists(outf) && file.info(outf)$size > 0) next
  message("  ", basename(inf), " -> ", basename(outf))
  clean_one_countfile(inf, outf)
}

message("Validating counts_clean ...")
for (f in clean_files) {
  cf <- count.fields(f, sep = "\t", quote = "", blank.lines.skip = FALSE, comment.char = "")
  bad <- which(cf != 2)
  if (length(bad) > 0) {
    i <- bad[1]
    ln <- readLines(f, n = i, warn = FALSE)
    stop(
      "Bad cleaned count file: ", f,
      "\nFirst bad line: ", i,
      " (NF=", cf[i], ")\nLine: ", ln[i]
    )
  }
}
message("counts_clean OK: ", counts_dir_clean)

# =========================
# PCA helper (exon usage)
# =========================
make_exon_usage_pca <- function(st_sub, countFiles_sub, flat_gff, outdir, tag) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  sampleTable <- data.frame(
    row.names = st_sub$sample,
    condition = factor(st_sub$condition),
    batch = factor(st_sub$batch)
  )

  dxd <- DEXSeqDataSetFromHTSeq(
    countfiles = countFiles_sub,
    sampleData = sampleTable,
    design = ~ sample + exon + condition:exon,
    flattenedfile = flat_gff
  )

  dxd <- estimateSizeFactors(dxd)

  cn <- counts(dxd, normalized = TRUE)
  this_cols <- colData(dxd)$exon == "this"
  nc <- cn[, this_cols, drop = FALSE]
  colnames(nc) <- as.character(colData(dxd)$sample[this_cols])

  gene <- as.character(rowRanges(dxd)$groupID)
  gene_sum <- rowsum(nc, group = gene, reorder = FALSE)
  idx <- match(gene, rownames(gene_sum))
  usage <- nc / gene_sum[idx, ]
  usage[!is.finite(usage)] <- NA

  mean_nc <- rowMeans(nc, na.rm = TRUE)
  sd_use <- apply(usage, 1, sd, na.rm = TRUE)

  keep <- mean_nc >= 10 & sd_use >= 0.02
  usage_f <- usage[keep, , drop = FALSE]

  eps <- 1e-5
  pp <- pmin(pmax(usage_f, eps), 1 - eps)
  logit <- log(pp / (1 - pp))
  logit[!is.finite(logit)] <- NA
  logit <- logit[complete.cases(logit), , drop = FALSE]

  if (nrow(logit) < 10) {
    warning("Too few exon bins after filtering for PCA in ", tag, " (n=", nrow(logit), "). Skipping PCA.")
    return(invisible(list(dxd = dxd, pca = NULL, kept_bins = nrow(logit))))
  }

  pca <- prcomp(t(logit), center = TRUE, scale. = TRUE)
  pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100

  pcs <- pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE]
  df <- data.frame(sample = rownames(pcs), pcs, stringsAsFactors = FALSE)
  df <- merge(df, st_sub[, c("sample", "condition", "batch")], by = "sample", all.x = TRUE)

  write.table(
    df,
    file.path(outdir, paste0("PCA_coordinates_", tag, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  plot_pair <- function(xpc, ypc) {
    ggplot(df, aes(x = .data[[xpc]], y = .data[[ypc]], color = condition, label = sample, shape = batch)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(max.overlaps = 100, size = 3) +
      theme_bw() +
      labs(
        title = paste0("PCA exon-usage: ", tag),
        x = sprintf("%s (%.1f%%)", xpc, pct[as.integer(sub("PC", "", xpc))]),
        y = sprintf("%s (%.1f%%)", ypc, pct[as.integer(sub("PC", "", ypc))])
      )
  }

  ggsave(file.path(outdir, paste0("PCA_", tag, "_PC1_PC2.pdf")), plot_pair("PC1", "PC2"), width = 9, height = 7)
  if (ncol(pca$x) >= 3) {
    ggsave(file.path(outdir, paste0("PCA_", tag, "_PC1_PC3.pdf")), plot_pair("PC1", "PC3"), width = 9, height = 7)
    ggsave(file.path(outdir, paste0("PCA_", tag, "_PC2_PC3.pdf")), plot_pair("PC2", "PC3"), width = 9, height = 7)
  }

  invisible(list(dxd = dxd, pca = pca, kept_bins = nrow(logit)))
}

# =========================
# Flatten DEXSeqResults so it can be written cleanly
# =========================
flatten_dxr <- function(dxr) {
  df <- as.data.frame(dxr, stringsAsFactors = FALSE, check.names = FALSE)
  df$exon_bin_id <- rownames(dxr)

  if ("transcripts" %in% names(df)) {
    tr <- dxr$transcripts
    df$transcripts <- vapply(tr, function(x) paste(x, collapse = ";"), character(1))
  }

  if ("genomicData" %in% names(df)) {
    gd <- dxr$genomicData
    df$genomicData.seqnames <- as.character(GenomicRanges::seqnames(gd))
    df$genomicData.start <- GenomicRanges::start(gd)
    df$genomicData.end <- GenomicRanges::end(gd)
    df$genomicData.width <- GenomicRanges::width(gd)
    df$genomicData.strand <- as.character(GenomicRanges::strand(gd))
    df$genomicData <- NULL
  }

  if ("countData" %in% names(df)) {
    cd <- dxr$countData
    cd_df <- as.data.frame(cd, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(cd_df) <- paste0("countData.", colnames(cd_df))
    df <- cbind(df, cd_df)
    df$countData <- NULL
  }

  is_list_col <- vapply(df, is.list, logical(1))
  if (any(is_list_col)) {
    for (nm in names(df)[is_list_col]) {
      df[[nm]] <- vapply(df[[nm]], function(x) paste(x, collapse = ";"), character(1))
    }
  }

  df
}

# =========================
# Run DEXSeq (NO batch correction)
# =========================
run_dexseq <- function(st_sub, countFiles_sub, flat_gff, outdir, tag, design_formula, BPPARAM, padj_cut, gene_q_cut) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  sampleTable <- data.frame(
    row.names = st_sub$sample,
    condition = factor(st_sub$condition),
    batch = factor(st_sub$batch)
  )

  dxd <- DEXSeqDataSetFromHTSeq(
    countfiles = countFiles_sub,
    sampleData = sampleTable,
    design = design_formula,
    flattenedfile = flat_gff
  )

  # light pre-filter to reduce memory before fitting
  cts_raw <- counts(dxd)
  keep <- rowSums(cts_raw) >= 10
  message("Exon bins before filter: ", nrow(dxd))
  message("Exon bins kept after filter: ", sum(keep))
  dxd <- dxd[keep, ]
  rm(cts_raw)
  gc()

  dxd <- estimateSizeFactors(dxd)
  gc()

  dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
  gc()

  dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
  gc()

  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = BPPARAM)
  gc()

  dxr <- DEXSeqResults(dxd)

  saveRDS(dxd, file.path(outdir, paste0("dxd_", tag, ".rds")))
  saveRDS(dxr, file.path(outdir, paste0("dxr_", tag, ".rds")))

  exons_df <- flatten_dxr(dxr)
  fwrite(
    as.data.table(exons_df),
    file.path(outdir, paste0("exons_all_", tag, ".tsv")),
    sep = "\t", quote = FALSE
  )

  exons_sig <- exons_df[!is.na(exons_df$padj) & exons_df$padj < padj_cut, , drop = FALSE]
  fwrite(
    as.data.table(exons_sig),
    file.path(outdir, paste0("exons_sig_padj", padj_cut, "_", tag, ".tsv")),
    sep = "\t", quote = FALSE
  )

  q <- perGeneQValue(dxr)
  genes_df <- data.frame(groupID = names(q), perGeneQValue = as.numeric(q), stringsAsFactors = FALSE)
  genes_df$pass_q <- !is.na(genes_df$perGeneQValue) & genes_df$perGeneQValue < gene_q_cut
  fwrite(
    as.data.table(genes_df),
    file.path(outdir, paste0("genes_qvalues_", tag, ".tsv")),
    sep = "\t", quote = FALSE
  )

  genes_sig <- genes_df[genes_df$pass_q, , drop = FALSE]
  fwrite(
    as.data.table(genes_sig),
    file.path(outdir, paste0("genes_sig_q", gene_q_cut, "_", tag, ".tsv")),
    sep = "\t", quote = FALSE
  )

  summary_df <- data.frame(
    tag = tag,
    design = paste(deparse(design_formula), collapse = ""),
    n_samples = nrow(st_sub),
    n_groups = length(unique(st_sub$condition)),
    n_exon_bins = nrow(exons_df),
    n_exons_sig_padj = nrow(exons_sig),
    n_genes = nrow(genes_df),
    n_genes_sig_q = sum(genes_df$pass_q),
    stringsAsFactors = FALSE
  )
  fwrite(
    as.data.table(summary_df),
    file.path(outdir, paste0("summary_", tag, ".tsv")),
    sep = "\t", quote = FALSE
  )

  pdf(file.path(outdir, paste0("dispersion_plot_", tag, ".pdf")))
  plotDispEsts(dxd)
  dev.off()

  # -------------------------------------------------------
  # HTML report using DEXSeqHTML default filtering
  # genes = NULL lets DEXSeqHTML choose significant hits at FDR
  # -------------------------------------------------------
  html_dir <- file.path(outdir, "DEXSeq_HTML_report")
  html_status <- file.path(outdir, "DEXSeq_HTML_status.txt")

  unlink(html_dir, recursive = TRUE, force = TRUE)
  dir.create(html_dir, showWarnings = FALSE, recursive = TRUE)

  tryCatch(
    {
      DEXSeqHTML(
        dxr,
        genes = NULL,
        FDR = padj_cut,
        color = c("#FF000080", "#0000FF80"),
        path = html_dir
      )
      writeLines(
        paste0(
          "DEXSeqHTML completed successfully using genes = NULL and FDR = ",
          padj_cut,
          "."
        ),
        html_status
      )
    },
    error = function(e) {
      msg <- paste0("DEXSeqHTML failed: ", conditionMessage(e))
      warning(msg)
      writeLines(msg, html_status)
    }
  )

  invisible(list(dxd = dxd, dxr = dxr))
}

# =========================
# RUN ANALYSIS
# =========================
if (force_results && dir.exists(outroot)) {
  message("FORCE_RESULTS=1: wiping outroot: ", outroot)
  unlink(outroot, recursive = TRUE, force = TRUE)
  dir.create(outroot, showWarnings = FALSE, recursive = TRUE)
}

countFiles_all <- file.path(counts_dir_clean, paste0(st$sample, ".dexseq.txt"))
stopifnot(all(file.exists(countFiles_all)))

BPP <- BiocParallel::SerialParam()

# 1) PCA (SELECT samples)
tag_pca <- "ALL_samples_select"
out_pca <- file.path(outroot, tag_pca, "PCA")
dir.create(out_pca, showWarnings = FALSE, recursive = TRUE)

pca_pdf <- file.path(out_pca, paste0("PCA_", tag_pca, "_PC1_PC2.pdf"))
if (!force_results && file.exists(pca_pdf)) {
  message("Skipping PCA (exists): ", pca_pdf)
} else {
  pca_res <- make_exon_usage_pca(st, countFiles_all, flat_gff, out_pca, tag_pca)
  message("PCA kept exon bins: ", pca_res$kept_bins)
}

# 2) Main DEXSeq: PPCD vs Control (SELECT; NO batch correction)
tag_main <- "PPCD_vs_Control_Select_noBatch"
out_main <- file.path(outroot, tag_main, "DEXSeq_noBatch")
dir.create(out_main, showWarnings = FALSE, recursive = TRUE)

sum_main <- file.path(out_main, paste0("summary_", tag_main, ".tsv"))
if (!force_results && file.exists(sum_main)) {
  message("Skipping main DEXSeq (exists): ", sum_main)
} else {
  message("\n=== Running DEXSeq: ", tag_main, " (NO batch; SELECT samples) ===")
  design_main <- ~ sample + exon + condition:exon
  run_dexseq(st, countFiles_all, flat_gff, out_main, tag_main, design_main, BPP, padj_cut, gene_q_cut)
  message("Done: ", tag_main)
}

writeLines(capture.output(sessionInfo()), file.path(outroot, "sessionInfo.txt"))

message("\nAll outputs in:\n", outroot)
RSCRIPT

chmod +x "$R_SCRIPT"

############################
# STEP 5: RUN R ANALYSIS
############################
echo "[4/5] R script written to: $R_SCRIPT"
echo "[5/5] Running DEXSeq R analysis ..."

export DEXROOT PADJ_CUT GENE_Q_CUT FORCE_CLEAN FORCE_RESULTS

Rscript "$R_SCRIPT"

echo
echo "Finished."
echo "Main output root:"
echo "  $RESULTS_ROOT"
echo
echo "Expected main DEXSeq result dir:"
echo "  $RESULTS_ROOT/PPCD_vs_Control_Select_noBatch/DEXSeq_noBatch"