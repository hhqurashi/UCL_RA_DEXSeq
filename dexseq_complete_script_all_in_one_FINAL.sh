#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# PPCD vs Control SELECT DEXSeq workflow (portable + PCA + batch model only)
#
# Runs:
#   1) Exon-usage PCA
#   2) DEXSeq batch-corrected model only
#
# Portable behavior:
#   - if counts_clean exists, use it
#   - if counts_clean is missing but raw counts exist, clean them
#   - if both are missing, create raw counts from BAMs and then clean them
#
# Outputs:
#   - PCA plots / coordinates
#   - saved dxr RDS for later exon-level export
#   - genes_qvalues TSV (annotated with gene names)
#   - genes_sig_q TSV (annotated with gene names)
#   - DEXSeq HTML report with extra gene_symbol column
#   - gene_id_to_name lookup table from the same GTF used for DEXSeq
#
# Notes:
#   - HTML keeps Ensembl IDs internally, but adds gene_symbol as an extra column
#   - TSVs are annotated with gene_name/gene_label from the same Ensembl GTF
#   - exon_sig export is done separately from the saved dxr object
###############################################################################

############################
# USER CONFIG
############################
STAR_SALMON_DIR="/media/pontikos_nas2/HarisQurashi/projects/12_Nihar_DEXSeq/PPCD_vs_Control/outputs/star_salmon"
GTF_GZ="/media/pontikos_nas2/HarisQurashi/refs/Homo_sapiens.GRCh38.115.gtf.gz"
OUT_BASE="/mnt/scratch/hqurashi/12_Nihar_DEXSeq_PPCD_vs_Control/dexseq_FINAL"

STRAND="reverse"         # reverse / yes / no
COUNT_THREADS=6          # used only if raw count generation is needed

FORCE_ANNOTATION=0       # 1 = remake local GTF + flattened GFF
FORCE_COUNTS=0           # 1 = rerun raw dexseq_count.py output
FORCE_CLEAN=0            # 1 = recreate counts_clean from raw counts
FORCE_RESULTS=0          # 1 = wipe DEXSeq results dir before R analysis

PADJ_CUT="0.05"
GENE_Q_CUT="0.05"

############################
# CHECKS
############################
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found. Activate dexseq_grch38."; exit 1; }
command -v python  >/dev/null 2>&1 || { echo "python not found. Activate dexseq_grch38."; exit 1; }
command -v awk     >/dev/null 2>&1 || { echo "awk not found."; exit 1; }

DEXSEQ_PREPARE="$(Rscript -e 'cat(system.file("python_scripts","dexseq_prepare_annotation.py", package="DEXSeq"))')"
[[ -n "$DEXSEQ_PREPARE" && -f "$DEXSEQ_PREPARE" ]] || {
  echo "Could not find dexseq_prepare_annotation.py via R DEXSeq. Got: '$DEXSEQ_PREPARE'"
  exit 1
}

DEXSEQ_COUNT="$(Rscript -e 'cat(system.file("python_scripts","dexseq_count.py", package="DEXSeq"))')"
[[ -n "$DEXSEQ_COUNT" && -f "$DEXSEQ_COUNT" ]] || {
  echo "Could not find dexseq_count.py via R DEXSeq. Got: '$DEXSEQ_COUNT'"
  exit 1
}

[[ -d "$STAR_SALMON_DIR" ]] || { echo "Missing STAR_SALMON_DIR: $STAR_SALMON_DIR"; exit 1; }
[[ -f "$GTF_GZ" ]] || { echo "Missing GTF_GZ: $GTF_GZ"; exit 1; }

DEXROOT="$OUT_BASE"
ANN_DIR="$DEXROOT/annotation"
COUNTS_DIR="$DEXROOT/counts"
COUNTS_CLEAN_DIR="$DEXROOT/counts_clean"
LOG_DIR="$DEXROOT/logs"
COUNT_LOG_DIR="$LOG_DIR/count_logs"
RESULTS_ROOT="$DEXROOT/results/DEXSeq_PPCD_vs_Control_Select"
SCRIPT_DIR="$DEXROOT/scripts_generated"

mkdir -p "$ANN_DIR" "$COUNTS_DIR" "$COUNTS_CLEAN_DIR" \
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
  echo "[1/4] Decompressing GTF to local annotation directory..."
  gzip -dc "$GTF_GZ" > "$GTF_LOCAL"
else
  echo "[1/4] Reusing local GTF: $GTF_LOCAL"
fi

echo -n "Local GTF exon count: "
awk '$3=="exon"{c++} END{print c+0}' "$GTF_LOCAL"

if [[ ! -f "$FLAT_GFF" ]]; then
  echo "[1/4] Creating flattened DEXSeq GFF with -r no ..."
  python "$DEXSEQ_PREPARE" -r no "$GTF_LOCAL" "$FLAT_GFF"
else
  echo "[1/4] Reusing flattened GFF: $FLAT_GFF"
fi
echo

############################
# STEP 2: BUILD SAMPLE TABLE + RESOLVE BAMs
############################
SAMPLE_TSV="$DEXROOT/sample_table.tsv"
COUNT_JOBS_TSV="$LOG_DIR/count_jobs.tsv"

echo "[2/4] Building sample table for SELECT samples ..."

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

    printf "%s\t%s\t%s\t%s\n" "$s" "${cond[$s]}" "${batch[$s]}" "$bam"
  done
} > "$SAMPLE_TSV"

echo "Wrote: $SAMPLE_TSV"
if command -v column >/dev/null 2>&1; then
  column -t "$SAMPLE_TSV"
else
  cat "$SAMPLE_TSV"
fi
echo

awk -v DEXDIR="$DEXROOT" 'BEGIN{OFS="\t"} NR>1 { print $1, $4, DEXDIR"/counts/"$1".dexseq.txt", DEXDIR"/counts_clean/"$1".dexseq.txt" }' \
  "$SAMPLE_TSV" > "$COUNT_JOBS_TSV"

echo "Wrote count job list: $COUNT_JOBS_TSV"
head -n 5 "$COUNT_JOBS_TSV" || true
echo

############################
# STEP 3: ENSURE CLEAN COUNT FILES EXIST
############################
clean_one_countfile() {
  local infile="$1"
  local outfile="$2"

  awk '
    BEGIN{OFS="\t"}
    {
      if (NF < 2) next

      feature=$1
      count=$2

      gsub(/"/, "", feature)
      gsub(/"/, "", count)

      if (feature == "" || count == "") next
      if (feature ~ /^_/) next
      if (index(feature, ":") == 0) next
      if (count !~ /^-?[0-9]+$/) next

      print feature, count
    }
  ' "$infile" > "$outfile"

  if [[ ! -s "$outfile" ]]; then
    echo "Cleaning produced empty output: $outfile" >&2
    return 1
  fi
}

validate_clean_countfile() {
  local f="$1"
  awk '
    BEGIN{ok=1}
    NF != 2 { ok=0; exit }
    $1 == "" { ok=0; exit }
    $2 == "" { ok=0; exit }
    $1 ~ /^_/ { ok=0; exit }
    index($1, ":") == 0 { ok=0; exit }
    $2 !~ /^-?[0-9]+$/ { ok=0; exit }
    END{
      if (NR == 0) ok=0
      exit(ok ? 0 : 1)
    }
  ' "$f"
}

echo "[3/4] Ensuring clean count files exist ..."

need_raw_counts=0
need_cleaning=0

if [[ "$FORCE_COUNTS" == "1" ]]; then
  need_raw_counts=1
fi

if [[ "$FORCE_CLEAN" == "1" ]]; then
  need_cleaning=1
fi

while IFS=$'\t' read -r sample bam raw clean; do
  if [[ ! -s "$clean" ]]; then
    need_cleaning=1
  fi
  if [[ ! -s "$raw" ]]; then
    need_raw_counts=1
  fi
done < "$COUNT_JOBS_TSV"

if [[ "$need_raw_counts" == "1" ]]; then
  command -v parallel >/dev/null 2>&1 || {
    echo "GNU parallel not found, but raw count generation is needed. Activate dexseq_grch38." >&2
    exit 1
  }

  echo "Raw count generation required."
  if [[ "$FORCE_COUNTS" == "1" ]]; then
    echo "FORCE_COUNTS=1, removing existing raw count files ..."
    while IFS=$'\t' read -r sample bam raw clean; do
      rm -f "$raw"
    done < "$COUNT_JOBS_TSV"
  fi

  export DEXSEQ_COUNT FLAT_GFF STRAND COUNT_LOG_DIR

  parallel -j "$COUNT_THREADS" --colsep '\t' '
    set -euo pipefail
    sample={1}
    bam={2}
    raw={3}

    if [[ -s "$raw" ]]; then
      echo "[skip raw] ${sample}"
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
      "$raw" \
      > "$COUNT_LOG_DIR/${sample}.log" 2>&1
  ' :::: "$COUNT_JOBS_TSV"

  echo "Raw counts ready in: $COUNTS_DIR"
  echo
else
  echo "Raw count generation not needed."
fi

if [[ "$need_cleaning" == "1" ]]; then
  echo "Cleaning raw counts -> counts_clean ..."

  if [[ "$FORCE_CLEAN" == "1" ]]; then
    echo "FORCE_CLEAN=1, removing existing clean count files ..."
    while IFS=$'\t' read -r sample bam raw clean; do
      rm -f "$clean"
    done < "$COUNT_JOBS_TSV"
  fi

  while IFS=$'\t' read -r sample bam raw clean; do
    if [[ -s "$clean" && "$FORCE_CLEAN" != "1" ]]; then
      echo "[skip clean] $sample"
      continue
    fi

    [[ -s "$raw" ]] || { echo "Missing raw count file for cleaning: $raw" >&2; exit 1; }

    echo "[clean] $sample"
    clean_one_countfile "$raw" "$clean"
    validate_clean_countfile "$clean" || {
      echo "Validation failed for cleaned count file: $clean" >&2
      exit 1
    }
  done < "$COUNT_JOBS_TSV"

  echo "Clean counts ready in: $COUNTS_CLEAN_DIR"
  echo
else
  echo "Cleaning not needed."
  echo
fi

echo "Final clean count file check:"
while IFS=$'\t' read -r sample bam raw clean; do
  [[ -s "$clean" ]] || { echo "Missing final clean count file: $clean" >&2; exit 1; }
done < "$COUNT_JOBS_TSV"
echo "All required clean count files found."
echo

############################
# STEP 4: WRITE + RUN R SCRIPT
############################
R_SCRIPT="$SCRIPT_DIR/run_DEXSeq_PPCD_vs_Control_Select_batchOnly.R"

cat > "$R_SCRIPT" << 'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DEXSeq)
  library(BiocParallel)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(rlang)
})

# =========================
# CONFIG
# =========================
DEXROOT <- Sys.getenv("DEXROOT")
if (nchar(DEXROOT) == 0) stop("DEXROOT env var not set")

sample_tsv <- file.path(DEXROOT, "sample_table.tsv")
counts_dir_clean <- file.path(DEXROOT, "counts_clean")
flat_gff <- file.path(DEXROOT, "annotation", "GRCh38.115.dexseq.r_no.gff")
gtf_local <- file.path(DEXROOT, "annotation", "Homo_sapiens.GRCh38.115.gtf")
outroot <- file.path(DEXROOT, "results", "DEXSeq_PPCD_vs_Control_Select")

padj_cut   <- as.numeric(Sys.getenv("PADJ_CUT", "0.05"))
gene_q_cut <- as.numeric(Sys.getenv("GENE_Q_CUT", "0.05"))
force_results <- as.integer(Sys.getenv("FORCE_RESULTS", "0")) == 1

message("DEXROOT          = ", DEXROOT)
message("sample_tsv       = ", sample_tsv)
message("counts_dir_clean = ", counts_dir_clean)
message("flat_gff         = ", flat_gff)
message("gtf_local        = ", gtf_local)
message("outroot          = ", outroot)
message("PADJ_CUT         = ", padj_cut)
message("GENE_Q_CUT       = ", gene_q_cut)
message("FORCE_RESULTS    = ", force_results)

stopifnot(file.exists(sample_tsv))
stopifnot(dir.exists(counts_dir_clean))
stopifnot(file.exists(flat_gff))
stopifnot(file.exists(gtf_local))

if (force_results && dir.exists(outroot)) {
  message("FORCE_RESULTS=1: wiping outroot: ", outroot)
  unlink(outroot, recursive = TRUE, force = TRUE)
}
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)

st <- read.delim(sample_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("sample", "condition", "batch", "bam") %in% colnames(st)))

st$condition <- factor(st$condition, levels = c("Control", "PPCD"))
st$batch <- factor(st$batch, levels = unique(st$batch))

message("Samples: ", paste(st$sample, collapse = ", "))
message("Condition x batch:\n", paste(capture.output(print(table(st$condition, st$batch))), collapse = "\n"))

countFiles_all <- file.path(counts_dir_clean, paste0(st$sample, ".dexseq.txt"))
stopifnot(all(file.exists(countFiles_all)))

BPP <- BiocParallel::SerialParam()

# =========================
# Gene ID -> gene name mapping from SAME GTF used for DEXSeq
# =========================
extract_gtf_attr <- function(x, key) {
  m <- regexec(paste0(key, ' "([^"]+)"'), x)
  reg <- regmatches(x, m)
  vapply(reg, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
}

make_gene_map <- function(gtf_file) {
  gtf <- data.table::fread(
    cmd = paste("grep -v '^#' ", shQuote(gtf_file)),
    sep = "\t",
    header = FALSE,
    quote = "",
    data.table = FALSE,
    showProgress = FALSE
  )

  if (ncol(gtf) < 9) stop("GTF parse failed: expected >= 9 columns in ", gtf_file)
  colnames(gtf)[1:9] <- c("seqname","source","feature","start","end","score","strand","frame","attribute")

  if (any(gtf$feature == "gene")) {
    gtf <- gtf[gtf$feature == "gene", , drop = FALSE]
  }

  gene_id <- extract_gtf_attr(gtf$attribute, "gene_id")
  gene_name <- extract_gtf_attr(gtf$attribute, "gene_name")

  gm <- data.frame(
    groupID_raw = gene_id,
    groupID_base = sub("\\.[0-9]+$", "", gene_id),
    gene_name = gene_name,
    stringsAsFactors = FALSE
  )

  gm <- gm[!is.na(gm$groupID_raw) & nzchar(gm$groupID_raw), , drop = FALSE]
  gm <- gm[!duplicated(gm$groupID_raw), , drop = FALSE]
  gm$gene_label <- ifelse(!is.na(gm$gene_name) & nzchar(gm$gene_name), gm$gene_name, gm$groupID_raw)

  gm
}

annotate_gene_ids <- function(df, gene_map, id_col = "groupID") {
  if (!id_col %in% colnames(df)) return(df)

  idx_raw <- match(df[[id_col]], gene_map$groupID_raw)
  idx_base <- match(sub("\\.[0-9]+$", "", df[[id_col]]), gene_map$groupID_base)
  idx <- ifelse(!is.na(idx_raw), idx_raw, idx_base)

  df$gene_name  <- gene_map$gene_name[idx]
  df$gene_label <- ifelse(!is.na(df$gene_name) & nzchar(df$gene_name), df$gene_name, df[[id_col]])

  front <- c("gene_label", "gene_name", id_col)
  rest <- setdiff(colnames(df), front)
  df[, c(front, rest), drop = FALSE]
}

make_html_extra_cols <- function(dxr, gene_map) {
  gene_ids <- unique(as.character(dxr$groupID))

  idx_raw <- match(gene_ids, gene_map$groupID_raw)
  idx_base <- match(sub("\\.[0-9]+$", "", gene_ids), gene_map$groupID_base)
  idx <- ifelse(!is.na(idx_raw), idx_raw, idx_base)

  gene_symbol <- gene_map$gene_name[idx]
  gene_symbol <- ifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, gene_ids)

  extra_cols <- data.frame(
    gene_symbol = gene_symbol,
    stringsAsFactors = FALSE,
    row.names = gene_ids
  )

  extra_cols
}

gene_map <- make_gene_map(gtf_local)
data.table::fwrite(
  data.table::as.data.table(gene_map[, c("gene_label","gene_name","groupID_raw","groupID_base")]),
  file.path(outroot, "gene_id_to_name_from_gtf.tsv"),
  sep = "\t",
  quote = FALSE
)

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
    return(invisible(nrow(logit)))
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

  kept_bins <- nrow(logit)

  rm(dxd, cn, nc, gene, gene_sum, idx, usage, usage_f, pp, logit, pca, pcs, df)
  gc()

  invisible(kept_bins)
}

# =========================
# DEXSeq runner (batch-corrected)
# =========================
run_dexseq_batch <- function(
  st_sub,
  countFiles_sub,
  flat_gff,
  outdir,
  tag,
  full_model,
  reduced_model,
  fit_var,
  BPPARAM,
  padj_cut,
  gene_q_cut,
  gene_map
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  sampleTable <- data.frame(
    row.names = st_sub$sample,
    condition = factor(st_sub$condition),
    batch = factor(st_sub$batch)
  )

  dxd <- DEXSeqDataSetFromHTSeq(
    countfiles = countFiles_sub,
    sampleData = sampleTable,
    design = full_model,
    flattenedfile = flat_gff
  )

  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, formula = full_model, BPPARAM = BPPARAM)
  dxd <- testForDEU(
    dxd,
    reducedModel = reduced_model,
    fullModel = full_model,
    BPPARAM = BPPARAM
  )
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = fit_var, BPPARAM = BPPARAM)

  dxr <- DEXSeqResults(dxd)
  saveRDS(dxr, file.path(outdir, paste0("dxr_", tag, ".rds")))
  rm(dxd)
  gc()

  q <- perGeneQValue(dxr)
  genes_df <- data.frame(
    groupID = names(q),
    perGeneQValue = as.numeric(q),
    stringsAsFactors = FALSE
  )
  genes_df$pass_q <- !is.na(genes_df$perGeneQValue) & genes_df$perGeneQValue < gene_q_cut
  genes_df <- annotate_gene_ids(genes_df, gene_map, id_col = "groupID")

  data.table::fwrite(
    data.table::as.data.table(genes_df),
    file.path(outdir, paste0("genes_qvalues_", tag, ".tsv")),
    sep = "\t",
    quote = FALSE
  )

  genes_sig <- genes_df[genes_df$pass_q, , drop = FALSE]
  data.table::fwrite(
    data.table::as.data.table(genes_sig),
    file.path(outdir, paste0("genes_sig_q", gene_q_cut, "_", tag, ".tsv")),
    sep = "\t",
    quote = FALSE
  )

  extra_cols <- make_html_extra_cols(dxr, gene_map)

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
        path = html_dir,
        extraCols = extra_cols
      )
      writeLines(
        paste0(
          "DEXSeqHTML completed successfully using genes = NULL, FDR = ",
          padj_cut,
          ", and extraCols gene_symbol."
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

  rm(dxr, q, genes_df, genes_sig, extra_cols)
  gc()

  invisible(NULL)
}

# =========================
# RUN PCA
# =========================
tag_pca <- "ALL_samples_select"
out_pca <- file.path(outroot, tag_pca, "PCA")
dir.create(out_pca, showWarnings = FALSE, recursive = TRUE)

pca_pdf <- file.path(out_pca, paste0("PCA_", tag_pca, "_PC1_PC2.pdf"))
if (!force_results && file.exists(pca_pdf)) {
  message("Skipping PCA (exists): ", pca_pdf)
} else {
  message("\n=== Running exon-usage PCA: ", tag_pca, " ===")
  kept_bins <- make_exon_usage_pca(st, countFiles_all, flat_gff, out_pca, tag_pca)
  message("PCA kept exon bins: ", kept_bins)
}

# =========================
# RUN BATCH-CORRECTED MODEL ONLY
# =========================
tag_batch <- "PPCD_vs_Control_Select_withBatch"
out_batch <- file.path(outroot, tag_batch, "DEXSeq_withBatch")
dir.create(out_batch, showWarnings = FALSE, recursive = TRUE)

html_status_batch <- file.path(out_batch, "DEXSeq_HTML_status.txt")
if (!force_results && file.exists(html_status_batch)) {
  message("Skipping batch-corrected DEXSeq (HTML status exists): ", html_status_batch)
} else {
  message("\n=== Running DEXSeq: ", tag_batch, " ===")
  full_batch <- ~ sample + exon + batch:exon + condition:exon
  reduced_batch <- ~ sample + exon + batch:exon

  run_dexseq_batch(
    st_sub = st,
    countFiles_sub = countFiles_all,
    flat_gff = flat_gff,
    outdir = out_batch,
    tag = tag_batch,
    full_model = full_batch,
    reduced_model = reduced_batch,
    fit_var = "condition",
    BPPARAM = BPP,
    padj_cut = padj_cut,
    gene_q_cut = gene_q_cut,
    gene_map = gene_map
  )
  message("Done: ", tag_batch)
}

message("\nAll outputs in:\n", outroot)
RSCRIPT

chmod +x "$R_SCRIPT"

echo "[4/4] R script written to: $R_SCRIPT"
echo "Running DEXSeq R analysis ..."

export DEXROOT PADJ_CUT GENE_Q_CUT FORCE_RESULTS

Rscript "$R_SCRIPT"

echo
echo "Finished."
echo "Main output root:"
echo "  $RESULTS_ROOT"
echo
echo "Expected result dirs:"
echo "  $RESULTS_ROOT/ALL_samples_select/PCA"
echo "  $RESULTS_ROOT/PPCD_vs_Control_Select_withBatch/DEXSeq_withBatch"