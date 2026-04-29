#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# PPCD vs Control SELECT DEXSeq workflow (lean + portable version)
#
# Lean main analysis:
#   - no PCA
#   - no batch correction
#   - no full exon export / flattening
#   - no RDS saves
#   - no dispersion plot
#   - no low-count filter
#
# Keeps only:
#   - main DEXSeq fit
#   - genes_qvalues TSV
#   - genes_sig_q TSV
#   - DEXSeq HTML report
#
# Portable behavior:
#   - if counts_clean exists, use it
#   - if counts_clean is missing but raw counts exist, clean them
#   - if both are missing, create raw counts from BAMs and then clean them
#
# Assumes:
#   - conda env "dexseq_grch38" is activated
#   - dexseq_prepare_annotation.py should use -r no
###############################################################################

############################
# USER CONFIG
############################
STAR_SALMON_DIR="/media/pontikos_nas2/HarisQurashi/projects/12_Nihar_DEXSeq/PPCD_vs_Control/outputs/star_salmon"
GTF_GZ="/media/pontikos_nas2/HarisQurashi/refs/Homo_sapiens.GRCh38.115.gtf.gz"
OUT_BASE="/mnt/scratch/hqurashi/12_Nihar_DEXSeq_PPCD_vs_Control/dexseq_NEW"

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
R_SCRIPT="$SCRIPT_DIR/run_DEXSeq_PPCD_vs_Control_Select_noBatch_lean.R"

cat > "$R_SCRIPT" << 'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DEXSeq)
  library(BiocParallel)
  library(data.table)
})

# =========================
# CONFIG
# =========================
DEXROOT <- Sys.getenv("DEXROOT")
if (nchar(DEXROOT) == 0) stop("DEXROOT env var not set")

sample_tsv <- file.path(DEXROOT, "sample_table.tsv")
counts_dir_clean <- file.path(DEXROOT, "counts_clean")
flat_gff <- file.path(DEXROOT, "annotation", "GRCh38.115.dexseq.r_no.gff")
outroot <- file.path(DEXROOT, "results", "DEXSeq_PPCD_vs_Control_Select")

padj_cut   <- as.numeric(Sys.getenv("PADJ_CUT", "0.05"))
gene_q_cut <- as.numeric(Sys.getenv("GENE_Q_CUT", "0.05"))
force_results <- as.integer(Sys.getenv("FORCE_RESULTS", "0")) == 1

message("DEXROOT          = ", DEXROOT)
message("sample_tsv       = ", sample_tsv)
message("counts_dir_clean = ", counts_dir_clean)
message("flat_gff         = ", flat_gff)
message("outroot          = ", outroot)
message("PADJ_CUT         = ", padj_cut)
message("GENE_Q_CUT       = ", gene_q_cut)
message("FORCE_RESULTS    = ", force_results)

stopifnot(file.exists(sample_tsv))
stopifnot(dir.exists(counts_dir_clean))
stopifnot(file.exists(flat_gff))

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

  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
  dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = BPPARAM)

  dxr <- DEXSeqResults(dxd)

  q <- perGeneQValue(dxr)
  genes_df <- data.frame(
    groupID = names(q),
    perGeneQValue = as.numeric(q),
    stringsAsFactors = FALSE
  )
  genes_df$pass_q <- !is.na(genes_df$perGeneQValue) & genes_df$perGeneQValue < gene_q_cut

  fwrite(
    as.data.table(genes_df),
    file.path(outdir, paste0("genes_qvalues_", tag, ".tsv")),
    sep = "\t",
    quote = FALSE
  )

  genes_sig <- genes_df[genes_df$pass_q, , drop = FALSE]
  fwrite(
    as.data.table(genes_sig),
    file.path(outdir, paste0("genes_sig_q", gene_q_cut, "_", tag, ".tsv")),
    sep = "\t",
    quote = FALSE
  )

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

  invisible(dxr)
}

tag_main <- "PPCD_vs_Control_Select_noBatch"
out_main <- file.path(outroot, tag_main, "DEXSeq_noBatch")
dir.create(out_main, showWarnings = FALSE, recursive = TRUE)

html_status <- file.path(out_main, "DEXSeq_HTML_status.txt")
if (!force_results && file.exists(html_status)) {
  message("Skipping main DEXSeq (HTML status exists): ", html_status)
} else {
  message("\n=== Running DEXSeq: ", tag_main, " (NO batch; SELECT samples) ===")
  design_main <- ~ sample + exon + condition:exon
  run_dexseq(st, countFiles_all, flat_gff, out_main, tag_main, design_main, BPP, padj_cut, gene_q_cut)
  message("Done: ", tag_main)
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
echo "Expected main DEXSeq result dir:"
echo "  $RESULTS_ROOT/PPCD_vs_Control_Select_noBatch/DEXSeq_noBatch"