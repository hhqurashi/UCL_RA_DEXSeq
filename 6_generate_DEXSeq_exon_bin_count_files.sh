#!/usr/bin/env bash
set -euo pipefail
: "${DEX:?DEX not set (source 2_set_paths_and_check_gtf.sh)}"

ANN="${ANN:-$DEX/annotation/GRCh38.dexseq.gff}"
STRAND="${STRAND:-reverse}"   # "reverse", "yes", or "no"
NCORES="${NCORES:-6}"

cd "$DEX"

[[ -f "$DEX/sample_table.tsv" ]] || { echo "Missing: $DEX/sample_table.tsv"; exit 1; }
[[ -f "$ANN" ]] || { echo "Missing annotation GFF: $ANN"; exit 1; }

# --- NEW: locate dexseq_count.py shipped with the R DEXSeq package ---
DEXSEQ_COUNT="$(Rscript -e 'cat(system.file("python_scripts","dexseq_count.py",package="DEXSeq"))')"
[[ -n "$DEXSEQ_COUNT" && -f "$DEXSEQ_COUNT" ]] || {
  echo "Could not find dexseq_count.py via R DEXSeq. Got: '$DEXSEQ_COUNT'"
  exit 1
}
echo "Using dexseq_count.py: $DEXSEQ_COUNT"

mkdir -p "$DEX/counts" "$DEX/logs/count_logs"

# job list: sample<TAB>bam<TAB>out
awk -v DEXDIR="$DEX" 'BEGIN{OFS="\t"} NR>1 { print $1, $4, DEXDIR"/counts/"$1".dexseq.txt" }' \
  "$DEX/sample_table.tsv" > "$DEX/logs/count_jobs.tsv"

echo "Wrote job list: $DEX/logs/count_jobs.tsv"
head -n 5 "$DEX/logs/count_jobs.tsv" || true
echo

# Export vars so GNU parallel jobs can see them (subshells)
export DEXSEQ_COUNT ANN STRAND DEX

parallel -j "$NCORES" --colsep '\t' '
  set -euo pipefail
  sample={1}; bam={2}; out={3};

  # --- NEW: skip samples with an existing non-empty count file ---
  if [[ -s "$out" ]]; then
    echo "[skip] ${sample} (counts already present: $out)"
    exit 0
  fi

  echo "[count] ${sample}"
  python "$DEXSEQ_COUNT" -p yes -r pos -s "$STRAND" "$ANN" "$bam" "$out" \
    > "$DEX/logs/count_logs/${sample}.log" 2>&1
' :::: "$DEX/logs/count_jobs.tsv"

echo
echo "Counts written to: $DEX/counts"
ls -lh "$DEX/counts" | head