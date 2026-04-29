#!/usr/bin/env bash
set -euo pipefail
: "${DEX:?DEX not set (source 2_set_paths_and_check_gtf.sh)}"
: "${GTF:?GTF not set (source 2_set_paths_and_check_gtf.sh)}"

OUT="$DEX/annotation/GRCh38.dexseq.gff"

# If you already have a known-good one elsewhere, just symlink it:
KNOWN="/media/pontikos_nas2/HarisQurashi/projects/12_Nihar_DEXSeq/RNAseqOct2025/outputs/DEXSeq/annotation/GRCh38.dexseq.gff"

if [[ -f "$OUT" ]]; then
  echo "Annotation already exists: $OUT"
  exit 0
fi

if [[ -f "$KNOWN" ]]; then
  echo "Linking existing flattened GFF:"
  echo "  $KNOWN"
  ln -sf "$KNOWN" "$OUT"
  exit 0
fi

echo "Generating flattened GFF (no existing one found)..."
dexseq_prepare_annotation.py "$GTF" "$OUT"
echo "Wrote: $OUT"