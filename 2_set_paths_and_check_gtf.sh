#!/usr/bin/env bash

# MASTER output dir for this analysis (everything goes here)
export BASE="/mnt/scratch/hqurashi/12_Nihar_DEXSeq_PPCD_vs_Control"

# Where the BAMs live for THIS analysis
export PIPE="/media/pontikos_nas2/HarisQurashi/projects/12_Nihar_DEXSeq/PPCD_vs_Control/outputs/star_salmon"

# Reuse the same GTF you used previously (GRCh38); change only if you truly have a different one
export GTF="/mnt/scratch/hqurashi/12_Nihar_DEXSeq/outputs/genome/genes.filtered.gtf"

# DEXSeq working directory (inside MASTER)
export DEX="$BASE/dexseq"

mkdir -p "$DEX"/{annotation,merged_bams,counts,counts_fmt3,results,logs}
mkdir -p "$DEX/logs"/count_logs

echo "BASE=$BASE"
echo "PIPE=$PIPE"
echo "GTF=$GTF"
echo "DEX=$DEX"
echo

[[ -d "$PIPE" ]] || { echo "Missing PIPE dir: $PIPE"; return 1 2>/dev/null || exit 1; }
[[ -f "$GTF"  ]]  || { echo "Missing GTF file: $GTF"; return 1 2>/dev/null || exit 1; }

echo -n "exon_features= "
awk '$3=="exon"{c++} END{print c+0}' "$GTF"

BAM_EXAMPLE=$(ls -1 "$PIPE"/*.markdup.sorted.bam 2>/dev/null | head -n 1 || true)
echo "BAM example: ${BAM_EXAMPLE:-NONE FOUND}"
if [[ -n "${BAM_EXAMPLE:-}" ]]; then
  samtools view -H "$BAM_EXAMPLE" | grep '^@SQ' | head -n 3
fi