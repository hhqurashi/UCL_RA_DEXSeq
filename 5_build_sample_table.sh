#!/usr/bin/env bash
set -euo pipefail
: "${DEX:?DEX not set (source 2_set_paths_and_check_gtf.sh)}"
: "${PIPE:?PIPE not set (source 2_set_paths_and_check_gtf.sh)}"

mkdir -p "$DEX/merged_bams" "$DEX/logs"

out="$DEX/sample_table.tsv"
echo -e "sample\tcondition\tbatch\tbam" > "$out"

# Samples
samples=(
  C23 C38 PPCD1S1 PPCD1S2 PPCD1S3 PPCD1S4
  C18 C56 C59 PPCD1S5P0 PPCD3S1P0
  PPCD1SUK
)

# Batch map
declare -A batch=(
  [C23]=B1 [C38]=B1 [PPCD1S1]=B1 [PPCD1S2]=B1 [PPCD1S3]=B1 [PPCD1S4]=B1
  [C18]=B2 [C56]=B2 [C59]=B2 [PPCD1S5P0]=B2 [PPCD3S1P0]=B2
  [PPCD1SUK]=B3
)

# Condition map
declare -A cond=(
  [C23]=Control [C38]=Control [C18]=Control [C56]=Control [C59]=Control
  [PPCD1S1]=PPCD [PPCD1S2]=PPCD [PPCD1S3]=PPCD [PPCD1S4]=PPCD
  [PPCD1S5P0]=PPCD [PPCD3S1P0]=PPCD [PPCD1SUK]=PPCD
)

THREADS="${THREADS:-8}"

for s in "${samples[@]}"; do
  [[ -n "${batch[$s]:-}" ]] || { echo "Missing batch for $s"; exit 1; }
  [[ -n "${cond[$s]:-}"  ]] || { echo "Missing condition for $s"; exit 1; }

  bam="$PIPE/${s}.markdup.sorted.bam"
  if [[ ! -f "$bam" ]]; then
    # fallback search if naming differs slightly
    shopt -s nullglob
    matches=("$PIPE/${s}"*.markdup.sorted.bam)
    shopt -u nullglob
    if (( ${#matches[@]} == 1 )); then
      bam="${matches[0]}"
    else
      echo "Could not uniquely find BAM for sample $s in $PIPE"
      echo "Tried: $PIPE/${s}.markdup.sorted.bam and $PIPE/${s}*.markdup.sorted.bam"
      exit 1
    fi
  fi

  link="$DEX/merged_bams/${s}.markdup.sorted.bam"
  ln -sf "$bam" "$link"

  # index: prefer existing .bai if present, otherwise create
  if [[ -f "${bam}.bai" ]]; then
    ln -sf "${bam}.bai" "${link}.bai"
  elif [[ -f "${bam%.bam}.bai" ]]; then
    ln -sf "${bam%.bam}.bai" "${link}.bai"
  else
    samtools index -@ "$THREADS" "$link"
  fi

  echo -e "${s}\t${cond[$s]}\t${batch[$s]}\t${link}" >> "$out"
done

echo "Wrote: $out"
column -t "$out"
echo
echo "Condition counts:"
cut -f2 "$out" | tail -n +2 | sort | uniq -c
echo
echo "Batch counts:"
cut -f3 "$out" | tail -n +2 | sort | uniq -c