cd $PIPE
ls *.markdup.sorted.bam > $DEX/logs/lane_bams.list

# map lane -> core sample (remove trailing _L<digit>)
> $DEX/logs/core_to_lane_map.tsv
while read -r bam; do
  s=$(basename "$bam" .markdup.sorted.bam)
  core=$(echo "$s" | sed -E 's/_L[0-9]+$//')
  echo -e "$core\t$s" >> $DEX/logs/core_to_lane_map.tsv
done < $DEX/logs/lane_bams.list

cut -f1 $DEX/logs/core_to_lane_map.tsv | sort -u > $DEX/logs/core_samples.list

THREADS=8

while read -r core; do
  lanes=$(awk -v c="$core" '$1==c{print $2}' $DEX/logs/core_to_lane_map.tsv)

  out=$DEX/merged_bams/${core}.markdup.sorted.bam

  n=$(echo "$lanes" | wc -w)
  if [ "$n" -gt 1 ]; then
    echo "[merge] $core ($n lanes)"
    samtools merge -@ $THREADS -o "$out" $(for l in $lanes; do echo $PIPE/${l}.markdup.sorted.bam; done)
  else
    echo "[link]  $core (1 lane)"
    ln -sf "$PIPE/${lanes}.markdup.sorted.bam" "$out"
  fi

  samtools index -@ $THREADS "$out"
done < $DEX/logs/core_samples.list