#!/bin/bash
# Download TOGA orthology + loss tables for every selected assembly (6 parallel downloads)
# Usage: s3_fetch_toga.sh <toga_species.tsv> <dl_dir>
B=https://genome.senckenberg.de/download/TOGA/human_hg38_reference
DL=$2
tail -n +2 "$1" | cut -f5 | xargs -P 6 -I{} bash -c '
d="{}"; name=$(basename "$d"); mkdir -p "'$DL'/toga/$name"
for f in orthologsClassification.tsv.gz loss_summ_data.tsv.gz; do
  out="'$DL'/toga/$name/$f"
  [ -s "$out" ] && continue
  curl -s -f --retry 3 -m 900 "'$B'/$d/$f" -o "$out.part" && mv "$out.part" "$out" || echo "FAIL $d/$f"
done'
echo DONE
