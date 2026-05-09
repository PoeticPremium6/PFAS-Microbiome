#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$PHYLOSEQ_DIR"

OTU_TABLE="$PHYLOSEQ_DIR/otu_table.tsv"
TAX_TABLE="$PHYLOSEQ_DIR/taxonomy_table.tsv"
TMP_DIR="$PHYLOSEQ_DIR/tmp"
mkdir -p "$TMP_DIR"

mapfile -t BRACKEN_FILES < <(find "$TAXONOMY_DIR" -type f -name "*_bracken*.tsv" | sort)
if [ "${#BRACKEN_FILES[@]}" -eq 0 ]; then
  echo "No Bracken files found in $TAXONOMY_DIR" >&2
  exit 1
fi

# Convert each Bracken file to a simple table: taxon, count, taxonomy_id.
for f in "${BRACKEN_FILES[@]}"; do
  sample=$(basename "$f")
  sample=${sample%%_bracken*}
  awk -F '\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} NR>1{print $1"\t"$(h["new_est_reads"] ? h["new_est_reads"] : 6)"\t"$(h["taxonomy_id"] ? h["taxonomy_id"] : 0)}' "$f" \
    > "$TMP_DIR/${sample}.species.tsv"
done

# Build taxa list.
cut -f1 "$TMP_DIR"/*.species.tsv | sort -u > "$TMP_DIR/all_taxa.txt"

# Header.
printf "Taxa" > "$OTU_TABLE"
for f in "$TMP_DIR"/*.species.tsv; do
  sample=$(basename "$f" .species.tsv)
  printf "\t%s" "$sample" >> "$OTU_TABLE"
done
printf "\n" >> "$OTU_TABLE"

# Counts.
while IFS= read -r taxon; do
  printf "%s" "$taxon" >> "$OTU_TABLE"
  for f in "$TMP_DIR"/*.species.tsv; do
    count=$(awk -F '\t' -v t="$taxon" '$1==t{print $2; found=1; exit} END{if(!found) print 0}' "$f")
    printf "\t%s" "$count" >> "$OTU_TABLE"
  done
  printf "\n" >> "$OTU_TABLE"
done < "$TMP_DIR/all_taxa.txt"

# Minimal taxonomy table. This can be expanded with GTDB-Tk/NCBI taxonomy if needed.
printf "Taxa\tSpecies\n" > "$TAX_TABLE"
while IFS= read -r taxon; do
  printf "%s\t%s\n" "$taxon" "$taxon" >> "$TAX_TABLE"
done < "$TMP_DIR/all_taxa.txt"

echo "Phyloseq-ready tables written to: $PHYLOSEQ_DIR"
echo " - $OTU_TABLE"
echo " - $TAX_TABLE"
