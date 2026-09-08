#!/usr/bin/env bash

REV="/path/to/PFAS_mSystems_revision"
PH6="${REV}/05_humann_functional_profiles"
INROOT="${PH6}/01_sample_outputs"
MERGED="${PH6}/02_merged"
LINKS="${MERGED}/input_links"
QC="${PH6}/03_qc"

mkdir -p "$MERGED" "$LINKS" "$QC"

module purge || true
module load gcc12-env/12.1.0 || true
module load micromamba/1.4.2 || true
export MAMBA_ROOT_PREFIX=/path/to/hpc_work/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate PFAS_HUMANN312

rm -rf "$LINKS"
mkdir -p "$LINKS"

echo "=== counting sample-level HUMAnN files ==="
for TYPE in genefamilies pathabundance pathcoverage; do
  N=$(find "$INROOT" -mindepth 2 -maxdepth 2 -type f -name "*_${TYPE}.tsv" -size +0c | wc -l)
  echo "${TYPE}: ${N}"
done

GF=$(find "$INROOT" -mindepth 2 -maxdepth 2 -type f -name "*_genefamilies.tsv" -size +0c | wc -l)
PA=$(find "$INROOT" -mindepth 2 -maxdepth 2 -type f -name "*_pathabundance.tsv" -size +0c | wc -l)
PC=$(find "$INROOT" -mindepth 2 -maxdepth 2 -type f -name "*_pathcoverage.tsv" -size +0c | wc -l)

if [ "$GF" != "65" ] || [ "$PA" != "65" ] || [ "$PC" != "65" ]; then
  echo "ERROR: expected 65 files for each HUMAnN output type."
  exit 2
fi

find "$INROOT" -mindepth 2 -maxdepth 2 -type f \
  \( -name "*_genefamilies.tsv" -o -name "*_pathabundance.tsv" -o -name "*_pathcoverage.tsv" \) \
  -size +0c | while read F; do
    ln -sf "$F" "${LINKS}/$(basename "$F")"
  done

echo "=== linked files ==="
find "$LINKS" -type l | wc -l

rm -f "${MERGED}"/Ronneby_HUMAnN_*.tsv

humann_join_tables --input "$LINKS" --file_name genefamilies \
  --output "${MERGED}/Ronneby_HUMAnN_genefamilies.tsv"
[ -s "${MERGED}/Ronneby_HUMAnN_genefamilies.tsv" ] || { echo "ERROR: genefamilies merge failed"; exit 10; }

humann_join_tables --input "$LINKS" --file_name pathabundance \
  --output "${MERGED}/Ronneby_HUMAnN_pathabundance.tsv"
[ -s "${MERGED}/Ronneby_HUMAnN_pathabundance.tsv" ] || { echo "ERROR: pathabundance merge failed"; exit 11; }

humann_join_tables --input "$LINKS" --file_name pathcoverage \
  --output "${MERGED}/Ronneby_HUMAnN_pathcoverage.tsv"
[ -s "${MERGED}/Ronneby_HUMAnN_pathcoverage.tsv" ] || { echo "ERROR: pathcoverage merge failed"; exit 12; }

humann_renorm_table \
  --input "${MERGED}/Ronneby_HUMAnN_pathabundance.tsv" \
  --output "${MERGED}/Ronneby_HUMAnN_pathabundance_cpm.tsv" \
  --units cpm
[ -s "${MERGED}/Ronneby_HUMAnN_pathabundance_cpm.tsv" ] || { echo "ERROR: pathabundance CPM renorm failed"; exit 13; }

{
  echo -e "file\trows\tcolumns"
  for F in "${MERGED}"/Ronneby_HUMAnN_*.tsv; do
    [ -s "$F" ] || continue
    rows=$(awk 'END{print NR}' "$F")
    cols=$(awk -F'\t' 'NR==1{print NF; exit}' "$F")
    echo -e "$(basename "$F")\t${rows}\t${cols}"
  done
} > "${QC}/PHASE6_merged_HUMAnN_table_dimensions.tsv"

echo "MERGE_COMPLETE"
cat "${QC}/PHASE6_merged_HUMAnN_table_dimensions.tsv"
