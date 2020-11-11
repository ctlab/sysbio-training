#!/usr/bin/env bash
# Author roman.cherniatchik@jetbrains.com.com

set -euo pipefail

ALL_DMRS_BED=$1
BASENAME=$2
DMR_MIN_MC_CHANGE=$3
DMR_MIN_CCOUNT=$4

DMR_FILTERED="${BASENAME}.bed"
DMR_FILTERED_UP="${BASENAME}.up.bed"
DMR_FILTERED_DOWN="${BASENAME}.down.bed"

# E.g. for GREAT
DMR_FILTERED_BED3="${BASENAME}.bed3"
DMR_FILTERED_UP_BED3="${BASENAME}.up.bed3"
DMR_FILTERED_DOWN_BED3="${BASENAME}.down.bed3"

awk -v MCH="$DMR_MIN_MC_CHANGE" -v N_C="$DMR_MIN_CCOUNT" \
    'function abs(v) {{return v < 0 ? -v : v}} {{if (abs($6) >= MCH && $5 >= N_C) print}}' \
     "$ALL_DMRS_BED" | sort -k1,1 -k2,3n > "$DMR_FILTERED"

echo "Done: $DMR_FILTERED"
echo "  Filtered $(wc -l "$DMR_FILTERED" | awk '{ print $1 }')/$(wc -l "$ALL_DMRS_BED" | awk '{ print $1 }') regions" \
    "where each region has >= $DMR_MIN_CCOUNT and" \
    "abs(methylation change) >= $DMR_MIN_MC_CHANGE"
cut -f 1,2,3 "$DMR_FILTERED" > "$DMR_FILTERED_BED3"

awk '{{ if ($6 > 0) print }}' "$DMR_FILTERED" > "$DMR_FILTERED_UP"
echo "  UP:   $(wc -l "$DMR_FILTERED_UP" | awk '{ print $1 }')"
# shellcheck disable=SC2086
cut -f 1,2,3 $DMR_FILTERED_UP > "$DMR_FILTERED_UP_BED3"

awk '{{ if ($6 < 0) print }}' "$DMR_FILTERED" > "$DMR_FILTERED_DOWN"
echo "  DOWN: $(wc -l "$DMR_FILTERED_DOWN" | awk '{ print $1 }')"
cut -f 1,2,3 "$DMR_FILTERED_DOWN" > "$DMR_FILTERED_DOWN_BED3"