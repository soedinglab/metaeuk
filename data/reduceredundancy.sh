#!/bin/bash -e

# reduce redundancy workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [[ "$1" == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
            echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <dp_protein_contig_strand_map> <dp_optimal_exon_sets> <metaeukBaseName> <tmpDir>" && exit 1;
# check if file exists
INPUT_MAP="$(abspath "$1")"
INPUT_OPTIMAL_EXON_SETS="$(abspath "$2")"

[ ! -f "${INPUT_MAP}" ] &&  echo "${INPUT_MAP} not found!" && exit 1;
[ ! -f "${INPUT_OPTIMAL_EXON_SETS}" ] &&  echo "${INPUT_OPTIMAL_EXON_SETS} not found!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $2 not found!" && mkdir -p "$4";


TMP_PATH="$(abspath "$4")"

# aggregate all TCS (Target + Contig + Strand) predictions by their CS
if notExists "${TMP_PATH}/dp_contig_strand_map"; then
    "$MMSEQS" swapdb "${INPUT_MAP}" "${TMP_PATH}/dp_contig_strand_map" \
        || fail "swapdb step died"
fi

# the CS map record contains 10 columns: proteinContigStrandId, proteinMMSeqs2Key, contigMMSeqs2Key, strand, combinedBitScore, combinedEvalue, numExons, lowContigCoord, highContigCoord, exonIDsStr.

## the following three steps create a primary, secondary and tertiary order:

# within each CS, sort by bitscore (decreasing order), tertiary order
if notExists "${TMP_PATH}/dp_contig_strand_map_sorted_by_bitscore"; then
    "$MMSEQS" filterdb "${TMP_PATH}/dp_contig_strand_map" "${TMP_PATH}/dp_contig_strand_map_sorted_by_bitscore" --sort-entries 2 --filter-column 5 \
        || fail "filterdb (to sort by bitscore in decreasing order) step died"
fi

# within each CS, sort by number of exons (decreasing order), secondary order
if notExists "${TMP_PATH}/dp_contig_strand_map_sorted_by_num_exons_subsorted_by_bit"; then
    "$MMSEQS" filterdb "${TMP_PATH}/dp_contig_strand_map_sorted_by_bitscore" "${TMP_PATH}/dp_contig_strand_map_sorted_by_num_exons_subsorted_by_bit" --sort-entries 2 --filter-column 7 \
        || fail "filterdb (to sort by number of exons in decreasing order) step died"
fi

# within each CS, sort by start position on the contig (increasing order), primary order
if notExists "${TMP_PATH}/dp_contig_strand_map_sorted_by_start_subsorted_by_num_and_bit"; then
    "$MMSEQS" filterdb "${TMP_PATH}/dp_contig_strand_map_sorted_by_num_exons_subsorted_by_bit" "${TMP_PATH}/dp_contig_strand_map_sorted_by_start_subsorted_by_num_and_bit" --sort-entries 1 --filter-column 8 \
        || fail "filterdb (to sort by start position on the contig in increasing order) step died"
fi

# greedy grouping of predictions based on the sorted map
if notExists "${TMP_PATH}/grouped_predictions"; then
    "$MMSEQS" grouppredictions "${TMP_PATH}/dp_contig_strand_map_sorted_by_start_subsorted_by_num_and_bit" "${TMP_PATH}/grouped_predictions" "${TMP_PATH}/grouped_predictions_no_overlap" \
        || fail "grouppredictions step died"
fi

mv -f "${TMP_PATH}/grouped_predictions" "$3_grouped_predictions" || fail "Could not move result to $3_grouped_predictions"
mv -f "${TMP_PATH}/grouped_predictions.index" "$3_grouped_predictions.index" || fail "Could not move result to $3_grouped_predictions.index"

# create a subdb of the dp files:
if notExists "$3_dp_protein_contig_strand_map"; then
    "$MMSEQS" createsubdb "$3_grouped_predictions" "${INPUT_MAP}" "$3_dp_protein_contig_strand_map" \
        || fail "createsubdb on INPUT_MAP step died"
fi

if notExists "$3_dp_optimal_exon_sets"; then
    "$MMSEQS" createsubdb "$3_grouped_predictions" "${INPUT_OPTIMAL_EXON_SETS}" "$3_dp_optimal_exon_sets" \
        || fail "createsubdb on INPUT_OPTIMAL_EXON_SETS step died"
fi


if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files from ${TMP_PATH}"
    rm -f "${TMP_PATH}"/dp_*
fi
