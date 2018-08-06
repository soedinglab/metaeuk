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
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d $(dirname "$1") ]; then
            echo "$(cd $(dirname "$1"); pwd)/$(basename "$1")"
    fi
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <metaeukBaseName> <outDB> <tmpDir>" && exit 1;
# check if file exists
INPUT_MAP="$(abspath "$1_dp_protein_contig_strand_map")"
INPUT_UNITED_EXONS="$(abspath "$1_united_exons_aa")"

[ ! -f "${INPUT_MAP}" ] &&  echo "${INPUT_MAP} not found!" && exit 1;
[ ! -f "${INPUT_UNITED_EXONS}" ] &&  echo "${INPUT_UNITED_EXONS} not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p $3;


TMP_PATH="$(abspath "$3")"

# aggregate all TCS (Target + Contig + Strand) predictions by their CS
if notExists "${TMP_PATH}/dp_contig_strand_map"; then
    $MMSEQS swapdb "${INPUT_MAP}" "${TMP_PATH}/dp_contig_strand_map" \
        || fail "swapdb step died"
fi

## the following three steps create a primary, secondary and tertiary order:

# within each CS, sort by bitscore (decreasing order), tertiary order
if notExists "${TMP_PATH}/dp_contig_strand_map_sorted_by_bitscore"; then
    $MMSEQS filterdb "${TMP_PATH}/dp_contig_strand_map" "${TMP_PATH}/dp_contig_strand_map_sorted_by_bitscore" --sort-entries 2 --filter-column 5 \
        || fail "filterdb (to sort by bitscore in decreasing order) step died"
fi

# within each CS, sort by number of exons (decreasing order), secondary order
if notExists "${TMP_PATH}/dp_contig_strand_map_sorted_by_num_exons_subsorted_by_bit"; then
    $MMSEQS filterdb "${TMP_PATH}/dp_contig_strand_map_sorted_by_bitscore" "${TMP_PATH}/dp_contig_strand_map_sorted_by_num_exons_subsorted_by_bit" --sort-entries 2 --filter-column 6 \
        || fail "filterdb (to sort by number of exons in decreasing order) step died"
fi

# within each CS, sort by start position on the contig (increasing order), primary order
if notExists "${TMP_PATH}/dp_contig_strand_map_sorted_by_start_subsorted_by_num_and_bit"; then
    $MMSEQS filterdb "${TMP_PATH}/dp_contig_strand_map_sorted_by_num_exons_subsorted_by_bit" "${TMP_PATH}/dp_contig_strand_map_sorted_by_start_subsorted_by_num_and_bit" --sort-entries 1 --filter-column 7 \
        || fail "filterdb (to sort by start position on the contig in increasing order) step died"
fi

# greedy grouping of predictions based on the sorted map
if notExists "${TMP_PATH}/grouped_predictions"; then
    $MMSEQS grouppredictions "${TMP_PATH}/dp_contig_strand_map_sorted_by_start_subsorted_by_num_and_bit" "${TMP_PATH}/grouped_predictions" \
        || fail "grouppredictions step died"
fi

mv -f "${TMP_PATH}/grouped_predictions" "$1_grouped_predictions" || fail "Could not move result to $1_grouped_predictions"
mv -f "${TMP_PATH}/grouped_predictions.index" "$1_grouped_predictions.index" || fail "Could not move result to $1_grouped_predictions.index"

# extracting a representative sequence from each group
if notExists "$1_grouped_predictions_rep"; then
    $MMSEQS result2repseq "${INPUT_UNITED_EXONS}" "$1_grouped_predictions" "$1_grouped_predictions_rep" \
        || fail "result2repseq step died"
fi

# create a symbolic link for the grouped predictions to the header and index file
ln -sf "${INPUT_UNITED_EXONS}_h" "$1_grouped_predictions_h" || fail "Could not create symbolic link for grouped predictions headers"
ln -sf "${INPUT_UNITED_EXONS}_h.index" "$1_grouped_predictions_h.index" || fail "Could not create symbolic link for grouped predictions headers index"
ln -sf "${INPUT_UNITED_EXONS}_h" "$1_grouped_predictions_rep_h" || fail "Could not create symbolic link for rep grouped predictions headers"
ln -sf "${INPUT_UNITED_EXONS}_h.index" "$1_grouped_predictions_rep_h.index" || fail "Could not create symbolic link for rep grouped predictions headers index"


# if [ -n "$REMOVE_TMP" ]; then
#     echo "Removing temporary files"
#     rm -f "${TMP_PATH}/nucl_*"
#     rm -f "${TMP_PATH}/aa_*"
#     rm -f "${TMP_PATH}/search*"
#     rm -r "${TMP_PATH}/tmp_search/"
# fi
