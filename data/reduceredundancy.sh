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

# check amount of input variables
[ "$#" -ne 3 ] && echo "Please provide <metaeukBaseName> <outDB> <tmpDir>" && exit 1;
# check if file exists
INPUT_MAP="$(abspath "$1_dp_protein_contig_strand_map")"

[ ! -f "${INPUT_MAP}" ] &&  echo "${INPUT_MAP} not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p $3;


TMP_PATH="$(abspath "$3")"

# aggregate all TCS (Target + Contig + Strand) predictions by their CS
if notExists "${TMP_PATH}/$1_dp_contig_strand_map"; then
    $MMSEQS swapdb "${INPUT_MAP}" "${TMP_PATH}/$1_dp_contig_strand_map" \
        || fail "swapdb step died"
fi

## the following two steps create a primary and secondary order:

# within each CS, sort by number of exons (decreasing order), secondary order
if notExists "${TMP_PATH}/$1_dp_contig_strand_map_sorted_by_num_exons"; then
    $MMSEQS filterdb "${TMP_PATH}/$1_dp_contig_strand_map" "${TMP_PATH}/$1_dp_contig_strand_map_sorted_by_num_exons" --sort-entries 2 --filter-column 6 \
        || fail "filterdb (to sort by number of exons in decreasing order) step died"
fi

# within each CS, sort by start position on the contig (increasing order), primary order
if notExists "${TMP_PATH}/$1_dp_contig_strand_map_sorted_by_start_subsorted_by_num_exons"; then
    $MMSEQS filterdb "${TMP_PATH}/$1_dp_contig_strand_map_sorted_by_num_exons" "${TMP_PATH}/$1_dp_contig_strand_map_sorted_by_start_subsorted_by_num_exons" --sort-entries 1 --filter-column 7 \
        || fail "filterdb (to sort by start position on the contig in increasing order) step died"
fi

# greedy grouping of predictions based on the sorted map
if notExists "${TMP_PATH}/$1_grouped_predictions"; then
    $MMSEQS grouppredictions "${TMP_PATH}/$1_dp_contig_strand_map_sorted_by_start_subsorted_by_num_exons" "${TMP_PATH}/$1_grouped_predictions" \
        || fail "grouppredictions step died"
fi


mv -f "${TMP_PATH}/$1_grouped_predictions" "$1_grouped_predictions" || fail "Could not move result to $1_grouped_predictions"

# if [ -n "$REMOVE_TMP" ]; then
#     echo "Removing temporary files"
#     rm -f "${TMP_PATH}/nucl_*"
#     rm -f "${TMP_PATH}/aa_*"
#     rm -f "${TMP_PATH}/search*"
#     rm -r "${TMP_PATH}/tmp_search/"
# fi
