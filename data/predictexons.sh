#!/bin/bash -e

# predict exons workflow script
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
[ "$#" -ne 4 ] && echo "Please provide <sequenceDB> <proteinTargetsDB> <outDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p $4;

INPUT_CONTIGS="$(abspath $1)"
INPUT_TARGET_PROTEINS="$(abspath $2)"
TMP_PATH="$(abspath $4)"

# extract coding fragments from input contigs (result in DNA)
if notExists "${TMP_PATH}/nucl_6f"; then
    $MMSEQS extractorfs "${INPUT_CONTIGS}" "${TMP_PATH}/nucl_6f" --orf-start-mode 1 \
        || fail "extractorfs step died"
fi

# translate each coding fragment (result in AA)
if notExists "${TMP_PATH}/aa_6f"; then
    $MMSEQS translatenucs "${TMP_PATH}/nucl_6f" "${TMP_PATH}/aa_6f" \
        || fail "translatenucs step died"
fi

# search with each aa fragment against a target DB (result has queries as implicit keys)
if notExists "${TMP_PATH}/search_res"; then
    $MMSEQS search "${TMP_PATH}/aa_6f" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/tmp_search" -a \
        || fail "search step died"
fi

# swap results (result has target as implicit keys)
if notExists "${TMP_PATH}/search_res_swap"; then
    $MMSEQS swapresults "${TMP_PATH}/aa_6f" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/search_res_swap" \
        || fail "swap step died"
fi

# join contig information to swapped results (result has additional info about the origin of the AA fragments)
if notExists "${TMP_PATH}/search_res_swap_w_contig_info"; then
    $MMSEQS filterdb "${TMP_PATH}/search_res_swap" "${TMP_PATH}/search_res_swap_w_contig_info" --join-db "${TMP_PATH}/nucl_6f_orf_lookup" --filter-column 1 \
        || fail "filterdb (to join contig info) step died"
fi

# sort joined swapped results by contig id
if notExists "${TMP_PATH}/search_res_swap_w_contig_info_sorted"; then
    $MMSEQS filterdb "${TMP_PATH}/search_res_swap_w_contig_info" "${TMP_PATH}/search_res_swap_w_contig_info_sorted" --sort-entries 1 --filter-column 12 \
        || fail "filterdb (to sort by contig) step died"
fi

# # for each target, with respect to each contig and each strand, find the optimal set of exons
# if notExists "${TMP_PATH}/optimal_set_for_each_target"; then
#     $MMSEQS collectoptimalset "${TMP_PATH}/search_res_swap_w_contig_info_sorted" "${TMP_PATH}/optimal_set_for_each_target" \
#         || fail "collectoptimalset step died"
# fi

# post processing (WILL BE CHANGED AFTER WORKFLOW IS COMPLETE)
mv -f "${TMP_PATH}/search_res_swap_w_contig_info_sorted" "$3" || fail "Could not move result to $3"
mv -f "${TMP_PATH}/search_res_swap_w_contig_info_sorted.index" "$3.index" || fail "Could not move result to $3.index"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/nucl_*"
    rm -f "${TMP_PATH}/aa_*"
    rm -f "${TMP_PATH}/search*"
    rm -r "${TMP_PATH}/tmp_search/"
fi
