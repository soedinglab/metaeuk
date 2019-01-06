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
[ "$#" -ne 4 ] && echo "Please provide <contigsDB> <proteinsDB> <predictexonsBaseName> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUT_CONTIGS="$(abspath "$1")"
INPUT_TARGET_PROTEINS="$(abspath "$2")"
TMP_PATH="$(abspath "$4")"

# extract coding fragments from input contigs (result in DNA)
if notExists "${TMP_PATH}/nucl_6f"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT_CONTIGS}" "${TMP_PATH}/nucl_6f" ${EXTRACTORFS_PAR} \
        || fail "extractorfs step died"
fi

# write extracted orfs locations on contig in alignment format
if notExists "${TMP_PATH}/nucl_6f_orf_aligned_to_contig"; then
    # shellcheck disable=SC2086
    "$MMSEQS" orftocontig "${INPUT_CONTIGS}" "${TMP_PATH}/nucl_6f" "${TMP_PATH}/nucl_6f_orf_aligned_to_contig" ${THREADS_PAR} \
        || fail "orftocontig step died"
fi

# translate each coding fragment (result in AA)
if notExists "${TMP_PATH}/aa_6f"; then
    # shellcheck disable=SC2086
    "$MMSEQS" translatenucs "${TMP_PATH}/nucl_6f" "${TMP_PATH}/aa_6f" ${TRANSLATENUCS_PAR} \
        || fail "translatenucs step died"
fi

# when running in null mode (to assess evalues), we reverse the AA fragments:
AA_FRAGS="${TMP_PATH}/aa_6f"
if [ -n "$REVERSE_FRAGMENTS" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" reverseseq "${AA_FRAGS}" "${AA_FRAGS}_reverse" ${THREADS_PAR} \
        || fail "reverseseq step died"
    AA_FRAGS="${AA_FRAGS}_reverse"
    echo "Will base search on ${AA_FRAGS}"
fi

# search with each aa fragment against a target DB (result has queries as implicit keys)
if notExists "${TMP_PATH}/search_res"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${AA_FRAGS}" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/tmp_search" ${SEARCH_PAR} \
        || fail "search step died"
fi

# swap results (result has targets as implicit keys)
if notExists "${TMP_PATH}/search_res_swap"; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapresults "${AA_FRAGS}" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/search_res_swap" -e 1000 ${THREADS_PAR} \
        || fail "swap step died"
fi

# join contig information to swapped results (result has additional info about the origin of the AA fragments)
if notExists "${TMP_PATH}/search_res_swap_w_contig_info"; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/search_res_swap" "${TMP_PATH}/search_res_swap_w_contig_info" --join-db "${TMP_PATH}/nucl_6f_orf_aligned_to_contig" --filter-column 1 ${THREADS_PAR} \
        || fail "filterdb (to join contig info) step died"
fi

# sort joined swapped results by contig id
if notExists "${TMP_PATH}/search_res_swap_w_contig_info_sorted"; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/search_res_swap_w_contig_info" "${TMP_PATH}/search_res_swap_w_contig_info_sorted" --sort-entries 1 --filter-column 11 ${THREADS_PAR} \
        || fail "filterdb (to sort by contig) step died"
fi

# for each target, with respect to each contig and each strand, find the optimal set of exons
if notExists "${TMP_PATH}/dp_protein_contig_strand_map"; then
    # shellcheck disable=SC2086
    "$MMSEQS" collectoptimalset "${TMP_PATH}/search_res_swap_w_contig_info_sorted" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/" ${COLLECTOPTIMALSET_PAR} \
        || fail "collectoptimalset step died"
fi

# post processing
mv -f "${TMP_PATH}/dp_protein_contig_strand_map" "$3_dp_protein_contig_strand_map" || fail "Could not move result to $3_dp_protein_contig_strand_map"
mv -f "${TMP_PATH}/dp_protein_contig_strand_map.index" "$3_dp_protein_contig_strand_map.index" || fail "Could not move result to $3_dp_protein_contig_strand_map.index"
mv -f "${TMP_PATH}/dp_optimal_exon_sets" "$3_dp_optimal_exon_sets" || fail "Could not move result to $3_dp_optimal_exon_sets"
mv -f "${TMP_PATH}/dp_optimal_exon_sets.index" "$3_dp_optimal_exon_sets.index" || fail "Could not move result to $3_dp_optimal_exon_sets.index"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files from ${TMP_PATH}"
    rm -f "${TMP_PATH}"/nucl_6f*
    rm -f "${TMP_PATH}"/aa_6f*
    rm -f "${TMP_PATH}"/search_res*
    rm -r "${TMP_PATH}"/tmp_search/
fi

