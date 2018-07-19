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

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <sequenceDB> <proteinTargetsDB> <outMetaeukBaseName> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p $4;

INPUT_CONTIGS="$(abspath "$1")"
INPUT_TARGET_PROTEINS="$(abspath "$2")"
TMP_PATH="$(abspath "$4")"

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
    $MMSEQS search "${TMP_PATH}/aa_6f" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/tmp_search" --alignment-mode 2\
        || fail "search step died"
fi

# swap results (result has targets as implicit keys)
if notExists "${TMP_PATH}/search_res_swap"; then
    $MMSEQS swapresults "${TMP_PATH}/aa_6f" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/search_res_swap" -e 0.1\
        || fail "swap step died"
fi

# join contig information to swapped results (result has additional info about the origin of the AA fragments)
if notExists "${TMP_PATH}/search_res_swap_w_contig_info"; then
    $MMSEQS filterdb "${TMP_PATH}/search_res_swap" "${TMP_PATH}/search_res_swap_w_contig_info" --join-db "${TMP_PATH}/nucl_6f_orf_lookup" --filter-column 1 \
        || fail "filterdb (to join contig info) step died"
fi

# sort joined swapped results by contig id
if notExists "${TMP_PATH}/search_res_swap_w_contig_info_sorted"; then
    $MMSEQS filterdb "${TMP_PATH}/search_res_swap_w_contig_info" "${TMP_PATH}/search_res_swap_w_contig_info_sorted" --sort-entries 1 --filter-column 11 \
        || fail "filterdb (to sort by contig) step died"
fi

# for each target, with respect to each contig and each strand, find the optimal set of exons
if notExists "${TMP_PATH}/dp_protein_contig_strand_map"; then
    $MMSEQS collectoptimalset "${TMP_PATH}/search_res_swap_w_contig_info_sorted" "${TMP_PATH}/dp" \
        || fail "collectoptimalset step died"
fi

# create sequence DBs from the results, these contain the original identifiers
if notExists "${TMP_PATH}/united_exons"; then
    $MMSEQS unitesetstosequencedb "${INPUT_CONTIGS}" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/dp" "${TMP_PATH}/united_exons" \
        || fail "unitesetstosequencedb step died"
fi

# post processing
mv -f "${TMP_PATH}/dp_protein_contig_strand_map" "$3_dp_protein_contig_strand_map" || fail "Could not move result to $3_dp_protein_contig_strand_map"
mv -f "${TMP_PATH}/dp_protein_contig_strand_map.index" "$3_dp_protein_contig_strand_map.index" || fail "Could not move result to $3_dp_protein_contig_strand_map.index"
mv -f "${TMP_PATH}/dp_optimal_exon_sets" "$3_dp_optimal_exon_sets" || fail "Could not move result to $3_dp_optimal_exon_sets"
mv -f "${TMP_PATH}/dp_optimal_exon_sets.index" "$3_dp_optimal_exon_sets.index" || fail "Could not move result to $3_dp_optimal_exon_sets.index"
mv -f "${TMP_PATH}/united_exons" "$3_united_exons" || fail "Could not move result to $3_united_exons"
mv -f "${TMP_PATH}/united_exons.index" "$3_united_exons.index" || fail "Could not move result to $3_united_exons.index"
mv -f "${TMP_PATH}/united_exons_h" "$3_united_exons_h" || fail "Could not move result to $3_united_exons_h"
mv -f "${TMP_PATH}/united_exons_h.index" "$3_united_exons_h.index" || fail "Could not move result to $3_united_exons_h.index"

# translate sequence DBs to AAs
if notExists "$3_united_exons_aa"; then
    $MMSEQS translatenucs "$3_united_exons" "$3_united_exons_aa" \
        || fail "translatenucs step died"
fi

# create a symbolic link for the map to the header and index file
ln -sf "$3_united_exons_aa_h" "$3_dp_protein_contig_strand_map_h" || fail "Could not create symbolic link for map headers"
ln -sf "$3_united_exons_aa_h.index" "$3_dp_protein_contig_strand_map_h.index" || fail "Could not create symbolic link for map headers index"


# if [ -n "$REMOVE_TMP" ]; then
#     echo "Removing temporary files"
#     rm -f "${TMP_PATH}/nucl_*"
#     rm -f "${TMP_PATH}/aa_*"
#     rm -f "${TMP_PATH}/search*"
#     rm -r "${TMP_PATH}/tmp_search/"
# fi
