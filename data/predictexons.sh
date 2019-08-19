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
[ "$#" -ne 4 ] && echo "Please provide <contigsDB> <targetDB> <outPredictionsDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUT_CONTIGS="$(abspath "$1")"
INPUT_TARGET_PROTEINS="$(abspath "$2")"
TMP_PATH="$(abspath "$4")"

# extract coding fragments from input contigs (result in DNA)
if notExists "${TMP_PATH}/nucl_6f.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT_CONTIGS}" "${TMP_PATH}/nucl_6f" ${EXTRACTORFS_PAR} \
        || fail "extractorfs step died"
fi

# translate each coding fragment (result in AA)
if notExists "${TMP_PATH}/aa_6f.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" translatenucs "${TMP_PATH}/nucl_6f" "${TMP_PATH}/aa_6f" ${TRANSLATENUCS_PAR} \
        || fail "translatenucs step died"
fi

# when running in null mode (to assess evalues), reverse the AA fragments:
AA_FRAGS="${TMP_PATH}/aa_6f"
if [ -n "$REVERSE_FRAGMENTS" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" reverseseq "${AA_FRAGS}" "${AA_FRAGS}_reverse" ${THREADS_PAR} \
        || fail "reverseseq step died"
    AA_FRAGS="${AA_FRAGS}_reverse"
    echo "Will base search on ${AA_FRAGS}"
fi

# search with each aa fragment against a target DB (result has queries as implicit keys)
if notExists "${TMP_PATH}/search_res.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${AA_FRAGS}" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/tmp_search" ${SEARCH_PAR} \
        || fail "search step died"
fi

# augment the search results with contig info and write a double alignment format where contigs are keys
if notExists "${TMP_PATH}/search_res_by_contig.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" resultspercontig "${INPUT_CONTIGS}" "${TMP_PATH}/nucl_6f" "${TMP_PATH}/search_res" "${TMP_PATH}/search_res_by_contig" ${THREADS_PAR} \
        || fail "resultspercontig step died"
fi

# for each target, with respect to each contig and each strand, find the optimal set of exons
if notExists "${TMP_PATH}/dp_predictions.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" collectoptimalset "${TMP_PATH}/search_res_by_contig" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/dp_predictions" ${COLLECTOPTIMALSET_PAR} \
        || fail "collectoptimalset step died"
fi

# post processing
"$MMSEQS" mvdb "${TMP_PATH}/dp_predictions" "$3" || fail "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files from ${TMP_PATH}"
    rm -f "${TMP_PATH}"/nucl_6f*
    rm -f "${TMP_PATH}"/aa_6f*
    rm -f "${TMP_PATH}"/search_res*
    rm -r "${TMP_PATH}"/tmp_search/
fi

