#!/bin/bash -e

# assign workflow script
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
[ "$#" -ne 4 ] && echo "Please provide <metaeukBaseName> <uniprotFasta> <outDB> <tmpDir>" && exit 1;
# check if file exists
INPUT_MAP="$(abspath "$1_dp_protein_contig_strand_map")"
INPUT_UNITED_EXONS="$(abspath "$1_united_exons_aa")"
INPUT_GROUPED_PREDICTIONS_REP="$(abspath "$1_grouped_predictions_rep")"
INPUT_UNIPROT_REF_FASTA="$(abspath "$2")"

[ ! -f "${INPUT_MAP}" ] &&  echo "${INPUT_MAP} not found!" && exit 1;
[ ! -f "${INPUT_UNITED_EXONS}" ] &&  echo "${INPUT_UNITED_EXONS} not found!" && exit 1;
[ ! -f "${INPUT_GROUPED_PREDICTIONS_REP}" ] &&  echo "${INPUT_GROUPED_PREDICTIONS_REP} not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

TMP_PATH="$(abspath "$4")"

# prepare reference taxonomy DB
if notExists "${TMP_PATH}/uniprotDB"; then
    "$MMSEQS" createdb "${INPUT_UNIPROT_REF_FASTA}" "${TMP_PATH}/uniprotDB" \
        || fail "createdb step died"
fi

if notExists "${TMP_PATH}/uniprotDB"; then
    "$MMSEQS" createtaxdb "${TMP_PATH}/uniprotDB" "${TMP_PATH}" \
        || fail "createtaxdb step died"
fi


# map each contig C to its predictions TCS
if notExists "${TMP_PATH}/pred_to_contig"; then
    "$MMSEQS" filterdb "${INPUT_MAP}" "${TMP_PATH}/pred_to_contig" --trim-to-one-column --filter-column 3 \
        || fail "filterdb step died"
fi

if notExists "${TMP_PATH}/contig_to_pred"; then
    "$MMSEQS" swapdb "${TMP_PATH}/pred_to_contig" "${TMP_PATH}/contig_to_pred" \
        || fail "swapdb step died"
fi

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files from ${TMP_PATH}"
    rm -f "${TMP_PATH}"/pred_to_contig*
    rm -f "${TMP_PATH}"/contig_to_pred*
fi
