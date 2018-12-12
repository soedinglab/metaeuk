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
[ "$#" -ne 6 ] && echo "Please provide <contigsDB> <proteinsDB> <dp_protein_contig_strand_map> <dp_optimal_exon_sets> <o:unitedexonsBasename> <tmpDir>" && exit 1;
# check if file exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[ ! -f "$3" ] &&  echo "$3 not found!" && exit 1;
[ ! -f "$4" ] &&  echo "$4 not found!" && exit 1;

INPUT_CONTIGS="$(abspath "$1")"
INPUT_TARGET_PROTEINS="$(abspath "$2")"
INPUT_MAP="$(abspath "$3")"
INPUT_OPTIMAL_EXON_SETS="$(abspath "$4")"

TMP_PATH="$(abspath "$6")"
[ ! -d "$6" ] &&  echo "tmp directory $6 not found!" && mkdir -p "$6";


UNITED_EXONS_NAME="$(abspath "$5_united_exons")"
if notExists "${UNITED_EXONS_NAME}"; then
    # shellcheck disable=SC2086
    "$MMSEQS" unitesetstosequencedb "${INPUT_CONTIGS}" "${INPUT_TARGET_PROTEINS}" "${INPUT_MAP}" "${INPUT_OPTIMAL_EXON_SETS}" "${UNITED_EXONS_NAME}" ${THREADS_PAR} \
            || fail "unitesetstosequencedb step died"

    # translate sequence DB to AAs
    if notExists "${UNITED_EXONS_NAME}_aa"; then
        # shellcheck disable=SC2086
        UNITED_EXONS_FILESIZE="$(stat -c%s "${UNITED_EXONS_NAME}")"
        if (( "${UNITED_EXONS_FILESIZE}" > 0 )); then
            # shellcheck disable=SC2086
            "$MMSEQS" translatenucs "${UNITED_EXONS_NAME}" "${UNITED_EXONS_NAME}_aa" ${TRANSLATENUCS_PAR} \
            || fail "translatenucs step died"
        fi
    fi
fi


if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files from ${TMP_PATH}"
fi
