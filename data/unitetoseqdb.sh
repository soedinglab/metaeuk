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
[ "$#" -ne 4 ] && echo "Please provide <contigsDB> <proteinTargetsDB> <metaeukBaseName> <tmpDir>" && exit 1;
# check if file exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;

INPUT_CONTIGS="$(abspath "$1")"
INPUT_TARGET_PROTEINS="$(abspath "$2")"
INPUT_MAP="$(abspath "$3_dp_protein_contig_strand_map")"
INPUT_OPTIMAL_EXON_SETS="$(abspath "$3_dp_optimal_exon_sets")"

[ ! -f "${INPUT_MAP}" ] &&  echo "${INPUT_MAP} not found!" && exit 1;
[ ! -f "${INPUT_OPTIMAL_EXON_SETS}" ] &&  echo "${INPUT_OPTIMAL_EXON_SETS} not found!" && exit 1;

TMP_PATH="$(abspath "$4")"
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";


# if only on representative sequence:
UNITED_EXONS_NAME="$3_united_exons"
if [ -n "$UNITE_EXONS_REP" ]; then
    # shellcheck disable=SC2086
    UNITED_EXONS_NAME="$3_united_exons_rep"
fi

if notExists "${UNITED_EXONS_NAME}"; then
    # shellcheck disable=SC2086
    if [ -n "$UNITE_EXONS_REP" ]; then
    # shellcheck disable=SC2086
        INPUT_GROUPED_PREDS="$(abspath "$3_grouped_predictions")"
        [ ! -f "${INPUT_GROUPED_PREDS}" ] &&  echo "${INPUT_GROUPED_PREDS} not found!" && exit 1;
        
        # create a subdb of the dp files:
        "$MMSEQS" createsubdb "${INPUT_GROUPED_PREDS}" "${INPUT_MAP}" "$3_dp_protein_contig_strand_map_rep" \
            || fail "createsubdb on map step died"

        "$MMSEQS" createsubdb "${INPUT_GROUPED_PREDS}" "${INPUT_OPTIMAL_EXON_SETS}" "$3_dp_optimal_exon_sets_rep" \
                || fail "createsubdb on exon sets step died"

        "$MMSEQS" unitesetstosequencedb "${INPUT_CONTIGS}" "${INPUT_TARGET_PROTEINS}" "$3_dp_protein_contig_strand_map_rep" "$3_dp_optimal_exon_sets_rep" "${UNITED_EXONS_NAME}" ${THREADS_PAR} \
            || fail "unitesetstosequencedb step died"
    else
    # here we work on all predictions (not reduced):
    # shellcheck disable=SC2086
        "$MMSEQS" unitesetstosequencedb "${INPUT_CONTIGS}" "${INPUT_TARGET_PROTEINS}" "${INPUT_MAP}" "${INPUT_OPTIMAL_EXON_SETS}" "${UNITED_EXONS_NAME}" ${THREADS_PAR} \
            || fail "unitesetstosequencedb step died"
    fi

    # the following part is joint to either procedure
    # translate sequence DBs to AAs
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
