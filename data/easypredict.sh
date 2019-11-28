#!/bin/sh -e

# easy predict workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <i:contigs> <i:targets> <o:predictionsFasta> <tmpDir>" && exit 1;
[   -f "$3" ] && echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

USER_INPUT_CONTIGS="$1"
USER_INPUT_TARGETS="$2"
TMP_PATH="$4"

INPUT_CONTIGS=""
INPUT_TARGETS=""

## handle contigs
# if input is a fasta
if notExists "${USER_INPUT_CONTIGS}.dbtype"; then
    # shellcheck disable=SC2086
    if notExists "${USER_INPUT_CONTIGS}"; then
        # shellcheck disable=SC2086
        echo "neither ${USER_INPUT_CONTIGS} nor ${USER_INPUT_CONTIGS}.dbtype was found!" && exit 1;
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${USER_INPUT_CONTIGS}" "${TMP_PATH}/contigs" --dbtype 2 ${VERBOSITY_COMP_PAR} \
            || fail "contigs createdb died"
    INPUT_CONTIGS="${TMP_PATH}/contigs"
else
    # db exists!
    INPUT_CONTIGS="${USER_INPUT_CONTIGS}"
fi


## handle targets
# if input is a fasta
if notExists "${USER_INPUT_TARGETS}.dbtype"; then
    if notExists "${USER_INPUT_TARGETS}"; then
        echo "neither ${USER_INPUT_TARGETS} nor ${USER_INPUT_TARGETS}.dbtype was found!" && exit 1;
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${USER_INPUT_TARGETS}" "${TMP_PATH}/targets" --dbtype 1 ${VERBOSITY_COMP_PAR} \
            || fail "targets createdb died"
    INPUT_TARGETS="${TMP_PATH}/targets"
else
    # db exists!
    INPUT_TARGETS="${USER_INPUT_TARGETS}"
fi

# produce MetaEuk calls by predictexons
if notExists "${TMP_PATH}/MetaEuk_calls.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" predictexons "${INPUT_CONTIGS}" "${INPUT_TARGETS}" "${TMP_PATH}/MetaEuk_calls" "${TMP_PATH}/tmp_predict" ${PREDICTEXONS_PAR} \
        || fail "predictexons step died"
fi

# reduce redundancy
if notExists "${TMP_PATH}/MetaEuk_preds.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" reduceredundancy "${TMP_PATH}/MetaEuk_calls" "${TMP_PATH}/MetaEuk_preds" "${TMP_PATH}/MetaEuk_preds_clust" ${REDUCEREDUNDANCY_PAR} \
        || fail "reduceredundancy step died"
fi

# create fasta
if notExists "$3"; then
    # shellcheck disable=SC2086
    "$MMSEQS" unitesetstofasta "${INPUT_CONTIGS}" "${INPUT_TARGETS}" "${TMP_PATH}/MetaEuk_preds" "$3" ${UNITESETSTOFASTA_PAR} \
        || fail "unitesetstofasta step died"
fi

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files from ${TMP_PATH}"
    "$MMSEQS" rmdb "${TMP_PATH}/contigs"
    "$MMSEQS" rmdb "${TMP_PATH}/contigs_h"
    "$MMSEQS" rmdb "${TMP_PATH}/targets"
    "$MMSEQS" rmdb "${TMP_PATH}/targets_h"
    "$MMSEQS" rmdb "${TMP_PATH}/MetaEuk_calls"
    "$MMSEQS" rmdb "${TMP_PATH}/MetaEuk_preds"
    "$MMSEQS" rmdb "${TMP_PATH}/MetaEuk_preds_clust"
    rm -r "${TMP_PATH}/tmp_predict"
    rm -f "${TMP_PATH}/easypredict.sh"
fi

