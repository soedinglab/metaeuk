#!/bin/sh -e

# tax to contig workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}
# check number of input variables
[ "$#" -ne 6 ] && echo "Please provide <i:contigsDb> <i:predictionsFasta> <i:predictionsFasta.headersMap.tsv> <i:taxAnnotTargetDb> <o:taxResult> <tmpDir>" && exit 1;
[   -f "$5" ] && echo "$5 exists already!" && exit 1;
[ ! -d "$6" ] && echo "tmp directory $6 not found!" && mkdir -p "$6";

CONTIGS_DB="$1"
PREDS_FASTA="$2"
PREDS_HEADERS_MAP="$3"
TAX_TARGET_DB="$4"
TAX_ASSIGNMENT_BASENAME="$5"
TMP_PATH="$6"

# convert fasta to db
if notExists "${TMP_PATH}/preds.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${PREDS_FASTA}" "${TMP_PATH}/preds" --shuffle 0 \
        || fail "createdb step died"
fi

# number tsv map lines
if notExists "${TMP_PATH}/preds_map_num.tsv"; then
    # shellcheck disable=SC2086
    awk '{print NR-1, "\t", $0}' "${PREDS_HEADERS_MAP}" > "${TMP_PATH}/preds_map_num.tsv" \
        || fail "awk step died"
fi

# create db from tsv map (pred -> contig)
if notExists "${TMP_PATH}/preds_map_num.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" tsv2db "${TMP_PATH}/preds_map_num.tsv" "${TMP_PATH}/preds_map_num" --output-dbtype 0 \
        || fail "tsv2db step died"
fi

# swap tsv map (contig -> pred)
if notExists "${TMP_PATH}/preds_map_num_swapped.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapdb "${TMP_PATH}/preds_map_num" "${TMP_PATH}/preds_map_num_swapped" \
        || fail "swapdb step died"
fi

# run taxonomy on the predictions
if notExists "${TMP_PATH}/tax_per_pred.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" taxonomy "${TMP_PATH}/preds" "${TAX_TARGET_DB}" "${TMP_PATH}/tax_per_pred" "${TMP_PATH}/tmp_taxonomy" ${TAXONOMY_PAR} \
        || fail "taxonomy died"
fi

# run aggregatetaxweights on the tax_per_pred and the mapping
if [ ! -e "${RESULTS}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" aggregatetaxweights "${TAX_TARGET_DB}" "${TMP_PATH}/preds_map_num_swapped" "${TMP_PATH}/tax_per_pred" "${TMP_PATH}/tax_per_pred_aln" "${TMP_PATH}/tax_per_contig" ${AGGREGATETAX_PAR} \
        || fail "aggregatetaxweights died"
fi

# create tsv for predictions
if notExists "${TAX_ASSIGNMENT_BASENAME}_tax_per_pred.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/preds" "${TMP_PATH}/tax_per_pred" "${TAX_ASSIGNMENT_BASENAME}_tax_per_pred.tsv" ${CREATETSV_PAR} \
        || fail "createtsv on predictions died"
fi

# create tsv for contigs
if notExists "${TAX_ASSIGNMENT_BASENAME}_tax_per_contig.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${CONTIGS_DB}" "${TMP_PATH}/tax_per_contig" "${TAX_ASSIGNMENT_BASENAME}_tax_per_contig.tsv" ${CREATETSV_PAR} \
        || fail "createtsv on contigs died"
fi

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files from ${TMP_PATH}"
    "$MMSEQS" rmdb "${TMP_PATH}/preds"
    "$MMSEQS" rmdb "${TMP_PATH}/preds_h"
    rm -f "${TMP_PATH}/preds_map_num.tsv"
    "$MMSEQS" rmdb "${TMP_PATH}/preds_map_num"
    "$MMSEQS" rmdb "${TMP_PATH}/preds_map_num_swapped"
    "$MMSEQS" rmdb "${TMP_PATH}/tax_per_pred_aln"
    rm -r "${TMP_PATH}/tmp_taxonomy"
    rm -f "${TMP_PATH}/taxtocontig.sh"
fi

