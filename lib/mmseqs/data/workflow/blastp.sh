#!/bin/sh -e
# Sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    # create link to data file which contains a list of all targets that should be aligned
    ln -s "$(abspath "${TDB}.index")" "${RES}"
    # create new index repeatedly pointing to same entry
    INDEX_SIZE="$(wc -c < "${TDB}.index")"
    awk -v size="$INDEX_SIZE" '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    # create dbtype (7)
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

notExists() {
	[ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUT="$1"
TARGET="$2"
TMP_PATH="$4"

STEP=0
STEPS="${STEPS:-1}"
ALN_RES_MERGE="$TMP_PATH/aln_0"
while [ "$STEP" -lt "$STEPS" ]; do
    SENS_PARAM=SENSE_${STEP}
    eval SENS="\$$SENS_PARAM"

    # 1. Prefilter hits
    if notExists "$TMP_PATH/pref_$STEP.dbtype"; then
      if [ "$PREFMODE" = "EXHAUSTIVE" ]; then
          fake_pref "${INPUT}" "${TARGET}" "$TMP_PATH/pref_$STEP"
      elif [ "$PREFMODE" = "UNGAPPED" ]; then
          # shellcheck disable=SC2086
          $RUNNER "$MMSEQS" ungappedprefilter "$INPUT" "$TARGET" "$TMP_PATH/pref_$STEP" $UNGAPPEDPREFILTER_PAR \
              || fail "Ungapped prefilter died"
      else
          # shellcheck disable=SC2086
          $RUNNER "$MMSEQS" prefilter "$INPUT" "$TARGET" "$TMP_PATH/pref_$STEP" $PREFILTER_PAR -s "$SENS" \
              || fail "Prefilter died"
      fi
    fi

    # 2. alignment module
    if [ "$STEPS" -eq 1 ]; then
        if notExists "$3.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$TARGET${ALIGNMENT_DB_EXT}" "$TMP_PATH/pref_$STEP" "$3" $ALIGNMENT_PAR  \
                || fail "Alignment died"
        fi
        break
    else
        if notExists "$TMP_PATH/aln_$STEP.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$TARGET${ALIGNMENT_DB_EXT}" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_$STEP" $ALIGNMENT_PAR  \
                || fail "Alignment died"
        fi
    fi

    # only merge results after first step
    if [ "$STEP" -gt 0 ]; then
        if notExists "$TMP_PATH/aln_${SENS}.hasmerged"; then
            if [ "$STEP" -lt $((STEPS-1)) ]; then
                # shellcheck disable=SC2086
                "$MMSEQS" mergedbs "$1" "$TMP_PATH/aln_merge_new" "$ALN_RES_MERGE" "$TMP_PATH/aln_$STEP" ${VERB_COMP_PAR} \
                    || fail "Mergedbs died"
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "$TMP_PATH/aln_merge" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" mvdb "$TMP_PATH/aln_merge_new" "$TMP_PATH/aln_merge" ${VERBOSITY}
            else
                # shellcheck disable=SC2086
                "$MMSEQS" mergedbs "$1" "$3" "$ALN_RES_MERGE" "$TMP_PATH/aln_$STEP" ${VERB_COMP_PAR} \
                    || fail "Mergedbs died"
                break
            fi
            touch "$TMP_PATH/aln_${STEP}.hasmerged"
        fi
    fi
    if [ "$STEP" -gt 0 ]; then
      ALN_RES_MERGE="$TMP_PATH/aln_merge"
    fi

    NEXTINPUT="$TMP_PATH/input_$STEP"
    #do not create subdb at last step
    if [ "$STEP" -lt "$((STEPS-1))" ]; then
        if notExists "$TMP_PATH/order_$STEP.dbtype"; then
            awk '$3 < 2 { print $1 }' "$TMP_PATH/aln_$STEP.index" > "$TMP_PATH/order_$STEP" \
                || fail "Awk step $STEP died"
        fi

        if [ ! -s "$TMP_PATH/order_$STEP" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" mvdb "$ALN_RES_MERGE" "$3" ${VERBOSITY}
            break
        fi

        if notExists "$NEXTINPUT.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "$TMP_PATH/order_$STEP" "$INPUT" "$NEXTINPUT" ${VERBOSITY} --subdb-mode 1 \
                || fail "Order step $STEP died"
        fi
    fi
    INPUT="$NEXTINPUT"
    STEP="$((STEP+1))"
done

if [ -n "$REMOVE_TMP" ]; then
    STEP=0
    while [ "$STEP" -lt "$STEPS" ]; do
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_$STEP" ${VERBOSITY}
        rm -f "${TMP_PATH}/order_$STEP"
        STEP="$((STEP+1))"
    done
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_merge" ${VERBOSITY}
    rm -f "$TMP_PATH/blastp.sh"
fi


