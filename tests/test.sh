#!/bin/bash -ex

METAEUK="$1"
DATAPATH="$2"
RESULTPATH="$3"

SETGAPEXTEND="$4"
MINTARGETCOV="$5"

mkdir -p "${RESULTPATH}"
"${METAEUK}" createdb "${DATAPATH}/contigs.fna" "${RESULTPATH}/contigs" --dont-split-seq-by-len --dbtype 2
"${METAEUK}" createdb "${DATAPATH}/proteins.faa" "${RESULTPATH}/proteins" --dbtype 1

"${METAEUK}" predictexons "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predEx" "${RESULTPATH}/tempFolder" --set-gap-extend "${SETGAPEXTEND}" --metaeuk-tcov "${MINTARGETCOV}"
"${METAEUK}" unitesetstofasta "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predEx" "${RESULTPATH}/predEx.fas" --protein 1
"${METAEUK}" reduceredundancy "${RESULTPATH}/predEx" "${RESULTPATH}/predRed" "${RESULTPATH}/predClust"
"${METAEUK}" unitesetstofasta "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predRed" "${RESULTPATH}/predRed.fas" --protein 1
"${METAEUK}" groupstoacc "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predClust" "${RESULTPATH}/predClust.tsv"

# check ungrouped and grouped output
perl compare_fasta_results.pl "${RESULTPATH}/predEx.fas" "${RESULTPATH}/predRed.fas" "${DATAPATH}/as_should_final_united_exons_aa.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep.fas"