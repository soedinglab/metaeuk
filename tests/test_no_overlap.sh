#!/bin/bash -ex

METAEUK="$1"
DATAPATH="$2"
RESULTPATH="$3"

mkdir -p "${RESULTPATH}"
"${METAEUK}" createdb "${DATAPATH}/contigs.fna" "${RESULTPATH}/contigs" --dont-split-seq-by-len --dbtype 2
"${METAEUK}" createdb "${DATAPATH}/proteins.faa" "${RESULTPATH}/proteins" --dbtype 1

"${METAEUK}" predictexons "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predEx" tempFolder --metaeuk-eval 0.0001 -e 100 --min-length 20 --metaeuk-tcov 0
"${METAEUK}" reduceredundancy "${RESULTPATH}/predEx" "${RESULTPATH}/predRedOverAllowed" "${RESULTPATH}/predClust" --overlap 1
"${METAEUK}" unitesetstofasta "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predRedOverAllowed" "${RESULTPATH}/predRedOverAllowed.fas" --protein 1
"${METAEUK}" reduceredundancy "${RESULTPATH}/predEx" "${RESULTPATH}/predRedNoOver" "${RESULTPATH}/predClust" --overlap 0
"${METAEUK}" unitesetstofasta "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predRedNoOver" "${RESULTPATH}/predRedNoOver.fas" --protein 1
"${METAEUK}" groupstoacc "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predClust" "${RESULTPATH}/predClust.tsv"

# check output #
perl compare_fasta_results.pl "${RESULTPATH}/predRedOverAllowed.fas" "${RESULTPATH}/predRedNoOver.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep_no_overlap.fas"