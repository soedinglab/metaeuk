#!/bin/bash -ex

METAEUK="$1"
DATAPATH="$2"
RESULTPATH="$3"

mkdir -p "${RESULTPATH}"

"${METAEUK}" easy-predict "${DATAPATH}/contigs.fna" "${DATAPATH}/proteins.faa" "${RESULTPATH}/predRedNoOver" tempFolder --overlap 0 --metaeuk-tcov 0
"${METAEUK}" easy-predict "${DATAPATH}/contigs.fna" "${DATAPATH}/proteins.faa" "${RESULTPATH}/predRedOverAllowed" tempFolder --overlap 1 --metaeuk-tcov 0

# check output #
perl compare_fasta_results.pl "${RESULTPATH}/predRedOverAllowed.fas" "${RESULTPATH}/predRedNoOver.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep_no_overlap.fas"