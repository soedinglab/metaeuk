#!/bin/bash -ex

METAEUK="$1"
DATAPATH="$2"
RESULTPATH="$3"

mkdir -p "${RESULTPATH}"

"${METAEUK}" easy-predict "${DATAPATH}/contigs.fna" "${DATAPATH}/proteins.faa" "${RESULTPATH}/predRedNoOver.fas" tempFolder --overlap 0 --metaeuk-tcov 0 --protein 1
"${METAEUK}" easy-predict "${DATAPATH}/contigs.fna" "${DATAPATH}/proteins.faa" "${RESULTPATH}/predRedOverAllowed.fas" tempFolder --overlap 1 --metaeuk-tcov 0 --protein 1

# check output #
perl compare_fasta_results.pl "${RESULTPATH}/predRedOverAllowed.fas" "${RESULTPATH}/predRedNoOver.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep_no_overlap.fas"