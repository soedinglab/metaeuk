#!/bin/bash -ex

METAEUK="$1"
DATAPATH="$2"
RESULTPATH="$3"

mkdir -p "${RESULTPATH}"
"${METAEUK}" createdb "${DATAPATH}/contigs.fna" "${RESULTPATH}/contigs" --dont-split-seq-by-len --dbtype 2
"${METAEUK}" createdb "${DATAPATH}/proteins.faa" "${RESULTPATH}/proteins" --dbtype 1

"${METAEUK}" predictexons "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predEx" "${RESULTPATH}/tempFolder" --set-gap-open -1 --set-gap-extend -3
"${METAEUK}" unitetoseqdbs "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predEx_dp_protein_contig_strand_map" "${RESULTPATH}/predEx_dp_optimal_exon_sets" "${RESULTPATH}/final" "${RESULTPATH}/temp2"
"${METAEUK}" convert2fasta "${RESULTPATH}/final_united_exons_aa" "${RESULTPATH}/final_united_exons_aa.fas"

"${METAEUK}" reduceredundancy "${RESULTPATH}/predEx_dp_protein_contig_strand_map" "${RESULTPATH}/predEx_dp_optimal_exon_sets" "${RESULTPATH}/redRed" "${RESULTPATH}/temp3"
"${METAEUK}" unitetoseqdbs "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/redRed_dp_protein_contig_strand_map" "${RESULTPATH}/redRed_dp_optimal_exon_sets" "${RESULTPATH}/final_grouped" "${RESULTPATH}/temp4"
"${METAEUK}" convert2fasta "${RESULTPATH}/final_grouped_united_exons_aa" "${RESULTPATH}/final_grouped_predictions_rep.fas"

# check output
perl compare_fasta_results.pl "${RESULTPATH}/final_united_exons_aa.fas" "${RESULTPATH}/final_grouped_predictions_rep.fas" "${DATAPATH}/as_should_final_united_exons_aa.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep.fas"