#!/bin/bash -ex

METAEUK="$1"
DATAPATH="$2"
RESULTPATH="$3"

mkdir -p "${RESULTPATH}"
"${METAEUK}" createdb "${DATAPATH}/contigs.fna" "${RESULTPATH}/contigs" --dont-split-seq-by-len
"${METAEUK}" createdb "${DATAPATH}/proteins.faa" "${RESULTPATH}/proteins"

"${METAEUK}" predictexons "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predEx" "${RESULTPATH}/tempFolder" --metaeuk-eval 0.0001 -e 100 --min-length 20
"${METAEUK}" unitetoseqdbs "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/predEx_dp_protein_contig_strand_map" "${RESULTPATH}/predEx_dp_optimal_exon_sets" "${RESULTPATH}/final" "${RESULTPATH}/temp2"
"${METAEUK}" convert2fasta "${RESULTPATH}/final_united_exons_aa" "${RESULTPATH}/final_united_exons_aa.fas"

"${METAEUK}" reduceredundancy "${RESULTPATH}/predEx_dp_protein_contig_strand_map" "${RESULTPATH}/predEx_dp_optimal_exon_sets" "${RESULTPATH}/redRed" "${RESULTPATH}/temp3"
"${METAEUK}" unitetoseqdbs "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/redRed_dp_protein_contig_strand_map" "${RESULTPATH}/redRed_dp_optimal_exon_sets" "${RESULTPATH}/final_grouped" "${RESULTPATH}/temp4"
"${METAEUK}" convert2fasta "${RESULTPATH}/final_grouped_united_exons_aa" "${RESULTPATH}/final_grouped_predictions_rep.fas"

"${METAEUK}" unitetoseqdbs "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/redRed_no_overlap_dp_protein_contig_strand_map" "${RESULTPATH}/redRed_no_overlap_dp_optimal_exon_sets" "${RESULTPATH}/final_grouped_no_overlap" "${RESULTPATH}/temp4"
"${METAEUK}" convert2fasta "${RESULTPATH}/final_grouped_no_overlap_united_exons_aa" "${RESULTPATH}/final_grouped_predictions_rep_no_overlap.fas"

"${METAEUK}" createtsv "${RESULTPATH}/final_united_exons_aa" "${RESULTPATH}/final_united_exons_aa" "${RESULTPATH}/redRed_grouped_predictions" "${RESULTPATH}/redRed_grouped_predictions.tsv"



# check output (grouped with / without per-strand overlaps)
perl compare_fasta_results.pl "${RESULTPATH}/final_grouped_predictions_rep.fas" "${RESULTPATH}/final_grouped_predictions_rep_no_overlap.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep_no_overlap.fas"