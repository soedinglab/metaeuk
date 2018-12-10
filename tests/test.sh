#!/bin/bash -ex

METAEUK="$1"
DATAPATH="$2"
RESULTPATH="$3"

mkdir -p "${RESULTPATH}"
"${METAEUK}" createdb "${DATAPATH}/contigs.fna" "${RESULTPATH}/contigs" --dont-split-seq-by-len
"${METAEUK}" createdb "${DATAPATH}/proteins.faa" "${RESULTPATH}/proteins"
"${METAEUK}" predictexons "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/final" "${RESULTPATH}/tempFolder"
"${METAEUK}" reduceredundancy "${RESULTPATH}/final" "${RESULTPATH}/final" "${RESULTPATH}/tempGroup"
"${METAEUK}" unitetoseqdb "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/final" "${RESULTPATH}/tempSeq"
"${METAEUK}" convert2fasta "${RESULTPATH}/final_united_exons_aa" "${RESULTPATH}/final_united_exons_aa.fas"
"${METAEUK}" unitetoseqdb "${RESULTPATH}/contigs" "${RESULTPATH}/proteins" "${RESULTPATH}/final" "${RESULTPATH}/tempSeq" --unite-exons-rep 1
"${METAEUK}" convert2fasta "${RESULTPATH}/final_united_exons_rep_aa" "${RESULTPATH}/final_grouped_predictions_rep.fas"
"${METAEUK}" createtsv "${RESULTPATH}/final_united_exons_aa" "${RESULTPATH}/final_united_exons_aa" "${RESULTPATH}/final_grouped_predictions" "${RESULTPATH}/final_grouped_predictions.tsv"

# check orf to contig alignment procedure #
perl check_orf_to_contig.pl "${RESULTPATH}/tempFolder" "${DATAPATH}/as_should_nucl_6f_orf_aligned_to_contig"

# check output
perl compare_fasta_results.pl "${RESULTPATH}/final_united_exons_aa.fas" "${RESULTPATH}/final_grouped_predictions_rep.fas" "${DATAPATH}/as_should_final_united_exons_aa.fas" "${DATAPATH}/as_should_final_grouped_predictions_rep.fas"