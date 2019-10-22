#!/bin/bash -ex

#if MPI;
#  export RUNNER="mpirun ..."


./test.sh "$1" minus_strand minus_strand_results -1 0
./test.sh "$1" multi_exon multi_exon_results -1 0
./test.sh "$1" two_contigs two_contigs_results -1 0
./test.sh "$1" target_overlap target_overlap_results -3 0
./test.sh "$1" cluster_rep cluster_rep_results -3 0
./test.sh "$1" target_cov target_cov_results -1 0.7

./test_no_overlap.sh "$1" no_overlap no_overlap_results

