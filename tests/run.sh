#!/bin/bash -ex

./test.sh "$1" minus_strand minus_strand_results
./test.sh "$1" multi_exon multi_exon_results
./test.sh "$1" two_contigs two_contigs_results
