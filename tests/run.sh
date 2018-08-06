#!/bin/bash -e

git clone git@github.com:soedinglab/metaeuk.git
cd metaeuk
git submodule init
git submodule update
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make -j 8
cd ../..

./test.sh metaeuk/build/src/metaeuk minus_strand minus_strand_results
./test.sh metaeuk/build/src/metaeuk multi_exon multi_exon_results
./test.sh metaeuk/build/src/metaeuk two_contigs two_contigs_results
