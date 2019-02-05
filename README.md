# metaeuk - fast gene discovery in large-scale eukaryotic metagenomics data

[ ![Codeship Status for soedinglab/metaeuk](https://app.codeship.com/projects/07a9f310-7bb9-0136-3e65-3e3f6cc64c07/status?branch=master)](https://app.codeship.com/projects/300789)

Metaeuk is a modular toolkit designed for fast and large-scale gene calling and annotation in eukaryotic metagenomics contigs. Metaeuk combines the fast and sensitive homology search capabilities of [MMseqs2](https://github.com/soedinglab/MMseqs2) with a dynamic programming procedure to recover optimal exons sets. It reduces redundancies in multiple discoveries of the same gene and resolves conflicting gene predictions on the same strand. Metaeuk is GPL-licensed open source software that is implemented in C++ and available for Linux and macOS. The software is designed to run on multiple cores. Metaeuk was used to create a catalog of protein coding genes of marine eukaryotes based on Tara Oceans 912 metagenomics samples.

[Levy Larin E, Mirdita M and Soeding J. Metaeuk - fast processing of large-scale eukaryotic metagenomics data. biorxiv, doi: doi.org/????/??? (2019)](https://www.biorxiv.org/content/???/2019/??/??/?????).

## Input 
Metaeuk will search for eukaryotic protein-coding genes in **contigs** based on similarity to a reference database of **proteins** or **protein profiles**. The starting point are Fasta files of sequences.

Read [here](https://github.com/soedinglab/mmseqs2/wiki#how-to-create-a-target-profile-database-from-pfam) to learn more on how to create a protein profile database using MMseqs2. Once created, this database can be used as referenceDB in the command below.

## Running metaeuk 
### Main Modules:

      predictexons      	Predict eukaryotic exons based on protein similarity
      reduceredundancy  	A greedy approach to group metaeuk predictions which share an exon
      unitetoseqdbs     	Unite the exons of predictions to sequence BD


### Converting Fasta files to databases:
    
    # create a database from contigs.fasta 
      metaeuk createdb contigs.fasta contigsDB --dont-split-seq-by-len --dbtype 2

    # create a database from reference.fasta  
      metaeuk createdb reference.fasta referenceDB --dbtype 1


### Important parameters: 

     --min-length        minimal number of codons in putative protein fragment
     -e                  maximal E-Value to retain a match between a putative protein fragment and a reference taraget 
     --metaeuk-eval      maximal combined E-Value to retain an optimal exon set
     --e-profile         if refernceDB is a profile database, maximal E-Value between a putative fragment and a reference profile
     --slice-search      if refernceDB is a profile database, should be added
     

### Collecting optimal exons sets:

This module will extract all putative protein fragments from each contig (**C**) and strand (**S**), query them against the reference targets (**T**) and use dynamic programming to retain for each **T** the optimal compatible exon set from each **C** & **S** (thus creating **TCS** predictions).
    
    metaeuk predictexons contigsDB referenceDB predExResult tempFolder --metaeuk-eval 0.0001 -e 100 --min-length 40
    
Since this step involves a search, it is the most time-demanding of all analyses steps. Upon completion, it will output: predExResult_dp_protein_contig_strand_map predExResult_dp_optimal_exon_sets. These contain information about each of the **TCS** prediction and its exons.


### Reducing redundancy:

If there are homologies in referenceDB (e.g., T1 is highly similar to T2), the same optimal exons set from a **C** & **S** combination will be predicted more than once. This module will group together **TCS** predictions that share and exon and will choose their representative. In addition, it will greedily obtain a subset of the **TCS** representatives, such that there is no overlap of predictions on the same contig and strand.
    
    metaeuk reduceredundancy predExResult_dp_protein_contig_strand_map predExResult_dp_optimal_exon_sets redRedResult tempFolder
    
Upon completion, it will output: redRedResult_dp_protein_contig_strand_map redRedResult_dp_optimal_exon_sets, and  redRedResult_grouped_predictions. The first two contain information about the **TCS** representatives. The third, maps from the representative to all **TCS** predictions that share an exon with it.
In addition, it will output: redRedResult_no_overlap_dp_protein_contig_strand_map redRedResult_no_overlap_dp_optimal_exon_sets, and  redRedResult_grouped_predictions_no_overlap. The first two contain information about the **TCS** representatives after resolving overlaps. The third, maps from the representative to all **TCS** representatives that overlap it (excluded).


### Creating sequence DBs:

The X_dp_protein_contig_strand_map and X_dp_optimal_exon_sets produced by the modules above, can be used to extract the sequences of the predicted protein-coding genes.
    
    metaeuk unitetoseqdbs contigsDB referenceDB X_dp_protein_contig_strand_map X_no_overlap_dp_optimal_exon_sets result tempFolder
    
It will output: result_united_exons, result_united_exons_aa, which correspond to the codon (DNA) and amino-acid sequences of the predictions.


### Converting to Fasta:

The X_dp_protein_contig_strand_map and X_dp_optimal_exon_sets produced by the modules above, can be used to extract the sequences of the predicted protein-coding genes.
    
    metaeuk convert2fasta result_united_exons_aa result_united_exons_aa.fas
    
#### The metaeuk header:

The header is composed of several sections, separated by pipes ('|'):

*>T_header|C_header|S|bitscore|E-Value|number_exons|low_coord|high_coord|exon1_coords|exon2_coords|...*

*coord* refers to the coordination on the contig. It is advisable to keep T_header and C_header short and without pipes. The exon_coords are of the structure:
*low[taken_low]:high[taken_high]:nucleotide_length[taken_nucleotide_length]*
Since metaeuk allows for a very overlap on T of two putative exons (see P1 and P2 in the illustartion below), when joining the sequences of the exons, one of them is shortened.

<p align="center"><img src="https://github.com/soedinglab/metaeuk/blob/master/imgs/small_overlap_allowed.png" height="150"/></p>

Example header (two exons on the minus strand):
*>ERR1719262_736507|ERR868377_k119_2347399|-|341|6.2e-93|2|54324|54855|54855[54855]:54754[54754]:102[102]|54656[54668]:54324[54324]:333[321]*



## Compile from source
Compiling metaeuk from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile PLASS `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the PLASS binary will be located in the `build/bin` directory.

      git clone git@github.com:soedinglab/metaeuk.git .
      git submodule init
      git submodule update
      mkdir build
      cd build
      cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_MPI=1 -DCMAKE_INSTALL_PREFIX=. ..
      make -j
      make install
      export PATH="$(pwd)/bin/:$PATH"
        
:exclamation: If you want to compile metaeuk on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and metaeuk will not be able to run multithreaded. Use the following cmake call:

      CXX="$(brew --prefix)/bin/g++-8" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

## Hardware requirements
Metaeuk needs roughly 1 byte of memory per residue to work efficiently. Metaeuk will scale its memory consumption based on the available main memory of the machine. Metaeuk needs a CPU with at least the SSE4.1 instruction set to run. 


