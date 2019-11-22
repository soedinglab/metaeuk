# MetaEuk - sensitive, high-throughput gene discovery and annotation for large-scale eukaryotic metagenomics

[ ![Codeship Status for soedinglab/metaeuk](https://app.codeship.com/projects/07a9f310-7bb9-0136-3e65-3e3f6cc64c07/status?branch=master)](https://app.codeship.com/projects/300789)

MetaEuk is a modular toolkit designed for large-scale gene discovery and annotation in eukaryotic metagenomic contigs. Metaeuk combines the fast and sensitive homology search capabilities of [MMseqs2](https://github.com/soedinglab/MMseqs2) with a dynamic programming procedure to recover optimal exons sets. It reduces redundancies in multiple discoveries of the same gene and resolves conflicting gene predictions on the same strand. MetaEuk is GPL-licensed open source software that is implemented in C++ and available for Linux and macOS. The software is designed to run on multiple cores.

Levy Karin E, Mirdita M and Soeding J. MetaEuk – sensitive, high-throughput gene discovery and annotation for large-scale eukaryotic metagenomics. submitted, 2019.

We are currently revising MetaEuk. Our ISMB/ECCB 2019 poster is available [here](https://f1000research.com/posters/8-1489)


<p align="center"><img src="https://github.com/soedinglab/metaeuk/blob/master/imgs/MetaEuk.png" height="250"/></p>

## Installation
MetaEuk can be used by compiling from source (see below) or downloading a [statically compiled version](https://mmseqs.com/metaeuk/). It requires a 64-bit system (check with `uname -a | grep x86_64`) with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on MacOS).
     
     # static build sse4.1
     wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz; tar xvfz metaeuk-linux-sse41.tar.gz; export PATH=$(pwd)/metaeuk/bin/:$PATH
     # static build AVX2
     wget https://mmseqs.com/metaeuk/metaeuk-linux-avx2.tar.gz; tar xvfz metaeuk-linux-avx2.tar.gz; export PATH=$(pwd)/metaeuk/bin/:$PATH

## Input 
MetaEuk will search for eukaryotic protein-coding genes in **contigs** based on similarity to a reference database of **proteins** or **protein profiles**. The starting point are Fasta files of sequences (you can use contigs.fna and proteins.faa from the tests/two_contigs directory as a small toy example). Convert the contigs.fna file to a nucleotide database by running the createdb command (--dbtype 2)

Read [here](https://github.com/soedinglab/mmseqs2/wiki#how-to-create-a-target-profile-database-from-pfam) to learn more on how to create a protein profile database using MMseqs2. Once created, this database can be used as referenceDB in the command below.

## Running MetaEuk 
### Main Modules:

      predictexons      	Predict eukaryotic exons based on protein similarity
      reduceredundancy  	A greedy approach to group metaeuk predictions which share an exon
      unitesetstofasta  	Create a fasta output from optimal exon sets
      groupstoacc     	Create a TSV output from representative to member


### Important parameters: 

     --min-length        minimal number of codons in putative protein fragment
     -e                  maximal E-Value to retain a match between a putative protein fragment and a reference taraget 
     --metaeuk-eval      maximal combined E-Value to retain an optimal exon set
     --metaeuk-tcov      minimal length ratio of combined set to target 
     --slice-search      if refernceDB is a profile database, should be added
     

### Collecting optimal exons sets:

This module will extract all putative protein fragments from each contig (**C**) and strand (**S**), query them against the reference targets (**T**) and use dynamic programming to retain for each **T** the optimal compatible exon set from each **C** & **S** (thus creating **TCS** predictions).
    
    metaeuk predictexons contigsDB referenceDB predExResultDB tempFolder --metaeuk-eval 0.0001 -e 100 --min-length 40
    
Since this step involves a search, it is the most time-demanding of all analyses steps. Upon completion, it will output a database (contigs are keys), where each line contains information about a **TCS** prediction and its exon (multi-exon predictions will span several lines).


### Reducing redundancy:

If there are homologies in referenceDB (e.g., T1 is highly similar to T2), the same optimal exons set from a **C** & **S** combination will be predicted more than once. This module will group together **TCS** predictions that share and exon and will choose their representative. By default, it will greedily obtain a subset of the **TCS** representatives, such that there is no overlap of predictions on the same contig and strand (to allow same-strand overlaps, run with ```--overlap 1```).
    
    metaeuk reduceredundancy predExResultDB predRedResultDB predGroupsDB
    
Upon completion, it will output: predRedResultDB and predGroupsDB. predRedResultDB contains information about the **TCS** representatives (same format as predExResultDB). Each line of predGroupsDB maps from a representative to all **TCS** predictions that share an exon with it.



### Converting to Fasta:

The predExResultDB/predRedResultDB produced by the modules above, can be used to extract the sequences of the predicted protein-coding genes. The parameter ```--protein``` controls whether to transalte the coding genes (1) or report in nucleotides (0, default)
    
    metaeuk unitesetstofasta contigsDB referenceDB predRedResultDB predRedResultProteins.fas --protein 1
    


#### The MetaEuk header:

The header is composed of several sections, separated by pipes ('|'):

*>T_acc|C_acc|S|bitscore|E-Value|number_exons|low_coord|high_coord|exon1_coords|exon2_coords|...*

*coord* refers to the coordination on the contig. It is advisable to keep T_acc and C_acc short and without pipes. The exon_coords are of the structure:
*low[taken_low]:high[taken_high]:nucleotide_length[taken_nucleotide_length]*

Since MetaEuk allows for a very overlap on T of two putative exons (see P2 and P3 in the illustartion below), when joining the sequences of the exons, one of them is shortened.

<p align="center"><img src="https://github.com/soedinglab/metaeuk/blob/master/imgs/small_overlap_allowed.png" height="150"/></p>

Example header (two exons on the minus strand):
*>ERR1719262_736507|ERR868377_k119_2347399|-|341|6.2e-93|2|54324|54855|54855[54855]:54754[54754]:102[102]|54656[54668]:54324[54324]:333[321]*


### Creating a TSV map of prediction representatives to their group members:

A TSV file, of lines of the format:

*T_acc_rep|C_acc|S    T_acc_member|C_acc|S*

can help mapping from each representative prediction after the redundancy reduction stage to all its group members. Since redundancy reduction is performed per contig and strand combination, there will always be agreement in these fields. Note, a representative also maps to itself.

    metaeuk groupstoacc contigsDB referenceDB predGroupsDB predGroups.tsv
    


## Compile from source
Compiling MetaEuk from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile MetaEuk `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the MetaEuk binary will be located in the `build/bin` directory.

      git clone git@github.com:soedinglab/metaeuk.git .
      git submodule init
      git submodule update
      mkdir build
      cd build
      cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_MPI=1 -DCMAKE_INSTALL_PREFIX=. ..
      make -j
      make install
      export PATH="$(pwd)/bin/:$PATH"
        
:exclamation: If you want to compile metaeuk on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and MetaEuk will not be able to run multithreaded. Use the following cmake call:

      CXX="$(brew --prefix)/bin/g++-8" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

## Hardware requirements
MetaEuk will scale its memory consumption based on the available main memory of the machine. MetaEuk needs a CPU with at least the SSE4.1 instruction set to run. 


