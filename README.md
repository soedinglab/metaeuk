# MetaEuk - sensitive, high-throughput gene discovery and annotation for large-scale eukaryotic metagenomics


[![BioConda Install](https://img.shields.io/conda/dn/bioconda/metaeuk.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/metaeuk)
[![Biocontainer Pulls](https://img.shields.io/endpoint?url=https%3A%2F%2Fmmseqs.com%2Fbiocontainer.php%3Fcontainer%3Dmetaeuk)](https://biocontainers.pro/#/tools/metaeuk)
[![Docker Pulls](https://img.shields.io/docker/pulls/soedinglab/metaeuk.svg)](https://hub.docker.com/r/soedinglab/metaeuk)
[![Build Status](https://dev.azure.com/elilevy/MetaEuk/_apis/build/status/soedinglab.metaeuk?branchName=master)](https://dev.azure.com/elilevy/MetaEuk/_build/latest?definitionId=2&branchName=master)
[![Chat Support](https://chat.mmseqs.com/api/v1/shield.svg?type=online&name=chat&icon=false)](https://chat.mmseqs.com/channel/metaeuk)

MetaEuk is a modular toolkit designed for large-scale gene discovery and annotation in eukaryotic metagenomic contigs. MetaEuk combines the fast and sensitive homology search capabilities of [MMseqs2](https://github.com/soedinglab/MMseqs2) with a dynamic programming procedure to recover optimal exons sets. It reduces redundancies in multiple discoveries of the same gene and resolves conflicting gene predictions on the same strand. MetaEuk is GPLv3-licensed open source software that is implemented in C++ and available for Linux and macOS. The software is designed to run efficiently on multiple cores.

## Publication

[Levy Karin E, Mirdita M and Soeding J. MetaEuk – sensitive, high-throughput gene discovery and annotation for large-scale eukaryotic metagenomics. Microbiome. 2020; 8:48](https://rdcu.be/b3ozK)

<p align="center"><img src="https://github.com/soedinglab/metaeuk/blob/master/imgs/MetaEuk.png" height="250"/></p>

## Installation
MetaEuk can be used by compiling from source (see below) or downloading a [statically compiled version](https://mmseqs.com/metaeuk/). It requires a 64-bit system (check with `uname -a | grep x86_64`) with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on MacOS).

```
# install via conda
conda install -c conda-forge -c bioconda metaeuk
# static Linux AVX2 build
wget https://mmseqs.com/metaeuk/metaeuk-linux-avx2.tar.gz; tar xzvf metaeuk-linux-avx2.tar.gz; export PATH=$(pwd)/metaeuk/bin/:$PATH
# static Linux SSE4.1 build
wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz; tar xzvf metaeuk-linux-sse41.tar.gz; export PATH=$(pwd)/metaeuk/bin/:$PATH
# static macOS build (universal binary with SSE4.1/AVX2/M1 NEON)
wget https://mmseqs.com/metaeuk/metaeuk-osx-universal.tar.gz; tar xzvf metaeuk-osx-universal.tar.gz; export PATH=$(pwd)/metaeuk/bin/:$PATH
```

Precompiled binaries for other architectures (ARM64, PPC64LE) and very old AMD/Intel CPUs (SSE2 only) are available at https://mmseqs.com/metaeuk.

## Input 
MetaEuk will search for eukaryotic protein-coding genes in **contigs** based on similarity to a reference set of **proteins** or **protein profiles**. The starting point are Fasta files of sequences (you can use contigs.fna and proteins.faa from the tests/two_contigs directory as a small toy example).

You could **either** use the ```easy-predict``` workflow directly on the Fasta files **or** convert them to databases by running the createdb command and later on specific MetaEuk modules.
Read [here](https://github.com/soedinglab/mmseqs2/wiki#how-to-create-a-target-profile-database-from-pfam) to learn more on how to create a protein profile database using MMseqs2. Once created, this database can be used as referenceDB in the commands below.

## Terminology
A **gene call** is an optimal set of exons predicted based on similarity to a specific target (**T**) in a specific contig (**C**) and strand (**S**). In the following it is referred to as a **TCS** or as a **call**. After redundancy reduction (see details below), the **representative TCS** is referred to as **prediction**.

## Running MetaEuk 
### Main Modules:

      easy-predict      	Predict proteins from contigs (fasta/db) based on similarities to targets (fasta/db) and return a fasta 
      predictexons      	Call optimal exon sets based on protein similarity
      reduceredundancy  	Cluster metaeuk calls which share an exon and select representative
      unitesetstofasta  	Create fasta output from optimal exon sets (and a TSV map between headers and internal identifiers)
      groupstoacc     	Create a TSV output from representative to calls
      taxtocontig     	Assign taxonomic labels to MetaEuk predictions and contigs by majority voting


### Important parameters: 

     --min-length        minimal number of codons in putative protein fragment
     -e                  maximal E-Value to retain a match between a putative protein fragment and a reference target 
     --metaeuk-eval      maximal combined E-Value to retain an optimal exon set
     --metaeuk-tcov      minimal length ratio of combined set to target 
     --slice-search      if referenceDB is a profile database, should be added
     

### easy-predict workflow:

This workflow combines the following MetaEuk modules into a single step: predictexons, reduceredundancy and unitesetstofasta (each of which is detailed below). Its inputs are contigs (either as a Fasta file or a previously created database) and targets (either as a FASTA file of protein sequences or a previously created database of proteins or protein profiles). It will run the modules and output the predictions in FASTA format.
    
    metaeuk easy-predict contigsFasta/contigsDB proteinsFasta/referenceDB predsResults tempFolder
    
It will result in **predsResults.fas** (protein sequences), **predsResults.codon.fas** and **predsResults.headersMap.tsv**


### Calling optimal exons sets:

This module will extract all putative protein fragments from each contig and strand, query them against the reference targets and use dynamic programming to retain for each **T** the optimal compatible exon set from each **C** & **S** (thus creating **TCS** calls).
    
    metaeuk predictexons contigsDB referenceDB callsResultDB tempFolder --metaeuk-eval 0.0001 -e 100 --min-length 40
    
Since this step involves a search, it is the most time-demanding of all analyses steps. Upon completion, it will output a database (contigs are keys), where each line contains information about a **TCS** and its exon (multi-exon **TCS**s will span several lines).


### Reducing redundancy:

If there are homologies in referenceDB (e.g., T1 is highly similar to T2), the same optimal exons set from a **C** & **S** combination will be called more than once. This module will group together **TCS**s that share an exon and will choose their representative **prediction**. By default, it will greedily obtain a subset of the **predictions**, such that there is no overlap of **predictions** on the same contig and strand (to allow same-strand overlaps, run with ```--overlap 1```).
    
    metaeuk reduceredundancy callsResultDB predsResultDB predGroupsDB
    
Upon completion, it will output: predsResultDB and predGroupsDB. predsResultDB contains information about the **predictions** (same format as callsResultDB). Each line of predGroupsDB maps from a **prediction** to all **TCS**s that share an exon with it.



### Converting to Fasta:

The callsResultDB/predsResultDB produced by the modules above, can be used to extract the sequences of the predicted protein-coding genes.
    
    metaeuk unitesetstofasta contigsDB referenceDB predsResultDB predsResults
    
It will result in **predsResults.fas** (protein sequences), **predsResults.codon.fas** and **predsResults.headersMap.tsv**


#### The MetaEuk header:

The header is composed of several sections, separated by pipes ('|'):

*>T_acc|C_acc|S|bitscore|E-Value|number_exons|low_coord|high_coord|exon1_coords|exon2_coords|...*

*coord* refers to the coordination on the contig (first base has coordinate 0). It is advisable to keep T_acc and C_acc short and without pipes. The exon_coords are of the structure:
*low[taken_low]:high[taken_high]:nucleotide_length[taken_nucleotide_length]*

Since MetaEuk allows for a very short overlap on T of two putative exons (see P2 and P3 in the illustration below), when joining the sequences of the exons, one of them is shortened. The coordinates of the codons taken from this exon will be in the square brackets (*[taken_low]*, *[taken_high]* and *[taken_nucleotide_length]*). These refer to the orange section of P3 below, while the coordinates outside the brackets refer to the yellow+orange section of P3.

<p align="center"><img src="https://github.com/soedinglab/metaeuk/blob/master/imgs/small_overlap_allowed.png" height="150"/></p>

Example header (two exons on the minus strand):

*>protein_acc|contig_acc|-|1146|0|2|3|1875|1875[1875]:970[970]:906[906]|893[869]:3[3]:891[867]*


### Creating a TSV map of predictions to their TCS group members:

A TSV file, of lines of the format:

*T_acc_rep|C_acc|S    T_acc_member|C_acc|S*

can help mapping from each representative prediction after the redundancy reduction stage to all its TCS group members. Since redundancy reduction is performed per contig and strand combination, there will always be agreement in these fields. Note, a representative also maps to itself.

    metaeuk groupstoacc contigsDB referenceDB predGroupsDB predGroups.tsv
    

### Taxonomic assignment with taxtocontig:

After obtaining MetaEuk predictions, the *taxtocontig* workflow allows assigning taxonomic labels to the predicted MetaEuk proteins and confer these predictions to their contigs. This workflow internally runs [*taxonomy*](https://github.com/soedinglab/MMseqs2/wiki#the-concept-of-lca) on the MetaEuk predictions, using any `--lca-mode`. It then performs majority voting among the taxonomically labeled predictions on a given contig to select a label for the contig. The parameter ```--majority``` indicates the minimal fraction of labeled predictions that agree in their taxonomic assignment (1.0 - consensus, 0.5 - at least 50%, etc.). The contig's label will be the last common ancestor (LCA) of the fraction of labeled predictions in agreement.

#### Example:
predictions' taxonomic labels: *Ostreococcus tauri*, *Ostreococcus mediterraneus*, *unclassified*, *Bathycoccus prasinos*
- contig label (`--majority 0.5`): *Ostreococcus* (genus), the LCA of 2 out of 3 labels
- contig label (`--majority 1`): *Bathycoccaceae* (family), the LCA of 3 out of 3 labels

#### Input:
- The output of a MetaEuk run: **contigsDB** (if you run MetaEuk with *easy-predict* you will find it at `<tmpDir>/latest/contigs`), **predsResults.fas** and **predsResults.headersMap.tsv**, which are produced by the *unitesetstofasta* module (called by *easy-predict*).
- A protein sequence database annotated with taxonomic information (**seqTaxDb**). See details [here](https://github.com/soedinglab/MMseqs2/wiki#creating-a-seqtaxdb). You could download such a resource with >88M entries [here](http://wwwuser.gwdg.de/~compbiol/metaeuk/2020_TAX_DB).

#### Command:
    metaeuk taxtocontig <i:contigsDB> <i:predsResults.fas> <i:predsResults.headersMap.tsv> <i:taxAnnotTargetDb> <o:taxResult> <tmpDir> --majority 0.5 --tax-lineage 1 --lca-mode 2
    
#### Output:
The run ends with two files: **taxResult_per_pred.tsv** and **taxResult_per_contig.tsv**, each of which is in [taxonomy result TSV format](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-output-and-tsv)

## Compile from source
Compiling MetaEuk from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile MetaEuk `git`, `g++` (4.9 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the MetaEuk binary will be located in the `build/bin` directory.

      git clone https://github.com/soedinglab/metaeuk.git .
      mkdir build
      cd build
      cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
      make -j
      make install
      export PATH="$(pwd)/bin/:$PATH"
        
:exclamation: If you want to compile metaeuk on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and MetaEuk will not be able to run multithreaded. Use the following cmake call:

      CXX="$(brew --prefix)/bin/g++-10" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

## Hardware requirements
MetaEuk will scale its memory consumption based on the available main memory of the machine. MetaEuk needs a CPU with at least the SSE4.1 instruction set to run. 


