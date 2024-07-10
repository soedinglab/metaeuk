# MetaEuk - sensitive, high-throughput gene discovery and annotation for large-scale eukaryotic metagenomics


[![BioConda Install](https://img.shields.io/conda/dn/bioconda/metaeuk.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/metaeuk)
[![Biocontainer Pulls](https://img.shields.io/endpoint?url=https%3A%2F%2Fmmseqs.com%2Fbiocontainer.php%3Fcontainer%3Dmetaeuk)](https://biocontainers.pro/#/tools/metaeuk)
[![Docker Pulls](https://img.shields.io/docker/pulls/soedinglab/metaeuk.svg)](https://hub.docker.com/r/soedinglab/metaeuk)
[![Build Status](https://dev.azure.com/elilevy/MetaEuk/_apis/build/status/soedinglab.metaeuk?branchName=master)](https://dev.azure.com/elilevy/MetaEuk/_build/latest?definitionId=2&branchName=master)

MetaEuk is a modular toolkit designed for large-scale gene discovery and annotation in eukaryotic metagenomic contigs. MetaEuk combines the fast and sensitive homology search capabilities of [MMseqs2](https://github.com/soedinglab/MMseqs2) with a dynamic programming procedure to recover optimal exons sets. It reduces redundancies in multiple discoveries of the same gene and resolves conflicting gene predictions on the same strand. MetaEuk is GPLv3-licensed open source software that is implemented in C++ and available for Linux and macOS. The software is designed to run efficiently on multiple cores.

<!--- TOC START -->
Table of Contents
-----------------
- [Publication](#publication)
- [Installation](#installation)
- [Input](#input)
- [Terminology](#terminology)
- [Running MetaEuk](#running-metaeuk)
    - [Main Modules](#main-modules)
    - [Using MMseqs2 commands within MetaEuk](#using-mmseqs2-commands-within-metaeuk)
    - [Important parameters](#important-parameters)
    - [easy-predict workflow](#easy-predict-workflow)
    - [Calling optimal exons sets](#calling-optimal-exons-sets)
    - [Reducing redundancy](#reducing-redundancy)
    - [Converting to Fasta and GFF](#converting-to-fasta-and-gff)
    - [Taxonomic assignment with taxtocontig](#taxonomic-assignment-with-taxtocontig)
- [Available reference databases](#available-reference-databases)
- [Compile from source](#compile-from-source)
- [Hardware requirements](#hardware-requirements)
<!--- TOC END -->



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
MetaEuk will search for eukaryotic protein-coding genes in **contigs** based on similarity to reference **proteins** or **protein profiles**. You could **either** use the ```easy-predict``` workflow directly on Fasta files **or** convert them to MMseqs2-formatted databases by running the `createdb` command and later on specific MetaEuk modules. Read [here](#available-reference-databases) about available reference database. You can use contigs.fna and proteins.faa from the tests/two_contigs directory as a small toy example.

## Terminology
A **gene call** is an optimal set of exons predicted based on similarity to a specific target (**T**) in a specific contig (**C**) and strand (**S**). In the following it is referred to as a **TCS** or as a **call**. After redundancy reduction (see details below), the **representative TCS** is referred to as **prediction**.

## Running MetaEuk 
### Main Modules:

      easy-predict      	Predict proteins from contigs (fasta/db) based on similarities to targets (fasta/db) and return a fasta & GFF
      predictexons      	Call optimal exon sets based on protein similarity
      reduceredundancy  	Cluster metaeuk calls which share an exon and select representative
      unitesetstofasta  	Create fasta output from optimal exon sets (and (1) a TSV map between headers and internal identifiers, (2) GFF summary)
      groupstoacc     	Create a TSV output from representative to calls
      taxtocontig     	Assign taxonomic labels to MetaEuk predictions and contigs by majority voting
      
 
### Using MMseqs2 commands within MetaEuk:
MMseqs2 commands are available through MetaEuk and no additional MMseqs2 installation is required.
For example, the MMseqs2 command `mmseqs createdb` can be replaced with `metaeuk createdb`, `mmseqs databases` with `metaeuk databases`, etc. Please see also the [MMseqs2 Wiki](https://github.com/soedinglab/MMseqs2/wiki) for more info about MMseqs2 commands.


### Important parameters: 

     --min-length        minimal number of codons in putative protein fragment
     -e                  maximal E-Value to retain a match between a putative protein fragment and a reference target 
     --metaeuk-eval      maximal combined E-Value to retain an optimal exon set
     --metaeuk-tcov      minimal length ratio of combined set to target 
     --exhaustive-search if referenceDB is a profile database, should be added (before version 4 called slice-search)
     --max-exon-sets     maximal number of exon sets on each contig and strand for a given target (from version 6)

### easy-predict workflow:

This workflow combines the following MetaEuk modules into a single step: predictexons, reduceredundancy and unitesetstofasta (each of which is detailed below). Its inputs are contigs (either as a Fasta file or a previously created database) and targets (either as a Fasta file of protein sequences or a previously created database of proteins or protein profiles). It will run the modules and output the predictions in Fasta format (as well as a GFF format).
    
    metaeuk easy-predict contigsFasta/contigsDB proteinsFasta/referenceDB predsResults tempFolder
    
It will result in **predsResults.fas** (protein sequences), **predsResults.codon.fas**, **predsResults.headersMap.tsv** and **predsResults.gff**.


### Calling optimal exons sets:

This module will extract all putative protein fragments from each contig and strand, query them against the reference targets and use dynamic programming to retain for each **T** the optimal compatible exon set from each **C** & **S** (thus creating **TCS** calls).
    
    metaeuk predictexons contigsDB referenceDB callsResultDB tempFolder --metaeuk-eval 0.0001 -e 100 --min-length 40
    
Since this step involves a search, it is the most time-demanding of all analyses steps. Upon completion, it will output a database (contigs are keys), where each line contains information about a **TCS** and its exon (multi-exon **TCS**s will span several lines).


#### OPTIONAL - calling of sub-optimal exon sets:

By default, MetaEuk calls a single and optimal compatible exon set from each **C** & **S** for each **T**. If you are interested in calling several matches to a certain **T** from each **C** & **S** (for example, to look for **gene duplications**), you can change the default value of ```max-exon-sets``` to the number of sets to look for (from version 6). A few important notes:

* If ```max-exon-sets``` > 1, then it is no longer guaranteed that ***TCS*** is a unique identifier. Therefore, when parsing the output of such runs, it is recommended to use ***TCS*** together with ***low_contig*** as the identifier (see details about the [MetaEuk header](https://github.com/soedinglab/metaeuk#the-metaeuk-header)).
* If I run with ```--max-exon-sets``` > 1, am I guaranteed to get ALL the predictions I get when running ```--max-exon-sets 1```? **No!** You most likely see all of them but this is not guaranteed because some complex cases can arise due to the redundancy reduction stage. You can see an example for such a case under tests/sub_opt/readme.txt.
* Running with ```max-exon-sets``` > 1 is mainly useful in case your contigs are long enough to contain several genes (less common in metagenomic data)


### Reducing redundancy:

If there are homologies in referenceDB (e.g., T1 is highly similar to T2), the same optimal exon set from a **C** & **S** combination will be called more than once. This module will group together **TCS**s that share an exon and will choose their representative **prediction**. By default, it will greedily obtain a subset of the **predictions**, such that there is no overlap of **predictions** on the same contig and strand (to allow same-strand overlaps, run with ```--overlap 1```).
    
    metaeuk reduceredundancy callsResultDB predsResultDB predGroupsDB
    
Upon completion, it will output: predsResultDB and predGroupsDB. predsResultDB contains information about the **predictions** (same format as callsResultDB). Each line of predGroupsDB maps from a **prediction** to all **TCS**s that share an exon with it.



### Converting to Fasta and GFF:

The callsResultDB/predsResultDB produced by the modules above, can be used to extract the sequences of the predicted protein-coding genes.
    
    metaeuk unitesetstofasta contigsDB referenceDB predsResultDB predsResults
    
It will result in **predsResults.fas** (protein sequences), **predsResults.codon.fas**, **predsResults.headersMap.tsv** and **predsResults.gff**


#### The MetaEuk header:

The basic header is composed of several sections, separated by pipes ('|'):

*>T_acc|C_acc|S|bitscore|E-Value|number_exons|low_coord|high_coord|exon1_coords|exon2_coords|...*

*coord* refers to the coordinates on the contig (first base has coordinate 0). It is advisable to keep T_acc and C_acc short and without pipes. The exon_coords are of the structure:
*low[taken_low]:high[taken_high]:nucleotide_length[taken_nucleotide_length]*

Since MetaEuk allows for a very short overlap on T of two putative exons (see P2 and P3 in the illustration below), when joining the sequences of the exons, one of them is shortened. The coordinates of the codons taken from this exon will be in square brackets (*[taken_low]*, *[taken_high]* and *[taken_nucleotide_length]*). These refer to the orange section of P3 below, while the coordinates outside the brackets refer to the yellow+orange section of P3.

<p align="center"><img src="https://github.com/soedinglab/metaeuk/blob/master/imgs/small_overlap_allowed.png" height="150"/></p>

Example header (two exons on the minus strand):

*>protein_acc|contig_acc|-|1146|0|2|3|1875|1875[1875]:970[970]:906[906]|893[869]:3[3]:891[867]*

##### OPTIONAL - adding information about stop codon positions:
By setting the flag `--write-frag-coords 1`, information about the position of stop codons will be added to the output. In this case the exon_coords will be given in the following structure:

*[fragment_low]low[taken_low]:[fragment_high]high[taken_high]:nucleotide_length[taken_nucleotide_length]*

In its initial stage, MetaEuk extracts putative coding fragments between stop codons. It later discovers exons within them by matching targets. The fragment coordinates in square brackets refer to the original fragment in which the exon was found. In addition to reporting these coordinates, MetaEuk will print the stop codon (`*` in the protein output) right at the end of the last exon, if it exists.

##### OPTIONAL - scanning for start codon before the first exon:
By default (`--len-scan-for-start 0`), MetaEuk only reports parts of the contig that match a target. In case of fragmented targets or very distantly-related targets, it can therefore produce predictions, which do not start with a methionine. By setting `--len-scan-for-start` to a positive number, e.g., 50, MetaEuk will scan up-to 50 nucleotides (16 codons) before the first exon of each prediction (upstream for predictions on the plus strand, downstream - for minus). 

The scan will be in the same frame as the first exon and not beyond its stop codon border. Within this "legal" window, the scan will finish at the closest methionine to the first exon's matched start. The fragment from the found ATG until the first exon's matched start will be padded to the reported sequence. In the case of predictions on the plus strand, the *low_coord* value (7th field) will be updated to a lower value and the length of the padded fragment will be reported in square brackets. If the scan is turned on but no padding occurred (if the prediction already started with methionine or if no ATG was found), then the *low_coord* value will remain the same, followed by 0 in square brackets. For predictions on the minus strand, the change will be to *high_coord* (8th field). All other fields, including the exon fields, will remain unchanged. Examples:

*>protein_acc|contig_acc|+|784|1.213e-233|4|100[18]|1444|...*
 Here, six codons including ATG, were padded before the first exon of a prediction on the plus strand, which starts at position 118 (100+18).

*>protein_acc|contig_acc|-|499|7.54e-148|2|100|911[12]|...*
 Here, four codons including ATG, were padded before the first exon of a prediction on the minus strand, which starts at position 899 (911-12).

*>protein_acc|contig_acc|-|499|7.54e-148|2|100|899[0]|...*
 Here, no padding occurred, but the scan option was set to a positive number.


Of note, for simplicity, MetaEuk considers only ATG as a start for this scan.


##### The MetaEuk GFF:

In addition to writing a Fasta file, MetaEuk writes a GFF file. Please note that GFF is not perfectly suitable for MetaEuk because MetaEuk doesn't predict non-coding regions. This means that by default the MetaEuk gene starts and ends where the first and last codons could be matched (or slightly padded if `--len-scan-for-start` is set to be positive, see section). The gene and mRNA categories are the same in the MetaEuk GFF (if `--len-scan-for-start` is set to be positive, these fields will reflect the padding, as explained). The exon and CDS coordinates will be the same unless a small target overlap was allowed, due to which, the MetaEuk exon was shortened (see above). In this case, the CDS will report the shortening. In the sixth column you can find their individual bitsocres. Unlike MetaEuk's native report in the Fasta header, the contig index starts at 1 and the start coordinate is always smaller than the end coordinate, as required by GFF. The last column contains the **TCS** identifier, followed by the low_coord of the prediction to support searching for sub-optimal exon sets (see section). Here is an example where a MetaEuk header of two exons is reported in GFF format:

*>protein_acc|contig_acc|-|508|1.15e-150|2|100|911|911[911]:582[582]:330[330]|501[501]:100[100]:402[402]*


    contig_acc    MetaEuk    gene    101     912     508     -       .       Target_ID=protein_acc;TCS_ID=protein_acc|contig_acc|-|low_coord
    contig_acc    MetaEuk    mRNA    101     912     508     -       .       Target_ID=protein_acc;TCS_ID=protein_acc|contig_acc|-|low_coord_mRNA;Parent=protein_acc|contig_acc|-|low_coord
    contig_acc    MetaEuk    exon    583     912     234     -       .       Target_ID=protein_acc;TCS_ID=protein_acc|contig_acc|-|low_coord_exon_0;Parent=protein_acc|contig_acc|-|low_coord_mRNA
    contig_acc    MetaEuk    CDS    583     912     234     -       .       Target_ID=protein_acc;TCS_ID=protein_acc|contig_acc|-|low_coord_CDS_0;Parent=protein_acc|contig_acc|-|low_coord_exon_0
    contig_acc    MetaEuk    exon    101     502     273     -       .       Target_ID=protein_acc;TCS_ID=protein_acc|contig_acc|-|low_coord_exon_1;Parent=protein_acc|contig_acc|-|low_coord_mRNA
    contig_acc    MetaEuk    CDS    101     502     273     -       .       Target_ID=protein_acc;TCS_ID=protein_acc|contig_acc|-|low_coord_CDS_1;Parent=protein_acc|contig_acc|-|low_coord_exon_1



### Creating a TSV map of predictions to their TCS group members:

A TSV file, of lines of the format (low_coord information added in version 6):

*T_acc_rep|C_acc|S|low_coord_rep    T_acc_member|C_acc|S|low_coord_member*

can help mapping from each representative prediction after the redundancy reduction stage to all its TCS group members. Since redundancy reduction is performed per contig and strand combination, there will always be agreement in these fields. Note, a representative also maps to itself.

    metaeuk groupstoacc contigsDB referenceDB predGroupsDB predGroups.tsv
    

### Taxonomic assignment with taxtocontig:

After obtaining MetaEuk predictions, the *taxtocontig* workflow allows assigning taxonomic labels to the predicted MetaEuk proteins and confer these predictions to their contigs. This workflow internally runs [*taxonomy*](https://github.com/soedinglab/MMseqs2/wiki#the-concept-of-lca) on the MetaEuk predictions, using any `--lca-mode`. It then performs majority voting among the taxonomically labeled predictions on a given contig to select a label for the contig. The parameter ```--majority``` indicates the minimal fraction of labeled predictions that agree in their taxonomic assignment (1.0 - consensus, 0.5 - at least 50%, etc.). The contig's label will be the last common ancestor (LCA) of the fraction of labeled predictions in agreement. Please note that MMseqs2 commands are avaialble through MetaEuk.

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

## Available reference databases
Any Fasta file containing protein sequences or MMseqs2-formatted database of proteins or protein profiles can be provided as a reference database to MetaEuk. 

Don't have one yet? Not a problem! Here is what you can do:
* Using the `databases` command, you can easily download several of the publicly available databases, as detailed [here](https://github.com/soedinglab/MMseqs2/wiki#downloading-databases). 
Conveniently, many of these databases will be downloaded with taxonomic information, which will both allow you to filter them according to your need (for example, retain only eukaryotic sequences), as detailed [here](https://github.com/soedinglab/MMseqs2/wiki#filtering-a-seqtaxdb) and use them for [taxonomic assignment with MetaEuk](#taxonomic-assignment-with-taxtocontig) at a later stage, if desired.

   Read [here](https://github.com/soedinglab/mmseqs2/wiki#how-to-create-a-target-profile-database-from-pfam) to learn more on how to create a protein profile database.

* Additional resources include two databases [released alongside the MetaEuk publication](https://wwwuser.gwdg.de/~compbiol/metaeuk/). These are focused on Eukaryotes in the marine environment. The first contains [~6 million proteins predicted by MetaEuk](https://wwwuser.gwdg.de/~compbiol/metaeuk/2019_11/MetaEuk_preds_Tara_vs_euk_profiles_uniqs.fas.gz). The second consists of [~88 protein profiles](https://wwwuser.gwdg.de/~compbiol/metaeuk/2019_11/MERC_MMETSP_Uniclust50_profiles.tar.gz) created from, among others, eukaryotic proteins from the marine environment. Of note, due to changes in the profile database format, the second resource has been updated (Dec 2022). If you for some reason need the old format, you can find it under the abovementioned release folder. 


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

      CC="$(brew --prefix)/bin/gcc-13" CXX="$(brew --prefix)/bin/g++-13" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

## Hardware requirements
MetaEuk will scale its memory consumption based on the available main memory of the machine. MetaEuk needs a CPU with at least the SSE2 instruction set to run. 


