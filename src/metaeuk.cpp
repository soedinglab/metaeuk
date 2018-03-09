#include "Command.h"
#include "CommandDeclarations.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "metaeuk";
const char* tool_name = "metaeuk";
const char* tool_introduction = "Metagenomic tool for Eukayotic data.";
const char* main_author = "Eli Levy Karin, eli.levy.karin@gmail.com";

LocalParameters& par = LocalParameters::getLocalInstance();

std::vector<struct Command> commands = {
        // Main tools  (for non-experts)
        {"predictexons",             predictexons,            &par.predictexonsworkflow,    COMMAND_MAIN,
                "Predict eukaryotic exons based on protein similarity.",
                "An analog of 6-frame translation to produce putative protein fragments. Search against protein DB. Compatible exon set identified with respect to each target.",
                "Eli Levy Karin <eli.levy.karin@gmail.com> ",
                "<i:sequenceDB> <i:proteinTargetsDB> <o:outDB> <tmpDir>",
                CITATION_MMSEQS2},
        {"collectoptimalset",             collectoptimalset,            &par.predictexonsworkflow,    COMMAND_MAIN,
                "Collect the optimal set of exons for a target protein",
                "A dynamic programming procedure on all candidates of each contig and strand combination",
                "Eli Levy Karin <eli.levy.karin@gmail.com> ",
                "<i:targetToQueryWithContigInfoDB> <o:outDB>",
                CITATION_MMSEQS2},
        {"createdb",             createdb,             &par.createdb,             COMMAND_MAIN,
                "Convert protein sequence set in a FASTA file to MMseqs sequence DB format",
                "converts a protein sequence flat/gzipped FASTA or FASTQ file to the MMseqs sequence DB format. This format is needed as input to mmseqs search, cluster and many other tools.",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <o:sequenceDB>",
                CITATION_MMSEQS2},
        {"extractorfs",          extractorfs,          &par.extractorfs,          COMMAND_DB,
                "Extract open reading frames from all six frames from nucleotide sequence DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2},
        {"translatenucs",        translatenucs,        &par.translatenucs,        COMMAND_DB,
                "Translate nucleotide sequence DB into protein sequence DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2},
        {"search",               search,               &par.searchworkflow,       COMMAND_MAIN,
                "Search with query sequence or profile DB (iteratively) through target sequence DB",
                "Searches with the sequences or profiles query DB through the target sequence DB by running the prefilter tool and the align tool for Smith-Waterman alignment. For each query a results file with sequence matches is written as entry into a database of search results (alignmentDB).\nIn iterative profile search mode, the detected sequences satisfying user-specified criteria are aligned to the query MSA, and the resulting query profile is used for the next search iteration. Iterative profile searches are usually much more sensitive than (and at least as sensitive as) searches with single query sequences.",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_MMSEQS2},
        // required for search:
        {"prefilter",            prefilter,            &par.prefilter,            COMMAND_EXPERT,
                "Search with query sequence / profile DB through target DB (k-mer matching + ungapped alignment)",
                "Searches with the sequences or profiles in query DB through the target sequence DB in two consecutive stages: a very fast k-mer matching stage (double matches on same diagonal) and a subsequent ungapped alignment stage. For each query a results file with sequence matches is written as entry into the prefilter DB.",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> & Maria Hauser",
                "<i:queryDB> <i:targetDB> <o:prefilterDB>",
                CITATION_MMSEQS2},
        // required for search:
        {"align",                align,                &par.align,                COMMAND_EXPERT,
                "Compute Smith-Waterman alignments for previous results (e.g. prefilter DB, cluster DB)",
                "Calculates Smith-Waterman alignment scores between all sequences in the query database and the sequences of the target database which passed the prefiltering.",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> & Maria Hauser",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:alignmentDB>",
                CITATION_MMSEQS2},
        {"swapresults",          swapresults,          &par.swapresult,          COMMAND_DB,
                "Reformat prefilter/alignment/cluster DB as if target DB had been searched through query DB",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> & Clovis Galiez",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2},
        {"filterdb",             filterdb,             &par.filterDb,             COMMAND_DB,
                "Filter a DB by conditioning (regex, numerical, ...) on one of its whitespace-separated columns",
                NULL,
                "Clovis Galiez & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2},
        // Special-purpose utilities
        {"convert2fasta",        convert2fasta,        &par.convert2fasta,        COMMAND_FORMAT_CONVERSION,
                "Convert sequence DB to FASTA format",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:sequenceDB> <o:fastaFile>",
                CITATION_MMSEQS2},
        {"shellcompletion",      shellcompletion,      &par.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_MMSEQS2},
};
