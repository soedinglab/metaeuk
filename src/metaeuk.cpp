#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const int NO_CITATION = 0;
const char* binary_name = "metaeuk";
const char* tool_name = "metaeuk";
const char* tool_introduction = "Eukaryotic gene-discovery in metagenomic contigs";
const char* main_author = "Eli Levy Karin, eli.levy.karin@gmail.com";
const char* show_extended_help = "1";
const char* show_bash_info = "1";
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<struct Command> commands = {
        // Main tools (workflows for non-experts)
        {"predictexons",             predictexons,            &localPar.predictexonsworkflow,    COMMAND_MAIN,
                "Call optimal exon sets based on protein similarity",
                "An analog of 6-frame translation to produce putative protein fragments. Search against protein DB. Compatible exon set identified with respect to each target.",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:targetsDB> <o:calledExonsDB> <tmpDir>",
                NO_CITATION, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"calledExonsDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL},
                                {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"reduceredundancy",             reduceredundancy,            &localPar.reduceredundancy,    COMMAND_MAIN,
                "Cluster metaeuk calls that share an exon and select representative prediction",
                "A greedy examination of calls according to their contig order, subordered by the number of exons. Calls in a cluster share an exon with the representative.",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:calledExonsDB> <o:predictionsExonsDB> <o:predToCall>",
                NO_CITATION, {{"calledExonsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"predictionsExonsDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"predToCall", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb}}},
        {"unitesetstofasta",             unitesetstofasta,            &localPar.unitesetstofasta,    COMMAND_MAIN,
                "Create a fasta output from optimal exon sets",
                "Each optimal set is joined to a single sequence of codons or amino-acids",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:targetsDB> <i:exonsDB> <o:unitedExonsFasta>",
                NO_CITATION, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"exonsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"unitedExonsFasta", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}},
        {"groupstoacc",             groupstoacc,            &localPar.onlythreads,    COMMAND_MAIN,
                "Create a TSV output from representative prediction to member",
                "Replace the internal contig, target and strand identifiers with accessions from the headers",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:targetsDB> <i:predToCall> <o:predToCallInfoTSV>",
                NO_CITATION, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"predToCall", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"predToCallInfoTSV", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}},
        // internal modules (COMMAND_EXPERT)
        {"resultspercontig",             resultspercontig,            &localPar.collectoptimalset,    COMMAND_EXPERT,
                "Swap fragments against targets search and join with contigs",
                "Each contig lists all target hits, mediated by the extracted putative fragments from that contig",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDb> <i:fragmentsDb> <i:fragmentToTargetSearchRes> <o:contigToSearchRes>",
                NO_CITATION, {{"contigsDb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"fragmentsDb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"fragmentToTargetSearchRes", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                {"contigToSearchRes", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}},
        {"collectoptimalset",             collectoptimalset,            &localPar.collectoptimalset,    COMMAND_EXPERT,
                "Collect the optimal set of exons for a target protein/profile",
                "A dynamic programming procedure on all candidates of each contig and strand combination",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigToSearchRes> <i:targetsDB> <o:calledExonsDB>",
                NO_CITATION,{{"contigToSearchRes", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"calledExonsDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}}
};
