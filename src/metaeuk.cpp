#include "Command.h"
#include "DownloadDatabase.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "Prefiltering.h"

const char* binary_name = "metaeuk";
const char* tool_name = "metaeuk";
const char* tool_introduction = "MetaEuk is homology-based strategy to efficiently query many contigs assembled from metagenomic samples against a comprehensive protein/profile target database to describe their protein repertoire. It does not require preliminary binning of the contigs and makes no assumption concerning the splicing signal when searching for multi-exon proteins.\n\nPlease cite:\nLevy Karin E, Mirdita M, Soding J: MetaEuk — sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics. Microbiome (2020) 8:48.";
const char* main_author = "Eli Levy Karin, eli.levy.karin@gmail.com";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
extern const char* MMSEQS_CURRENT_INDEX_VERSION;
const char* index_version_compatible = MMSEQS_CURRENT_INDEX_VERSION;
bool hide_base_commands = true;
bool hide_base_downloads = false;
void (*validatorUpdate)(void) = 0;
std::vector<KmerThreshold> externalThreshold = {};

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<struct Command> commands = {
        // Main tools (workflows for non-experts)
        {"predictexons",             predictexons,            &localPar.predictexonsworkflow,    COMMAND_MAIN,
                "Call optimal exon sets based on protein similarity",
                "An analog of 6-frame translation to produce putative protein fragments. Search against protein DB. Compatible exon set identified with respect to each target.",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:targetsDB> <o:calledExonsDB> <tmpDir>",
                CITATION_METAEUK, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                   {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                   {"calledExonsDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL},
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"easy-predict",             easypredict,            &localPar.easypredictworkflow,    COMMAND_MAIN,
                "Predict protein-coding genes from contigs (fasta/database) based on similarities to targets (fasta/database) and return a fasta of the predictions in a single step",
                "Combines the following MetaEuk modules into a single step: predictexons, reduceredundancy and unitesetstofasta",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigs> <i:targets> <o:predictionsFasta> <tmpDir>",
                CITATION_METAEUK, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                   {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                   {"predictionsFasta", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"taxtocontig",             taxtocontig,            &localPar.taxpercontigworkflow,    COMMAND_MAIN,
                "Assign taxonomic labels to predictions and aggregate them per contig",
                "The LCA of a majority of predictions will be assigned to their contig",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:predictionsFasta> <i:predictionsFasta.headersMap.tsv> <i:taxAnnotTargetDb> <o:taxResult> <tmpDir>",
                CITATION_METAEUK, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                   {"predictionsFasta", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                   {"predictionsFasta.headersMap.tsv", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                   {"taxAnnotTargetDb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb},
                                   {"taxResult", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxResult},
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"reduceredundancy",             reduceredundancy,            &localPar.reduceredundancy,    COMMAND_MAIN,
                "Cluster metaeuk calls that share an exon and select representative prediction",
                "A greedy examination of calls according to their contig order, subordered by the number of exons. Calls in a cluster share an exon with the representative.",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:calledExonsDB> <o:predictionsExonsDB> <o:predToCall>",
                CITATION_METAEUK, {{"calledExonsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                   {"predictionsExonsDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                   {"predToCall", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb}}},
        {"unitesetstofasta",             unitesetstofasta,            &localPar.unitesetstofasta,    COMMAND_MAIN,
                "Create a fasta output from optimal exon sets",
                "Each optimal set is joined to a single sequence of codons or amino-acids. In addition, a TSV map for each header to internal idenfiers.",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:targetsDB> <i:exonsDB> <o:unitedExonsFasta>",
                CITATION_METAEUK, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                   {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                   {"exonsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                   {"unitedExonsFasta", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}},
        {"groupstoacc",             groupstoacc,            &localPar.onlythreads,    COMMAND_MAIN,
                "Create a TSV output from representative prediction to member",
                "Replace the internal contig, target and strand identifiers with accessions from the headers",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:targetsDB> <i:predToCall> <o:predToCallInfoTSV>",
                CITATION_METAEUK, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                   {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                   {"predToCall", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                   {"predToCallInfoTSV", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}},
        // internal modules (COMMAND_EXPERT)
        {"resultspercontig",             resultspercontig,            &localPar.collectoptimalset,    COMMAND_EXPERT,
                "Swap fragments against targets search and join with contigs",
                "Each contig lists all target hits, mediated by the extracted putative fragments from that contig",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDb> <i:fragmentsDb> <i:fragmentToTargetSearchRes> <o:contigToSearchRes>",
                CITATION_METAEUK, {{"contigsDb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                   {"fragmentsDb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                   {"fragmentToTargetSearchRes", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                   {"contigToSearchRes", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}},
        {"collectoptimalset",             collectoptimalset,            &localPar.collectoptimalset,    COMMAND_EXPERT,
                "Collect the optimal set of exons for a target protein/profile",
                "A dynamic programming procedure on all candidates of each contig and strand combination",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigToSearchRes> <i:targetsDB> <o:calledExonsDB>",
                CITATION_METAEUK,{{"contigToSearchRes", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                  {"targetsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                  {"calledExonsDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}}
};

std::vector<DatabaseDownload> externalDownloads = {};
