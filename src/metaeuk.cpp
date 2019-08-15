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
                "Predict eukaryotic exons based on protein similarity",
                "An analog of 6-frame translation to produce putative protein fragments. Search against protein DB. Compatible exon set identified with respect to each target. Two dynamic programming outputs: predictexonsBaseName_dp_protein_contig_strand_map and predictexonsBaseName_dp_optimal_exon_sets",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:proteinsDB> <o:predictexonsBaseName> <tmpDir>",
                NO_CITATION, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                {"proteinsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"predictexonsBaseName", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL},
                                {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"reduceredundancy",             reduceredundancy,            &localPar.reduceredundancyworkflow,    COMMAND_MAIN,
                "A greedy approach to group metaeuk predictions which share an exon",
                "A protein coding gene can be predicted more than once due to target DB homologies. A cluster representative is selected. Predictions in a cluster share an exon with the representative. Outputs: reduceRedundBaseName_grouped_predictions, reduceRedundBaseName_dp_protein_contig_strand_map, reduceRedundBaseName_dp_optimal_exon_sets, reduceRedundBaseName_grouped_predictions_no_overlap, reduceRedundBaseName_no_overlap_dp_protein_contig_strand_map and reduceRedundBaseName_no_overlap_dp_optimal_exon_sets",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:dp_protein_contig_strand_map> <i:dp_optimal_exon_sets> <o:reduceRedundBaseName> <tmpDir>",
                NO_CITATION, {{"dp_protein_contig_strand_map", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"dp_optimal_exon_sets", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"reduceRedundBaseName", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL},
                                {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        // taxonomy workflow is under construction...
        {"assigntaxonomy",             assigntaxonomy,            &localPar.assigntaxonomyworkflow,    COMMAND_HIDDEN,
                "Assign taxonomy to each representative prediction using 2bLCA",
                "Assign taxonomy to each representative prediction using 2bLCA against a uniprot-based refernce database. Contigs are assigned based on their predictions",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:metaeukBaseName> <i:uniprotFasta> <tmpDir>",
                NO_CITATION, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"unitetoseqdbs",             unitetoseqdbs,            &localPar.unitetoseqdbsworkflow,    COMMAND_MAIN,
                "Unite the exons of predictions to sequence BD",
                "Unite and translate the exons of predictions (either before or after redundancy reduction - depending on the dp files) to sequence BD. Outputs: unitedexonsBasename_united_exons, unitedexonsBasename_united_exons_aa",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:proteinsDB> <i:dp_protein_contig_strand_map> <i:dp_optimal_exon_sets> <o:unitedexonsBasename> <tmpDir>",
                NO_CITATION, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                {"proteinsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"dp_protein_contig_strand_map", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"dp_optimal_exon_sets", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"unitedexonsBasename", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL},
                                {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        // internal modules (COMMAND_EXPERT)
        {"resultspercontig",             resultspercontig,            &localPar.collectoptimalset,    COMMAND_EXPERT,
                "resultspercontig",
                "resultspercontig",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDb> <i:fragmentsDb> <i:fragmentToTargetSearchRes> <o:contigToSearchRes>",
                NO_CITATION,{{"contigsDb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"fragmentsDb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"fragmentToTargetSearchRes", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                {"contigToSearchRes", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}},
        {"collectoptimalset",             collectoptimalset,            &localPar.collectoptimalset,    COMMAND_EXPERT,
                "Collect the optimal set of exons for a target protein",
                "A dynamic programming procedure on all candidates of each contig and strand combination",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:proteinToPotentialExonsWithContigInfoDB> <i:proteinsDB> <o:predictexonsBaseName>",
                NO_CITATION,{{"proteinToPotentialExonsWithContigInfoDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                {"proteinsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"predictexonsBaseName", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, NULL}}},
        {"grouppredictions",             grouppredictions,            &localPar.onlythreads,    COMMAND_EXPERT,
                "Assignment of predictions to cluster",
                "A greedy examination of predictions accoridng to their contig order, subordered by the number of exons. Predictions in a cluster share an exon with the representative.",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigStrandSortedMap> <o:groupedPredictionsDB> <o:groupedPredictionsDBNoOverlap>",
                NO_CITATION, {{"contigStrandSortedMap", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"groupedPredictionsDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb},
                                {"groupedPredictionsDBNoOverlap", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb}}},
        {"unitesetstosequencedb",             unitesetstosequencedb,            &localPar.onlythreads,    COMMAND_EXPERT,
                "Create a sequence DB from optimal exon sets",
                "Each optimal set is joined to a single sequence",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:contigsDB> <i:proteinsDB> <i:dp_protein_contig_strand_map> <i:dp_optimal_exon_sets> <o:unitedexons>",
                NO_CITATION, {{"contigsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb},
                                {"proteinsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"dp_protein_contig_strand_map", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"dp_optimal_exon_sets", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                {"unitedexons", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb}}},
        // overwrite MMseqs2 swapdb because of its input type validation
        {"swapdb",               swapdb,               &localPar.swapdb,               COMMAND_HIDDEN,
                "Create a DB where the key is from the first column of the input result DB",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>, Clovis Galiez & Eli Levy Karin",
                "<i:resultDB> <o:resultDB>",
                NO_CITATION, {{"resultDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                   {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb}}}
};
