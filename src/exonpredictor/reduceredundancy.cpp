#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <sys/time.h>

#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include <limits>
#include <cstdint>
#include <queue>
#include <string.h>
#include "itoa.h"
#include "PredictionParser.h"

#ifdef OPENMP
#include <omp.h>
#endif

const size_t EXPECTED_NUM_PREDICTIONS = 100000;

void clusterPredictions (std::vector<Prediction> &contigPredictions, std::vector<Prediction> &repContigPredictions) {
    // the index i iterates over cluster tmp_representatives.
    // after member collection is done, tmp_representative is replaced by the member with the highest bitscore
    for (size_t i = 0; i < contigPredictions.size(); ++i) {
        // if i is already assigned - skip it, it is not a cluster representative
        if (contigPredictions[i].isClustered) {
            continue;
        }

        // collect cluster members - i is the tmp_representative:
        size_t finalClusterId = contigPredictions[i].targetKey;
        contigPredictions[i].clusterId = contigPredictions[i].targetKey;
        unsigned int maxScore = contigPredictions[i].totalBitscore;
        contigPredictions[i].isClustered = true;
        size_t max_j = i;
        
        for (size_t j = (i + 1); j < contigPredictions.size(); ++j) {
            if (contigPredictions[j].lowContigCoord >= contigPredictions[i].highContigCoord) {
                // overlap is over - no need to compare to other j - finish and move to the next i
                break;
            }
            
            bool doIandJshareAnExon = false;
            for (size_t exon_id_i = 0; exon_id_i < contigPredictions[i].optimalExonSet.size(); ++exon_id_i) {
                for (size_t exon_id_j = 0; exon_id_j < contigPredictions[j].optimalExonSet.size(); ++exon_id_j) {
                    if (contigPredictions[i].optimalExonSet[exon_id_i].exonKey == contigPredictions[j].optimalExonSet[exon_id_j].exonKey) {
                        doIandJshareAnExon = true;
                        goto endExonComparison;
                    }
                }
            }

            // assign j to the cluster of i if it has a mutual exon and if it wasn't already assigned:
            endExonComparison:
            if (doIandJshareAnExon && (! contigPredictions[j].isClustered)) {
                contigPredictions[j].isClustered = true;
                // tmp representative is i:
                contigPredictions[j].clusterId = contigPredictions[i].targetKey;

                // if bitscore of j is better - update the cluster identifier to be the target key of j
                if (contigPredictions[j].totalBitscore > maxScore) {
                    maxScore = contigPredictions[j].totalBitscore;
                    finalClusterId = contigPredictions[j].targetKey;
                }

                // will be used to assign finalClusterId later
                if (max_j < j) {
                    max_j = j;
                }
            }
        }

        // collecting all j members for tmp_representative i is finished.
        // assign finalClusterId (target key of best scoring representative)
        for (size_t j = i; j <= max_j; ++j) {
            if (contigPredictions[j].clusterId == contigPredictions[i].targetKey) {
                contigPredictions[j].clusterId = finalClusterId;
            }
            // if by now the clusterId is equal to finalClusterId --> this is a representative
            if (contigPredictions[j].clusterId == contigPredictions[j].targetKey) {
                repContigPredictions.emplace_back(contigPredictions[j]);
            }
        }
    }
}

void excludeSameStrandOverlaps (std::vector<Prediction> &repContigPredictions) {
    // sort vector by E-value:
    std::stable_sort(repContigPredictions.begin(), repContigPredictions.end(), Prediction::comparePredictionsByEvalue);
    // go over sorted representatives - unassigned representatives are non-overlap representatives
    for (size_t i = 0; i < repContigPredictions.size(); ++i) {
        // if i is already assigned - skip it, it is not a non-overlap cluster representative
        if (repContigPredictions[i].isNoOverlapClustered) {
            continue;
        }
        // initialize the new cluster:
        repContigPredictions[i].isNoOverlapClustered = true;
        repContigPredictions[i].noOverlapClusterId = repContigPredictions[i].targetKey;

        // collect cluster members:
        for (size_t j = (i + 1); j < repContigPredictions.size(); ++j) {
            if (
                ((repContigPredictions[j].highContigCoord < repContigPredictions[i].highContigCoord) && 
                (repContigPredictions[j].highContigCoord > repContigPredictions[i].lowContigCoord)) ||

                ((repContigPredictions[j].lowContigCoord < repContigPredictions[i].highContigCoord) && 
                (repContigPredictions[j].lowContigCoord > repContigPredictions[i].lowContigCoord)) ||

                ((repContigPredictions[j].highContigCoord < repContigPredictions[i].highContigCoord) && 
                (repContigPredictions[j].lowContigCoord > repContigPredictions[i].lowContigCoord)) ||

                ((repContigPredictions[j].highContigCoord > repContigPredictions[i].highContigCoord) && 
                (repContigPredictions[j].lowContigCoord < repContigPredictions[i].lowContigCoord))
                ) {
                    // i begins in the middle of j or j begins in the middle of i ==> overlap
                    repContigPredictions[j].isNoOverlapClustered = true;
                    repContigPredictions[j].noOverlapClusterId = repContigPredictions[i].targetKey;
            }
        }
    }
}

void writeRepPredsInDPFormat (std::vector<Prediction> &repContigPredictions, char * predBuff, bool allowOverlaps,
                                DBWriter &repWriter, unsigned int thread_idx) {
    for (size_t i = 0; i < repContigPredictions.size(); ++i) {
        // if same strand overlaps are not allowed, skip predictions that were worse than another representatives
        if ((allowOverlaps == false) && (repContigPredictions[i].noOverlapClusterId != repContigPredictions[i].targetKey)) {
            continue;
        }
        size_t predLen = Prediction::predictionToBuffer(predBuff, repContigPredictions[i]);
        repWriter.writeAdd(predBuff, predLen, thread_idx);
    }
}

void writePredsClusters (std::vector<Prediction> &predictions, char * predBuff, DBWriter &repWriter, unsigned int thread_idx) {
    for (size_t i = 0; i < predictions.size(); ++i) {
        size_t predLen = Prediction::predictionClusterToBuffer(predBuff, predictions[i]);
        repWriter.writeAdd(predBuff, predLen, thread_idx);
    }
}

int reduceredundancy(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, true, 0, 0);

    // db1 = input, predictions per contig
    DBReader<unsigned int> predsPerContig(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    predsPerContig.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // db2 = output, DP format of representative predictions (par.overlapAllowed will exclude overlaps by default)
    DBWriter writerGroupedPredictions(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writerGroupedPredictions.open();

    // db3 = output, grouping of predictions: T,S of representatives to T,S of prediction
    DBWriter writerRepToMembers(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writerRepToMembers.open();

    Debug::Progress progress(predsPerContig.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        // per thread variables
        const char *entry[255];

        std::vector<Prediction> plusContigPredictions;
        plusContigPredictions.reserve(EXPECTED_NUM_PREDICTIONS);
        std::vector<Prediction> plusContigRepPreds;

        std::vector<Prediction> minusContigPredictions;
        minusContigPredictions.reserve(EXPECTED_NUM_PREDICTIONS);
        std::vector<Prediction> minusContigRepPreds;

        char predictionBuffer[5000];

#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < predsPerContig.getSize(); id++) {
            progress.updateProgress();

            unsigned int contigKey = predsPerContig.getDbKey(id);

            char *results = predsPerContig.getData(id, thread_idx);
            
            // for verifying legal input
            unsigned int prevTargetKey = 0;
            // for predictions collection
            unsigned int prevTargetKeyPlus = 0;
            unsigned int prevTargetKeyMinus = 0;
            bool isFirstIterationPlus = true;
            bool isFirstIterationMinus = true;

            // keep track of offset when a contig starts
            writerRepToMembers.writeStart(thread_idx);
            writerGroupedPredictions.writeStart(thread_idx);

            // process a specific contig
            while (*results != '\0') {
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // each line informs of a prediction and a single exon
                // the first 7 columns describe the entire prediction
                // the last 10 columns describe a single exon
                if (columns != 17) {
                    Debug(Debug::ERROR) << "There should be 17 columns in the input file. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int targetKey = Prediction::getTargetKey(entry);
                int strand = Prediction::getStrand(entry);

                // verify legal input
                if (prevTargetKey > targetKey) {
                    Debug(Debug::ERROR) << "Predictions are assumed to be sorted by their target keys. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }
                prevTargetKey = targetKey;

                if (strand == PLUS) {
                    if ((isFirstIterationPlus == true) || (prevTargetKeyPlus != targetKey)) {
                        plusContigPredictions.emplace_back(Prediction());
                        plusContigPredictions.back().setByDPRes(entry);
                        isFirstIterationPlus = false;
                    }
                    // add the exon key
                    plusContigPredictions.back().addExon(entry);
                    prevTargetKeyPlus = targetKey;
                } else {
                    if ((isFirstIterationMinus == true) || (prevTargetKeyMinus != targetKey)) {
                        minusContigPredictions.emplace_back(Prediction());
                        minusContigPredictions.back().setByDPRes(entry);
                        isFirstIterationMinus = false;
                    }
                    // add the exon key
                    minusContigPredictions.back().addExon(entry);
                    prevTargetKeyMinus = targetKey;
                }

                results = Util::skipLine(results);
            }

            // finished collecting all preds from current contig
            clusterPredictions(plusContigPredictions, plusContigRepPreds);
            excludeSameStrandOverlaps(plusContigRepPreds);

            clusterPredictions(minusContigPredictions, minusContigRepPreds);
            excludeSameStrandOverlaps(minusContigRepPreds);

            // write clusters
            writePredsClusters(plusContigPredictions, predictionBuffer, writerRepToMembers, thread_idx);
            writePredsClusters(minusContigPredictions, predictionBuffer, writerRepToMembers, thread_idx);

            // join representatives from both strands and sort by targetKey to comply with expectd order of DP format
            plusContigRepPreds.insert(plusContigRepPreds.end(), minusContigRepPreds.begin(), minusContigRepPreds.end());
            std::stable_sort(plusContigRepPreds.begin(), plusContigRepPreds.end(), Prediction::comparePredictionsByTarget);
            writeRepPredsInDPFormat(plusContigRepPreds, predictionBuffer, par.overlapAllowed, writerGroupedPredictions, thread_idx);

            // close the contig entry with a null byte
            writerRepToMembers.writeEnd(contigKey, thread_idx);
            writerGroupedPredictions.writeEnd(contigKey, thread_idx);

            // move to another contig:
            plusContigPredictions.clear();
            minusContigPredictions.clear();
            plusContigRepPreds.clear();
            minusContigRepPreds.clear();
        }
    }
    writerRepToMembers.close(true);
    writerGroupedPredictions.close(true);
    return EXIT_SUCCESS;
}
