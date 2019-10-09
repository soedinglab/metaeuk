#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include "itoa.h"
#include "PredictionParser.h"

#include <limits>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm> 
#include <sstream>

#ifdef OPENMP
#include <omp.h>
#endif

struct dpMatrixRow {
    // constructor
    dpMatrixRow(size_t iPrevPotentialExonId, int iPathScore, size_t iNumExonsInPath) : 
        prevPotentialExonId(iPrevPotentialExonId), pathScore(iPathScore), numExonsInPath(iNumExonsInPath) {

    }
    // the prevPotentialExonId refers to the row Id (i.e., the sorted order)
    size_t prevPotentialExonId;
    int pathScore;
    size_t numExonsInPath;
};

bool isPairCompatible(const PotentialExon & firstPotentialExonOnContig, const PotentialExon & secondPotentialExonOnContig, const size_t minIntronLength, const size_t maxIntronLength, const size_t maxAaOvelap) {
	// it is assumed firstPotentialExonOnContig comes before secondPotentialExonOnContig on contig, i.e.:
    // firstPotentialExonOnContig.contigStart <= secondPotentialExonOnContig.contigStart
    // because negative coordinates are used for the minus strand, the logic works

    // check same strand:
    if (firstPotentialExonOnContig.strand != secondPotentialExonOnContig.strand) {
        return false;
    }
    
    // check the first one does not contain the second one:
    if (secondPotentialExonOnContig.contigEnd < firstPotentialExonOnContig.contigEnd) {
        return false;
    }

    // check overlap on contig:
    int diffOnContig = secondPotentialExonOnContig.contigStart - firstPotentialExonOnContig.contigEnd - 1;
    if (diffOnContig < 0) {
        // overlap on the contig is not allowed:
        return false;
    }
    // if no contig overlap - check legal intron length:
    size_t diffOnContigNonNeg = abs(diffOnContig);
    if ((diffOnContigNonNeg < minIntronLength) || (diffOnContigNonNeg > maxIntronLength)) {
        return false;
    }

    // check overlap on target:
    int diffAAs = secondPotentialExonOnContig.targetMatchStart - firstPotentialExonOnContig.targetMatchEnd - 1;
    // if diffAAs is negative - there is some target overlap:
    if (diffAAs < 0) {
        size_t diffAAsAbs = abs(diffAAs);
        // check overlap is not too long
        if (diffAAsAbs > maxAaOvelap) {
            return false;
        }
    }

    // check contig order is as target order:
    if (secondPotentialExonOnContig.targetMatchStart < firstPotentialExonOnContig.targetMatchStart) {
        return false;
    }

    return true;
}

int getPenaltyForProtCoords(const PotentialExon & prevPotentialExon, const PotentialExon & currPotentialExon, const int setGapOpenPenalty, const int setGapExtendPenalty) {
    // this function is called on a compatible pair
    int diffAAs = currPotentialExon.targetMatchStart - prevPotentialExon.targetMatchEnd - 1;
    if (diffAAs < 0) {
        // legal overlap that should be penalized:
        // by default setGapOpenPenalty = setGapExtendPenalty so this is linear penalty
        int penalty = setGapOpenPenalty + setGapExtendPenalty * (abs(diffAAs) - 1);
        return penalty;
    }
    else if (diffAAs <= 1) {
        // no penalty for missing up to one AA in between exons:
        return 0;
    }
    else {
        // penalize for missing protein fragment:
        // by default setGapOpenPenalty = setGapExtendPenalty so this is linear penalty
        int penalty = setGapOpenPenalty + setGapExtendPenalty * (diffAAs - 1);
        return penalty;
    }

    // never reached:
    return -1000;
}

int findoptimalsetbydp(std::vector<PotentialExon> & potentialExonCandidates, std::vector<PotentialExon> & optimalExonSet, const size_t minIntronLength, const size_t maxIntronLength, const size_t maxAaOvelap, const int setGapOpenPenalty, const int setGapExtendPenalty) {
    size_t numPotentialExonCandidates = potentialExonCandidates.size();
    if (numPotentialExonCandidates == 0) {
        // nothing to do here!
        return (0);
    }

    // sort vector by start on contig:
    std::stable_sort(potentialExonCandidates.begin(), potentialExonCandidates.end(), PotentialExon::comparePotentialExons);

    // prevIdsAndScoresBestPath will hold the DP computation results
    // Each row i represents a potentialExon. They are sorted according to the start on the contig.
    // Each row i is of the struct dpMatrixRow, which works as follows:
    // prevPotentialExonId keeps the id j such that j is the previous potentialExon
    // on the best path ending with the potentialExon i. It will allow for the trace back.
    // pathScore contains the score itself and numExonsInPath contans the number of exons in the path (including i)
    std::vector<dpMatrixRow> prevIdsAndScoresBestPath;
    prevIdsAndScoresBestPath.reserve(numPotentialExonCandidates);
    // initialize:
    for (size_t id = 0; id < numPotentialExonCandidates; ++id) {
        prevIdsAndScoresBestPath.emplace_back(dpMatrixRow(id, potentialExonCandidates[id].bitScore, 1));
    }

    int bestPathScore = 0;
    size_t lastPotentialExonInBestPath = 0;
    // dynamic programming to fill in the matrix, go over all rows - previous values have been computed:
    for (size_t currPotentialExonId = 0; currPotentialExonId < numPotentialExonCandidates; ++currPotentialExonId) {
        for (size_t prevPotentialExonId = 0; prevPotentialExonId < currPotentialExonId; ++prevPotentialExonId) {
            if (isPairCompatible(potentialExonCandidates[prevPotentialExonId], potentialExonCandidates[currPotentialExonId], minIntronLength, maxIntronLength, maxAaOvelap)) {
                int bestScorePathPrevIsLast = prevIdsAndScoresBestPath[prevPotentialExonId].pathScore;
                size_t numExonsInPathPrevIsLast = prevIdsAndScoresBestPath[prevPotentialExonId].numExonsInPath;
                int costOfPrevToCurrTransition = getPenaltyForProtCoords(potentialExonCandidates[prevPotentialExonId], potentialExonCandidates[currPotentialExonId], setGapOpenPenalty, setGapExtendPenalty);
                
                size_t currNumExonsWithPrev = numExonsInPathPrevIsLast + 1;
                int bonusForAddingAnExon = (int) log2(currNumExonsWithPrev); // not the most accurate...
                int currScoreWithPrev = bestScorePathPrevIsLast + costOfPrevToCurrTransition + potentialExonCandidates[currPotentialExonId].bitScore + bonusForAddingAnExon;

                // update row of currPotentialExon in case of improvement:
                if (currScoreWithPrev > prevIdsAndScoresBestPath[currPotentialExonId].pathScore) {
                    prevIdsAndScoresBestPath[currPotentialExonId].prevPotentialExonId = prevPotentialExonId;
                    prevIdsAndScoresBestPath[currPotentialExonId].pathScore = currScoreWithPrev;
                    prevIdsAndScoresBestPath[currPotentialExonId].numExonsInPath = currNumExonsWithPrev;
                }
            }
        }

        // update the global max in case of improvement:
        if (prevIdsAndScoresBestPath[currPotentialExonId].pathScore > bestPathScore) {
            lastPotentialExonInBestPath = currPotentialExonId;
            bestPathScore = prevIdsAndScoresBestPath[currPotentialExonId].pathScore;
        }
    }

    // traceback:
    size_t currExonId = lastPotentialExonInBestPath;
    while (prevIdsAndScoresBestPath[currExonId].prevPotentialExonId != currExonId) {
        optimalExonSet.emplace_back(potentialExonCandidates[currExonId]);
        currExonId = prevIdsAndScoresBestPath[currExonId].prevPotentialExonId;
    }
    // include in the optimal set
    optimalExonSet.emplace_back(potentialExonCandidates[currExonId]);

    // after the traceback, the first exon is in the last place in the vector
    std::reverse(optimalExonSet.begin(), optimalExonSet.end()); 

    return (bestPathScore);
}

int collectoptimalset(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, true, 0, 0);

    if (par.minExonAaLength < par.maxAaOverlap) {
        Debug(Debug::ERROR) << "minExonAaLength was set to be smaller than maxAaOverlap. This can cause trouble for very short exons...\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> resultPerContigReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultPerContigReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> targetsData(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    targetsData.open(DBReader<unsigned int>::NOSORT);
    // get number of AAs in target DB for an E-Value computation
    size_t totNumOfAAsInTargetDb = targetsData.getAminoAcidDBSize(); // method now returns db size for proteins and for profiles by checking dbtype
    targetsData.close();
    double dMetaeukEvalueThr = (double)par.metaeukEvalueThr; // converting to double for precise comparisons

    DBWriter predWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    predWriter.open();

    Debug::Progress progress(resultPerContigReader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<PotentialExon> plusStrandPotentialExons;
        plusStrandPotentialExons.reserve(10000);
        std::vector<PotentialExon> minusStrandPotentialExons;
        minusStrandPotentialExons.reserve(10000);
        std::vector<PotentialExon> plusStrandOptimalExonSet;
        plusStrandOptimalExonSet.reserve(100);
        std::vector<PotentialExon> minusStrandOptimalExonSet;
        minusStrandOptimalExonSet.reserve(100);

        const char *entry[255];

        // each exon line within a prediction has 17 columns
        char exonLineBuffer[2048];
        // this buffer will hold a single prediction with all its exons
        std::string predictionBuffer;
        predictionBuffer.reserve(10000);


#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < resultPerContigReader.getSize(); id++) {
            progress.updateProgress();

            unsigned int contigKey = resultPerContigReader.getDbKey(id);

            char *results = resultPerContigReader.getData(id, thread_idx);
            
            unsigned int currTargetKey = 0;
            bool isFirstIteration = true;

            // keep track of offset when a contig starts
            predWriter.writeStart(thread_idx);

            // process a specific contig
            while (*results != '\0') {
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // each line is a concatentaion of two alignemnts: target<-->potentialExon and potentialExon<-->contig
                if (columns != 20) {
                    Debug(Debug::ERROR) << "there should be 20 columns in the input file. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                PotentialExon currExon;
                currExon.setByAln(entry);

                unsigned int targetKey = currExon.targetKey;
                
                if (isFirstIteration) {
                    currTargetKey = targetKey;
                    isFirstIteration = false;
                }

                // after collecting all the exons for the current target - find optimal set on each strand:
                if (targetKey != currTargetKey) {
                    if (targetKey < currTargetKey) {
                        Debug(Debug::ERROR) << "the targets are assumed to be sorted in increasing order. This doesn't seem to be the case.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    // sort + dynamic programming to find the optimals set:
                    int totalBitScorePlus = findoptimalsetbydp(plusStrandPotentialExons, plusStrandOptimalExonSet, par.minIntronLength, par.maxIntronLength, par.maxAaOverlap, par.setGapOpenPenalty, par.setGapExtendPenalty);
                    int totalBitScoreMinus = findoptimalsetbydp(minusStrandPotentialExons, minusStrandOptimalExonSet, par.minIntronLength, par.maxIntronLength, par.maxAaOverlap, par.setGapOpenPenalty, par.setGapExtendPenalty);

                    // write optimal sets to result file:
                    if (plusStrandOptimalExonSet.size() > 0) {
                        // compute E-Values of the optimal set:
                        // Evalue = m X n * 2^(-S), where m = totNumOfAAsInTargetDb, n = twoStrands, S = combinedNormalizedAlnBitScore
                        double log2EvaluePlus = log2(totNumOfAAsInTargetDb) + log2(2) - totalBitScorePlus;
                        double combinedEvaluePlus = pow(2, log2EvaluePlus);
                        if (combinedEvaluePlus <= dMetaeukEvalueThr) {
                            Prediction predToWrite(currTargetKey, PLUS, totalBitScorePlus, combinedEvaluePlus, plusStrandOptimalExonSet);
                            Prediction::predictionToBuffer(predictionBuffer, exonLineBuffer, predToWrite);
                            predWriter.writeAdd(predictionBuffer.c_str(), predictionBuffer.size(), thread_idx);
                            predictionBuffer.clear();
                        }
                    }
                    if (minusStrandOptimalExonSet.size() > 0) {
                        // compute E-Values of the optimal set:
                        // Evalue = m X n * 2^(-S), where m = totNumOfAAsInTargetDb, n = twoStrands, S = combinedNormalizedAlnBitScore
                        double log2EvalueMinus = log2(totNumOfAAsInTargetDb) + log2(2) - totalBitScoreMinus;
                        double combinedEvalueMinus = pow(2, log2EvalueMinus);
                        if (combinedEvalueMinus <= dMetaeukEvalueThr) {
                            Prediction predToWrite(currTargetKey, MINUS, totalBitScoreMinus, combinedEvalueMinus, minusStrandOptimalExonSet);
                            Prediction::predictionToBuffer(predictionBuffer, exonLineBuffer, predToWrite);
                            predWriter.writeAdd(predictionBuffer.c_str(), predictionBuffer.size(), thread_idx);
                            predictionBuffer.clear();
                        }
                    }

                    // empty vectors between targets:
                    plusStrandPotentialExons.clear();
                    minusStrandPotentialExons.clear();
                    plusStrandOptimalExonSet.clear();
                    minusStrandOptimalExonSet.clear();
                    currTargetKey = targetKey;
                }
                
                // push current potentialExon to vector:
                size_t potentialExonAALen = std::abs(currExon.nucleotideLen) / 3;
                if (potentialExonAALen >= par.minExonAaLength) {
                    if (currExon.strand == PLUS) {
                        plusStrandPotentialExons.emplace_back(currExon);
                    } else {
                        minusStrandPotentialExons.emplace_back(currExon);
                    }
                }
                results = Util::skipLine(results);
            }

            // one last time - required for the matches of the contig against the last target
            // sort + dynamic programming to find the optimals set:
            int totalBitScorePlus = findoptimalsetbydp(plusStrandPotentialExons, plusStrandOptimalExonSet, par.minIntronLength, par.maxIntronLength, par.maxAaOverlap, par.setGapOpenPenalty, par.setGapExtendPenalty);
            int totalBitScoreMinus = findoptimalsetbydp(minusStrandPotentialExons, minusStrandOptimalExonSet, par.minIntronLength, par.maxIntronLength, par.maxAaOverlap, par.setGapOpenPenalty, par.setGapExtendPenalty);
            
            // write optimal sets to result file:
            if (plusStrandOptimalExonSet.size() > 0) {
                // compute E-Values of the optimal set:
                // Evalue = m X n * 2^(-S), where m = totNumOfAAsInTargetDb, n = twoStrands, S = combinedNormalizedAlnBitScore
                double log2EvaluePlus = log2(totNumOfAAsInTargetDb) + log2(2) - totalBitScorePlus;
                double combinedEvaluePlus = pow(2, log2EvaluePlus);
                if (combinedEvaluePlus <= dMetaeukEvalueThr) {
                    Prediction predToWrite(currTargetKey, PLUS, totalBitScorePlus, combinedEvaluePlus, plusStrandOptimalExonSet);
                    Prediction::predictionToBuffer(predictionBuffer, exonLineBuffer, predToWrite);
                    predWriter.writeAdd(predictionBuffer.c_str(), predictionBuffer.size(), thread_idx);
                    predictionBuffer.clear();
                }
            }
            if (minusStrandOptimalExonSet.size() > 0) {
                // compute E-Values of the optimal set:
                // Evalue = m X n * 2^(-S), where m = totNumOfAAsInTargetDb, n = twoStrands, S = combinedNormalizedAlnBitScore
                double log2EvalueMinus = log2(totNumOfAAsInTargetDb) + log2(2) - totalBitScoreMinus;
                double combinedEvalueMinus = pow(2, log2EvalueMinus);
                if (combinedEvalueMinus <= dMetaeukEvalueThr) {
                    Prediction predToWrite(currTargetKey, MINUS, totalBitScoreMinus, combinedEvalueMinus, minusStrandOptimalExonSet);
                    Prediction::predictionToBuffer(predictionBuffer, exonLineBuffer, predToWrite);
                    predWriter.writeAdd(predictionBuffer.c_str(), predictionBuffer.size(), thread_idx);
                    predictionBuffer.clear();
                }
            }
            // close the contig entry with a null byte
            predWriter.writeEnd(contigKey, thread_idx);

            // empty vectors between contigs:
            plusStrandPotentialExons.clear();
            minusStrandPotentialExons.clear();
            plusStrandOptimalExonSet.clear();
            minusStrandOptimalExonSet.clear();
        }
    }

    predWriter.close();
    resultPerContigReader.close();
    return EXIT_SUCCESS;
}