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
    dpMatrixRow(size_t iPrevPotentialExonId, int iPathScore, size_t iNumExonsInPath, double iPathTargetCov, int iPathAALen) : 
        prevPotentialExonId(iPrevPotentialExonId), pathScore(iPathScore), numExonsInPath(iNumExonsInPath), 
        pathTargetCov(iPathTargetCov), pathAALen(iPathAALen) {

    }
    // the prevPotentialExonId refers to the row Id (i.e., the sorted order)
    size_t prevPotentialExonId;
    int pathScore;
    size_t numExonsInPath;
    double pathTargetCov;
    int pathAALen;
};

bool isPairCompatible(const PotentialExon & firstPotentialExonOnContig, const PotentialExon & secondPotentialExonOnContig, 
                      const size_t minIntronLength, const size_t maxIntronLength, const size_t maxAaOvelap, size_t & aaOverlapTarget) {
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
    aaOverlapTarget = 0;
    // if diffAAs is negative - there is some target overlap:
    if (diffAAs < 0) {
        aaOverlapTarget = abs(diffAAs);
        // check overlap is not too long
        if (aaOverlapTarget > maxAaOvelap) {
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

int findoptimalsetbydp(std::vector<PotentialExon> & potentialExonCandidates, std::vector<PotentialExon> & optimalExonSet, 
                        const size_t minIntronLength, const size_t maxIntronLength, const size_t maxAaOvelap, const int setGapOpenPenalty, 
                        const int setGapExtendPenalty, const double dMetaeukTargetCovThr) {
    size_t numPotentialExonCandidates = potentialExonCandidates.size();
    if (numPotentialExonCandidates == 0) {
        // nothing to do here!
        return (0);
    }

    // sort vector by not-used(0)/used(1) and then by start on contig:
    std::stable_sort(potentialExonCandidates.begin(), potentialExonCandidates.end(), PotentialExon::comparePotentialExons);

    // find the first used exon (placed at the end)
    size_t firstIndexOfUsed = numPotentialExonCandidates;
    for (size_t i = 0; i < potentialExonCandidates.size(); ++i) {
        if (potentialExonCandidates[i].isUsed == true) {
            firstIndexOfUsed = i;
            break;
        }
    }
    // remove all the used exons
    potentialExonCandidates.resize(firstIndexOfUsed);
    numPotentialExonCandidates = potentialExonCandidates.size();

    // prevIdsAndScoresBestPath will hold the DP computation results
    // Each row i represents a potentialExon. They are sorted according to the start on the contig.
    // Each row i is of the struct dpMatrixRow, which works as follows:
    // prevPotentialExonId keeps the id j such that j is the previous potentialExon
    // on the best path ending with the potentialExon i. It will allow for the trace back.
    // pathScore contains the score itself and numExonsInPath contans the number of exons in the path (including i)
    // pathTargetCov contains the proportion of the target the path covers
    // pathAALen contains the total number of AAs in the path
    int targetLength = potentialExonCandidates[0].targetLen;
    if (targetLength == 0) {
        Debug(Debug::ERROR) << "target length is 0 and this cannot be.\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<dpMatrixRow> prevIdsAndScoresBestPath;
    prevIdsAndScoresBestPath.reserve(numPotentialExonCandidates);
    // initialize:
    for (size_t id = 0; id < numPotentialExonCandidates; ++id) {
        prevIdsAndScoresBestPath.emplace_back(dpMatrixRow(id, potentialExonCandidates[id].bitScore, 1, 
                                              potentialExonCandidates[id].targetCov, potentialExonCandidates[id].aaLen));
        
        // sanity check - all exons refer to the same target
        if (potentialExonCandidates[id].targetLen != targetLength) {
            Debug(Debug::ERROR) << "two exons are analyzed in the context of different targets.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    int bestPathScore = 0;
    size_t lastPotentialExonInBestPath = 0;
    
    // dynamic programming to fill in the matrix, go over all rows - previous values have been computed:
    for (size_t currPotentialExonId = 0; currPotentialExonId < numPotentialExonCandidates; ++currPotentialExonId) {
        for (size_t prevPotentialExonId = 0; prevPotentialExonId < currPotentialExonId; ++prevPotentialExonId) {
            size_t pairAaOverlapTarget = 0;
            if (isPairCompatible(potentialExonCandidates[prevPotentialExonId], potentialExonCandidates[currPotentialExonId], 
                                 minIntronLength, maxIntronLength, maxAaOvelap, pairAaOverlapTarget)) {
                
                int bestScorePathPrevIsLast = prevIdsAndScoresBestPath[prevPotentialExonId].pathScore;
                size_t numExonsInPathPrevIsLast = prevIdsAndScoresBestPath[prevPotentialExonId].numExonsInPath;
                int costOfPrevToCurrTransition = getPenaltyForProtCoords(potentialExonCandidates[prevPotentialExonId], potentialExonCandidates[currPotentialExonId], setGapOpenPenalty, setGapExtendPenalty);
                
                size_t currNumExonsWithPrev = numExonsInPathPrevIsLast + 1;
                int bonusForAddingAnExon = (int) log2(currNumExonsWithPrev); // not the most accurate...
                int currScoreWithPrev = bestScorePathPrevIsLast + costOfPrevToCurrTransition + potentialExonCandidates[currPotentialExonId].bitScore + bonusForAddingAnExon;
                
                // add curr candidate contribution to tcov and path length
                double currCovWithPrev = prevIdsAndScoresBestPath[prevPotentialExonId].pathTargetCov + potentialExonCandidates[currPotentialExonId].targetCov;
                int currAALenWithPrev = prevIdsAndScoresBestPath[prevPotentialExonId].pathAALen + potentialExonCandidates[currPotentialExonId].aaLen - pairAaOverlapTarget;

                // update row of currPotentialExon in case of improvement:
                if (currScoreWithPrev > prevIdsAndScoresBestPath[currPotentialExonId].pathScore) {
                    prevIdsAndScoresBestPath[currPotentialExonId].prevPotentialExonId = prevPotentialExonId;
                    prevIdsAndScoresBestPath[currPotentialExonId].pathScore = currScoreWithPrev;
                    prevIdsAndScoresBestPath[currPotentialExonId].numExonsInPath = currNumExonsWithPrev;
                    prevIdsAndScoresBestPath[currPotentialExonId].pathTargetCov = currCovWithPrev;
                    prevIdsAndScoresBestPath[currPotentialExonId].pathAALen = currAALenWithPrev;
                }
            }
        }

        // update the global max in case of improvement that covers the target:
        //if (prevIdsAndScoresBestPath[currPotentialExonId].pathTargetCov >= dMetaeukTargetCovThr) {
        if ((double)prevIdsAndScoresBestPath[currPotentialExonId].pathAALen / (double)targetLength >= dMetaeukTargetCovThr) {
            if (prevIdsAndScoresBestPath[currPotentialExonId].pathScore > bestPathScore) {
                lastPotentialExonInBestPath = currPotentialExonId;
                bestPathScore = prevIdsAndScoresBestPath[currPotentialExonId].pathScore;
            }
        }
    }

    // bestPathScore is 0 when no path covers enough of the target
    if (bestPathScore == 0) {
        return (0);
    }

    // traceback:
    size_t currExonId = lastPotentialExonInBestPath;
    while (prevIdsAndScoresBestPath[currExonId].prevPotentialExonId != currExonId) {
        optimalExonSet.emplace_back(potentialExonCandidates[currExonId]);
        potentialExonCandidates[currExonId].isUsed = true;
        currExonId = prevIdsAndScoresBestPath[currExonId].prevPotentialExonId;
    }
    // include in the optimal set
    optimalExonSet.emplace_back(potentialExonCandidates[currExonId]);
    potentialExonCandidates[currExonId].isUsed = true;

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
    double dMetaeukTargetCovThr = (double)par.metaeukTargetCovThr;

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
                    
                    size_t numIters = 0;
                    while ((numIters < par.maxExonSets) && (plusStrandPotentialExons.size() > 0 || minusStrandPotentialExons.size() > 0)) {
                        // sort + dynamic programming to find the optimals set:
                        int totalBitScorePlus = findoptimalsetbydp(plusStrandPotentialExons, plusStrandOptimalExonSet, par.minIntronLength, par.maxIntronLength, par.maxAaOverlap, par.setGapOpenPenalty, par.setGapExtendPenalty, dMetaeukTargetCovThr);
                        int totalBitScoreMinus = findoptimalsetbydp(minusStrandPotentialExons, minusStrandOptimalExonSet, par.minIntronLength, par.maxIntronLength, par.maxAaOverlap, par.setGapOpenPenalty, par.setGapExtendPenalty, dMetaeukTargetCovThr);

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

                        // empty optimal sets between writes:
                        plusStrandOptimalExonSet.clear();
                        minusStrandOptimalExonSet.clear();
                        numIters++;
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
            size_t numIters = 0;
            while ((numIters < par.maxExonSets) && (plusStrandPotentialExons.size() > 0 || minusStrandPotentialExons.size() > 0)) {
                    
                // sort + dynamic programming to find the optimals set:
                int totalBitScorePlus = findoptimalsetbydp(plusStrandPotentialExons, plusStrandOptimalExonSet, par.minIntronLength, par.maxIntronLength, par.maxAaOverlap, par.setGapOpenPenalty, par.setGapExtendPenalty, dMetaeukTargetCovThr);
                int totalBitScoreMinus = findoptimalsetbydp(minusStrandPotentialExons, minusStrandOptimalExonSet, par.minIntronLength, par.maxIntronLength, par.maxAaOverlap, par.setGapOpenPenalty, par.setGapExtendPenalty, dMetaeukTargetCovThr);
                
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

                // empty optimal sets between writes:
                plusStrandOptimalExonSet.clear();
                minusStrandOptimalExonSet.clear();
                numIters++;
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
