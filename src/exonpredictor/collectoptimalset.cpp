#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>

#include "LocalParameters.h"
#include "Matcher.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include <limits>
#include <cstdint>
#include <queue>
#include "itoa.h"

#ifdef OPENMP
#include <omp.h>
#endif

const int PLUS = 1;
const int MINUS = -1;

struct potentialExon {
    // constructor
    potentialExon(unsigned int iMMSeqs2Key, int iAlnScore, int iContigStart, int iContigEnd, int iStrand, int iTargetMatchStart, int iTargetMatchEnd, int iTargetLen, float iPotentialExonSequenceIdentity, double iPotentialExonEval, const double scoreBias) :
        MMSeqs2Key(iMMSeqs2Key), alnScore(iAlnScore), contigStart(iContigStart), contigEnd(iContigEnd), strand(iStrand), targetMatchStart(iTargetMatchStart), targetMatchEnd(iTargetMatchEnd) {
            // update the result_t object:
            float proteinCover = float(targetMatchEnd - targetMatchStart + 1) / iTargetLen;
            int potentialExonLengthInNucleotides = iContigEnd - iContigStart + 1;
            int alignemntLength = potentialExonLengthInNucleotides / 3;
            if ((potentialExonLengthInNucleotides % 3 != 0) || ((alignemntLength * 3) != potentialExonLengthInNucleotides)) {
                std::string potentialExonStr = potentialExonToStr();
                Debug(Debug::ERROR) << "seems like the coordiantes do not dictate a legal length for a codon segment. " << potentialExonStr << "\n";
                EXIT(EXIT_FAILURE);
            }

            // since the search was swapped, the "db" is the potentialExon and the "q" is the target
            potentialExonAlignemntRes.score      = iAlnScore;
            potentialExonAlignemntRes.alnLength  = alignemntLength;
            if (!(MathUtil::AreSame(0,scoreBias))) {
                double addedScoreBias = alignemntLength * scoreBias;
                double alnScoreNoBias = iAlnScore - addedScoreBias;
                potentialExonAlignemntRes.score = (alnScoreNoBias < 0.0) ? int(alnScoreNoBias - 0.5) : int(alnScoreNoBias + 0.5);
            }
            potentialExonAlignemntRes.qcov       = proteinCover;
            potentialExonAlignemntRes.dbcov      = 1.0;
            potentialExonAlignemntRes.seqId      = iPotentialExonSequenceIdentity;
            potentialExonAlignemntRes.eval       = iPotentialExonEval;
            

            potentialExonAlignemntRes.qStartPos  = iTargetMatchStart;
            potentialExonAlignemntRes.qEndPos    = iTargetMatchEnd;
            potentialExonAlignemntRes.qLen       = iTargetLen;

            potentialExonAlignemntRes.dbKey      = iMMSeqs2Key;
            potentialExonAlignemntRes.dbStartPos = iContigStart;
            potentialExonAlignemntRes.dbEndPos   = iContigEnd;
            potentialExonAlignemntRes.dbLen      = potentialExonLengthInNucleotides;
            
            potentialExonAlignemntRes.backtrace  = "";
    }

    // information extracted from MMSeqs2 local alignment
    unsigned int MMSeqs2Key;
    int alnScore;
    // contig start and end refer to the first (and last) nucleotides to participate in the alignment
    // the coordinates are with respect to the contig start (5', plus strand) and are negative
    // in case of the minus strand. This way, in both strands, start < end.
    int contigStart;
    int contigEnd;
    int strand;
    int targetMatchStart;
    int targetMatchEnd;

    // will assist in printing to result file:
    Matcher::result_t potentialExonAlignemntRes;

    // to allow sorting a vector of potentialExon by their start on the contig
    static bool comparePotentialExons (const potentialExon & aPotentialExon, const potentialExon & anotherPotentialExon) {
        if(aPotentialExon.contigStart < anotherPotentialExon.contigStart)
            return true;
        if(aPotentialExon.contigStart > anotherPotentialExon.contigStart)
            return false;
        // the following lines will break even cases in a consistent way
        if(aPotentialExon.contigEnd < anotherPotentialExon.contigEnd)
            return true;
        if(aPotentialExon.contigEnd > anotherPotentialExon.contigEnd)
            return false;
        // if this line is reached, it is the same potentialExon (same start & same end)...
        return false;
    }

    std::string potentialExonToStr() {
        std::stringstream ss;
        ss << "MMSeqs2Key: " << MMSeqs2Key << ", alnScore: " << alnScore << ", contigStart: " << contigStart << ", contigEnd: " << contigEnd << ", targetMatchStart: " << targetMatchStart << ", targetMatchEnd: " << targetMatchEnd;
        return (ss.str());
    }
};

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

bool isPairCompatible(const potentialExon & firstPotentialExonOnContig, const potentialExon & secondPotentialExonOnContig, const size_t minIntronLength, const size_t maxIntronLength, const size_t maxAaOvelap) {
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

int getPenaltyForProtCoords(const potentialExon & prevPotentialExon, const potentialExon & currPotentialExon, const int setGapOpenPenalty, const int setGapExtendPenalty) {
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

int findoptimalsetbydp(std::vector<potentialExon> & potentialExonCandidates, std::vector<potentialExon> & optimalExonSet, const size_t minIntronLength, const size_t maxIntronLength, const size_t maxAaOvelap, const int setGapOpenPenalty, const int setGapExtendPenalty) {
    size_t numPotentialExonCandidates = potentialExonCandidates.size();
    if (numPotentialExonCandidates == 0) {
        // nothing to do here!
        return (0);
    }

    // sort vector by start on contig:
    std::stable_sort(potentialExonCandidates.begin(), potentialExonCandidates.end(), potentialExon::comparePotentialExons);

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
        prevIdsAndScoresBestPath.emplace_back(dpMatrixRow(id, potentialExonCandidates[id].alnScore, 1));
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
                int currScoreWithPrev = bestScorePathPrevIsLast + costOfPrevToCurrTransition + potentialExonCandidates[currPotentialExonId].alnScore + bonusForAddingAnExon;

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

    // trace back:
    size_t currExonId = lastPotentialExonInBestPath;
    while (prevIdsAndScoresBestPath[currExonId].prevPotentialExonId != currExonId) {
        optimalExonSet.emplace_back(potentialExonCandidates[currExonId]);
        currExonId = prevIdsAndScoresBestPath[currExonId].prevPotentialExonId;
    }
    // include in the optimal set
    optimalExonSet.emplace_back(potentialExonCandidates[currExonId]);

    return (bestPathScore);
}

void getOptimalSetContigCoords (std::vector<potentialExon> & optimalExonSet, int & lowContigCoord, int & highContigCoord) {
    int numExonsInOptimalSet = optimalExonSet.size();
    potentialExon firstExon = optimalExonSet[numExonsInOptimalSet - 1]; // the first exon is in the last place in the vector
    potentialExon lastExon = optimalExonSet[0]; // the last exon is in the 0 place in the vector

    // since contigStart and contigEnd are negative on the MINUS strand, we multiply by (-1)
    // to assure ContigCoords are always positive
    lowContigCoord = (firstExon.strand == PLUS) ? firstExon.contigStart : (-1 * lastExon.contigEnd);
    highContigCoord = (firstExon.strand == PLUS) ? lastExon.contigEnd : (-1 * firstExon.contigStart);
}

size_t fillPredictionBuffer (char * predictionBuffer, unsigned int targetKey, int strand, int totalBitScore, double combinedEvalue, std::vector<potentialExon> & optimalExonSet) {
    // the buffer contains one or more lines
    // each line starts with the prediciotns fields: targetKey, contigKey, strand, sumBitScore, ...
    // then each line contains a single exon information
    int lowContigCoord;
    int highContigCoord;
    getOptimalSetContigCoords (optimalExonSet, lowContigCoord, highContigCoord);
    size_t numExons = optimalExonSet.size();

    char * basePos = predictionBuffer;
    char * tmpBuff = basePos;

    // go over vector in reverse order (the last exon is in place 0 in the vector)
    for (int i = (numExons - 1); i >= 0; --i) {
        // add the columns that are joint for all exons
        tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(targetKey), tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(static_cast<uint32_t>(strand), tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(totalBitScore), tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff += sprintf(tmpBuff, "%.3E", combinedEvalue);
        tmpBuff++;
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(numExons), tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(lowContigCoord), tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(highContigCoord), tmpBuff);
        *(tmpBuff-1) = '\t';

        // add the exon information
        size_t len = Matcher::resultToBuffer(tmpBuff, optimalExonSet[i].potentialExonAlignemntRes, false);
        tmpBuff += len;

        // add a new line after each exon
        *(tmpBuff-1) = '\n';
    }

    // close the buffer
    *(tmpBuff) = '\0';
    return (tmpBuff - basePos);
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
        std::vector<potentialExon> plusStrandPotentialExons;
        plusStrandPotentialExons.reserve(10000);
        std::vector<potentialExon> minusStrandPotentialExons;
        minusStrandPotentialExons.reserve(10000);
        std::vector<potentialExon> plusStrandOptimalExonSet;
        plusStrandOptimalExonSet.reserve(100);
        std::vector<potentialExon> minusStrandOptimalExonSet;
        minusStrandOptimalExonSet.reserve(100);

        const char *entry[255];
        
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < resultPerContigReader.getSize(); id++) {
            progress.updateProgress();

            unsigned int contigKey = resultPerContigReader.getDbKey(id);

            char *results = resultPerContigReader.getData(id, thread_idx);
            
            unsigned int currTargetKey = 0;
            bool isFirstIteration = true;

            // this buffer will hold a single prediction with all its exons
            // 20 columns * num_exons = 40 * 100
            char predictionBuffer[5000]; 

            // keep track of offset when a contig starts
            predWriter.writeStart(thread_idx);

            // process a specific contig
            while (*results != '\0') {
                const size_t columns = Util::getWordsOfLine(results, entry, 255);

                // each line is a concatentaion of two alignemnts
                // the first alignemnt is with respect to target<-->potentialExon
                // the second alignemnt is with respect to potentialExon<-->contig (constructed to contain potentialExonId)
                size_t firstColumnOfSecondAlignemnt = 10;
                if (columns != 20) {
                    Debug(Debug::ERROR) << "there should be 20 columns in the input file. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int targetKey = Util::fast_atoi<int>(entry[0]);
                int potentialExonToTargetAlnScore = Util::fast_atoi<int>(entry[1]);
                double potentialExonSequenceIdentity = atof(entry[2]);
                double potentialExonEvalue = atof(entry[3]);
                int potentialExonMatchStart = Util::fast_atoi<int>(entry[4]);
                int potentialExonMatchEnd = Util::fast_atoi<int>(entry[5]);
                int targetMatchStart = Util::fast_atoi<int>(entry[7]);
                int targetMatchEnd = Util::fast_atoi<int>(entry[8]);
                int targetLen = Util::fast_atoi<int>(entry[9]);

                // take relavant info from potentialExon<-->contig alignment:
                // the potentialExonId is there thanks to a hack by resultstocontig.cpp
                unsigned int potentialExonId = Util::fast_atoi<int>(entry[firstColumnOfSecondAlignemnt]);
                int potentialExonContigStartBeforeTrim = Util::fast_atoi<int>(entry[firstColumnOfSecondAlignemnt + 7]);
                int potentialExonContigEndBeforeTrim = Util::fast_atoi<int>(entry[firstColumnOfSecondAlignemnt + 8]);

                // "trim" the potentialExon, i.e., adjust start and end on the contig:
                int potentialExonContigStart;
                int potentialExonContigEnd;
                int potentialExonStrand;

                // plus strand:
                if (potentialExonContigStartBeforeTrim < potentialExonContigEndBeforeTrim) {
                    potentialExonContigStart = potentialExonContigStartBeforeTrim + (potentialExonMatchStart * 3);
                    potentialExonContigEnd = potentialExonContigStartBeforeTrim + (potentialExonMatchEnd * 3) + 2;
                    potentialExonStrand = PLUS;
                }
                // minus strand:
                else {
                    // multiplying by minus allows carrying out same logic as with plus strand:
                    potentialExonContigStart = -1 * (potentialExonContigStartBeforeTrim - (potentialExonMatchStart * 3));
                    potentialExonContigEnd = -1 * (potentialExonContigStartBeforeTrim - (potentialExonMatchEnd * 3) - 2);
                    potentialExonStrand = MINUS;
                }

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
                            size_t mapCombinationLen = fillPredictionBuffer(predictionBuffer, currTargetKey, PLUS, totalBitScorePlus, combinedEvaluePlus, plusStrandOptimalExonSet);
                            predWriter.writeAdd(predictionBuffer, mapCombinationLen, thread_idx);
                        }
                    }
                    if (minusStrandOptimalExonSet.size() > 0) {
                        // compute E-Values of the optimal set:
                        // Evalue = m X n * 2^(-S), where m = totNumOfAAsInTargetDb, n = twoStrands, S = combinedNormalizedAlnBitScore
                        double log2EvalueMinus = log2(totNumOfAAsInTargetDb) + log2(2) - totalBitScoreMinus;
                        double combinedEvalueMinus = pow(2, log2EvalueMinus);
                        if (combinedEvalueMinus <= dMetaeukEvalueThr) {
                            size_t mapCombinationLen = fillPredictionBuffer(predictionBuffer, currTargetKey, MINUS, totalBitScoreMinus, combinedEvalueMinus, minusStrandOptimalExonSet);
                            predWriter.writeAdd(predictionBuffer, mapCombinationLen, thread_idx);
                        }
                    }

                    // empty vectors between targets:
                    plusStrandPotentialExons.clear();
                    minusStrandPotentialExons.clear();
                    plusStrandOptimalExonSet.clear();
                    minusStrandOptimalExonSet.clear();
                    currTargetKey = targetKey;
                }
                
                // push current potentialExon struct to vector:
                size_t potentialExonAALen = (std::abs(potentialExonContigEnd - potentialExonContigStart) + 1) / 3;
                if (potentialExonAALen >= par.minExonAaLength) {
                    if (potentialExonStrand == PLUS) {
                        plusStrandPotentialExons.emplace_back(potentialExonId, potentialExonToTargetAlnScore, potentialExonContigStart, potentialExonContigEnd, potentialExonStrand, targetMatchStart, targetMatchEnd, targetLen, potentialExonSequenceIdentity, potentialExonEvalue, par.scoreBias);
                    } else {
                        minusStrandPotentialExons.emplace_back(potentialExonId, potentialExonToTargetAlnScore, potentialExonContigStart, potentialExonContigEnd, potentialExonStrand, targetMatchStart, targetMatchEnd, targetLen, potentialExonSequenceIdentity, potentialExonEvalue, par.scoreBias);
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
                    size_t mapCombinationLen = fillPredictionBuffer(predictionBuffer, currTargetKey, PLUS, totalBitScorePlus, combinedEvaluePlus, plusStrandOptimalExonSet);
                    predWriter.writeAdd(predictionBuffer, mapCombinationLen, thread_idx);
                }
            }
            if (minusStrandOptimalExonSet.size() > 0) {
                // compute E-Values of the optimal set:
                // Evalue = m X n * 2^(-S), where m = totNumOfAAsInTargetDb, n = twoStrands, S = combinedNormalizedAlnBitScore
                double log2EvalueMinus = log2(totNumOfAAsInTargetDb) + log2(2) - totalBitScoreMinus;
                double combinedEvalueMinus = pow(2, log2EvalueMinus);
                if (combinedEvalueMinus <= dMetaeukEvalueThr) {
                    size_t mapCombinationLen = fillPredictionBuffer(predictionBuffer, currTargetKey, MINUS, totalBitScoreMinus, combinedEvalueMinus, minusStrandOptimalExonSet);
                    predWriter.writeAdd(predictionBuffer, mapCombinationLen, thread_idx);
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

    predWriter.close(true);
    resultPerContigReader.close();

    return EXIT_SUCCESS;
}