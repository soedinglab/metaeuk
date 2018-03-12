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

#ifdef OPENMP
#include <omp.h>
#endif

const int PLUS = 1;
const int MINUS = -1;

const int MINIMAL_INTRON_LENGTH = 15;
int const MAX_AA_OVERLAP = 10;

// simple consts, to be changed in the future:
int const CONST_LEGAL_OVERLAP_PENALTY = -5;
int const GAP_OPEN_PENALTY = -2;
int const GAP_EXTEND_PENALTY = -1;

struct potentialExon {
	// constructor
	potentialExon(int iMMSeqs2Key, int iAlnScore, int iContigStart, int iContigEnd, int iStrand, int iProteinMatchStart, int iProteinMatchEnd) :
		MMSeqs2Key(iMMSeqs2Key), alnScore(iAlnScore), contigStart(iContigStart), contigEnd(iContigEnd), strand(iStrand), proteinMatchStart(iProteinMatchStart), proteinMatchEnd(iProteinMatchEnd) {
	}
	
	// information extracted from MMSeqs2 local alignment
    int MMSeqs2Key;
    int alnScore;
	int contigStart; // the first nucleotide to participate in the alignment times the strand
	int contigEnd; // the last nucleotide to participate in the alignment times the strand
	int strand;
    int proteinMatchStart;
	int proteinMatchEnd;

    // define operator to allow sorting a vector of potentialExon by their start on the contig
    bool operator < (const potentialExon & anotherPotentialExon) const {
        return (contigStart < anotherPotentialExon.contigStart);
    }

    std::string potentialExonToStr() {
        std::stringstream ss;
        ss << "MMSeqs2Key: " << MMSeqs2Key << ", alnScore: " << alnScore << ", contigStart: " << contigStart << ", contigEnd: " << contigEnd << ", proteinMatchStart: " << proteinMatchStart << ", proteinMatchEnd: " << proteinMatchEnd;
        return (ss.str());
    }
};

bool isPairCompatible(const potentialExon & firstPotentialExonOnContig, const potentialExon & secondPotentialExonOnContig) {
	// check same strand:
	if (firstPotentialExonOnContig.strand != secondPotentialExonOnContig.strand) {
		return false;
	}
    
    // check one does not contain the other:
	if (secondPotentialExonOnContig.contigEnd < firstPotentialExonOnContig.contigEnd) {
		return false;
	}

	// check gap/overlap on contig:
	int diffOnContig = secondPotentialExonOnContig.contigStart - firstPotentialExonOnContig.contigEnd - 1;
	if (diffOnContig < MINIMAL_INTRON_LENGTH) {
		return false;
	}

	// check gap/overlap on target (also that contig order is as target order):
	int numAAsInGap = secondPotentialExonOnContig.proteinMatchStart - firstPotentialExonOnContig.proteinMatchEnd - 1;
	if (numAAsInGap < -(MAX_AA_OVERLAP)) {
		return false;
	}

	return true;
}

int getPenaltyForProtCoords(const potentialExon & prevPotentialExon, const potentialExon & currPotentialExon) {
	int numAAsInGap = currPotentialExon.proteinMatchStart - prevPotentialExon.proteinMatchEnd - 1;
	if (numAAsInGap < 0) {
		// legal overlap that should be penalized:
		return CONST_LEGAL_OVERLAP_PENALTY;
	}
	else if (numAAsInGap <= 1) {
		// no penalty for missing up to one AA in between fragments:
		return 0;
	}
	else {
        // penalize for missing protein fragment:
		int penalty = GAP_OPEN_PENALTY + GAP_EXTEND_PENALTY * (numAAsInGap - 1);
		return penalty;
	}

	// never reached:
	return 1;
}

void findoptimalsetbydp(std::vector<potentialExon> & potentialExonCandidates, std::vector<potentialExon> & optimalExonSet) {   
    size_t numPotentialExonCandidates = potentialExonCandidates.size();
    if (numPotentialExonCandidates == 0) {
        return;
    }

    // sort vectors by start on contig:
    std::sort(potentialExonCandidates.begin(), potentialExonCandidates.end());

    // prevIdsAndScoresBestPath will hold the DP computation results
	// the first column will be used for backtracing and works as follows:
	// for entry i it keeps the id j such that j is the previous potentialExon
	// on the best path ending with the potentialExon i
	// the second column contains the score itself
    // "id" below refers to the order of the start on the contig. It will allow for the trace back.
    std::vector<std::vector<int>> prevIdsAndScoresBestPath;
    prevIdsAndScoresBestPath.reserve(numPotentialExonCandidates);
	// initialize:
	for (size_t id = 0; id < numPotentialExonCandidates; ++id) {
        std::vector<int> rowOfId;
        rowOfId.emplace_back(id);
        rowOfId.emplace_back(potentialExonCandidates[id].alnScore);
		prevIdsAndScoresBestPath.emplace_back(rowOfId);
	}

    int bestPathScore = 0;
	int lastPotentialExonInBestPath = 0;
	// dynamic programming to fill in the matrix, go over all rows - previous values have been computed:
    for (size_t currPotentialExonId = 0; currPotentialExonId < numPotentialExonCandidates; ++currPotentialExonId) {
        int currPotentialExonAlnScore = potentialExonCandidates[currPotentialExonId].alnScore;
        // initialize:
		int currBestScore = currPotentialExonAlnScore;
        for (size_t prevPotentialExonId = 0; prevPotentialExonId < currPotentialExonId; ++prevPotentialExonId) {
            if (isPairCompatible(potentialExonCandidates[prevPotentialExonId], potentialExonCandidates[currPotentialExonId])) {
                int bestScorePathPrevIsLast = prevIdsAndScoresBestPath[prevPotentialExonId][1];
				int costOfPrevToCurrTransition = getPenaltyForProtCoords(potentialExonCandidates[prevPotentialExonId], potentialExonCandidates[currPotentialExonId]);
				int currScoreWithPrev = bestScorePathPrevIsLast + costOfPrevToCurrTransition + currPotentialExonAlnScore;

                // update row of current potentialExon in case of improvement:
				if (currScoreWithPrev > currBestScore) {
					currBestScore = currScoreWithPrev;
					prevIdsAndScoresBestPath[currPotentialExonId][0] = prevPotentialExonId;
					prevIdsAndScoresBestPath[currPotentialExonId][1] = currScoreWithPrev;
				}
            }
        }

        // update the global max in case of improvement:
		if (currBestScore > bestPathScore) {
			lastPotentialExonInBestPath = currPotentialExonId;
			bestPathScore = currBestScore;
		}
    }

    // trace back (for now, just with print):
    std::cout << "--- best path (last to first) with score: " << bestPathScore << " ---" << std::endl;
	int currExonId = lastPotentialExonInBestPath;
	while (prevIdsAndScoresBestPath[currExonId][0] != currExonId) {
		std::string currExonStr = potentialExonCandidates[currExonId].potentialExonToStr();
        std::cout << currExonStr << std::endl << "is after:" << std::endl;
		currExonId = prevIdsAndScoresBestPath[currExonId][0];

        // include in the optimal set
        optimalExonSet.emplace_back(potentialExonCandidates[currExonId]);
	}
	std::string currExonStr = potentialExonCandidates[currExonId].potentialExonToStr();
    std::cout << currExonStr << std::endl;
    std::cout << "-------" << std::endl;
    // include in the optimal set
    optimalExonSet.emplace_back(potentialExonCandidates[currExonId]);

}

int collectoptimalset(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, 2, true, true);

    DBReader<unsigned int> resultReader(par.db1.c_str(), par.db1Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    resultWriter.open();

    // analyze each entry of the result DB, this is a swapped DB
    // so the original targets play the role of queries...
    // i.e., the results are from protein --> potentialExon
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

        char potentialExonMMSeqs2KeyStr[255 + 1];
        char *entry[255];           
        
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            Debug::printProgress(id);

            // unsigned int proteinID = resultReader.getDbKey(id);

            char *results = resultReader.getData(id);
            
            int currContigId = -1;
            bool isFirstIteration = true;
            
            while (*results != '\0') {
                Util::parseKey(results, potentialExonMMSeqs2KeyStr);
                const unsigned int potentialExonMMSeqs2Key = (unsigned int) strtoul(potentialExonMMSeqs2KeyStr, NULL, 10);
                const size_t columns = Util::getWordsOfLine(results, entry, 255);

                // the structure of entry is a concatentaion of two alignemnts
                // the first alignemnt is with respect to protein<-->potentialExon
                // the second alignemnt is with respect to potentialExon<-->contig
                // the number of columns can be 20 or 22 (depending on the "-a" option)
                size_t firstColumnOfSecondAlignemnt = 10;
                if (columns == 22) {
                    firstColumnOfSecondAlignemnt = 11;
                }

                // take relavant info from protein<-->potentialExon alignment:
                int potentialExonToProteinAlnScore = Util::fast_atoi<int>(entry[1]);
                int proteinMatchStart =  Util::fast_atoi<int>(entry[4]);
                int proteinMatchEnd = Util::fast_atoi<int>(entry[5]);
                int potentialExonMatchStart = Util::fast_atoi<int>(entry[7]);
                int potentialExonMatchEnd = Util::fast_atoi<int>(entry[8]);

                // take relavant info from potentialExon<-->contig alignment:
                int potentialExonContigId = Util::fast_atoi<int>(entry[firstColumnOfSecondAlignemnt]);
                int potentialExonContigStartBeforeTrim = Util::fast_atoi<int>(entry[firstColumnOfSecondAlignemnt + 7]);
                int potentialExonContigEndBeforeTrim = Util::fast_atoi<int>(entry[firstColumnOfSecondAlignemnt + 8]);

                // "trim" the potentialExon, i.e., adjust start and end on the contig
                int potentialExonContigStart;
                int potentialExonContigEnd;
                int potentialExonStrand;

                // plus strand:
                if (potentialExonContigStartBeforeTrim < potentialExonContigEndBeforeTrim) {
                    potentialExonContigStart = potentialExonContigStartBeforeTrim + (potentialExonMatchStart * 3);
                    potentialExonContigEnd = potentialExonContigStartBeforeTrim + (potentialExonMatchEnd * 3) + 2;
                    potentialExonStrand = PLUS;
                }
                // mimus strand:
                else {
                    // multiplying by minus allows carrying out same logic as with plus strand
                    potentialExonContigStart = -1 * (potentialExonContigStartBeforeTrim - (potentialExonMatchStart * 3));
                    potentialExonContigEnd = -1 * (potentialExonContigStartBeforeTrim - (potentialExonMatchEnd * 3) - 2);
                    potentialExonStrand = MINUS;
                }


                if (isFirstIteration) {
                    currContigId = potentialExonContigId;
                    isFirstIteration = false;
                }

                if (potentialExonContigId != currContigId) {
                    if (potentialExonContigId < currContigId) {
                        Debug(Debug::ERROR) << "ERROR: the contigs are assumed to be sorted in increasing order. This doesn't seem to be the case.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    // sort + dynamic programming to find the optimals set:
                    findoptimalsetbydp(plusStrandPotentialExons, plusStrandOptimalExonSet);
                    findoptimalsetbydp(minusStrandPotentialExons, minusStrandOptimalExonSet);
                    
                    // TO DO: write optimal sets to result file:
                    // ...

                    // empty vectors:
                    plusStrandPotentialExons.clear();
                    minusStrandPotentialExons.clear();
                    plusStrandOptimalExonSet.clear();
                    minusStrandOptimalExonSet.clear();
                    currContigId = potentialExonContigId;
                }
                
                // push current potentialExon struct to vector:
                if (potentialExonStrand == PLUS) {
                    plusStrandPotentialExons.emplace_back(potentialExonMMSeqs2Key, potentialExonToProteinAlnScore, potentialExonContigStart, potentialExonContigEnd, potentialExonStrand, proteinMatchStart, proteinMatchEnd);
                } else {
                    minusStrandPotentialExons.emplace_back(potentialExonMMSeqs2Key, potentialExonToProteinAlnScore, potentialExonContigStart, potentialExonContigEnd, potentialExonStrand, proteinMatchStart, proteinMatchEnd);
                }

                results = Util::skipLine(results);
            }

            // one last time - required for the matches of the last contig against the protein
            // sort + dynamic programming to find the optimals set:
            findoptimalsetbydp(plusStrandPotentialExons, plusStrandOptimalExonSet);
            findoptimalsetbydp(minusStrandPotentialExons, minusStrandOptimalExonSet);
            
            // TO DO: write optimal sets to result file:
            // ...

            // empty vectors:
            plusStrandPotentialExons.clear();
            minusStrandPotentialExons.clear();
            plusStrandOptimalExonSet.clear();
            minusStrandOptimalExonSet.clear();
            
            bool stam = false;

        }
    }

    // cleanup
    resultWriter.close();
    resultReader.close();
    
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}



