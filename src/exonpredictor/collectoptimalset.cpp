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

static const int PLUS = 1;
static const int MINUS = -1;

const size_t MINIMAL_INTRON_LENGTH = 15;
int const MAX_AA_OVERLAP = 10;

// simple consts, to be changed in the future:
int const CONST_LEGAL_OVERLAP_PENALTY = -5;
int const GAP_OPEN_PENALTY = -2;
int const GAP_EXTEND_PENALTY = -1;

struct potentialExon {
	// constructor:
	potentialExon(size_t alnScore, int contigStart, int contigEnd, int strand, size_t proteinMatchStart, size_t proteinMatchEnd) :
		_alnScore(alnScore), _contigStart(contigStart), _contigEnd(contigEnd), _strand(strand), _proteinMatchStart(proteinMatchStart), _proteinMatchEnd(proteinMatchEnd) {
	}
	
	// information extracted from MMSeqs2 local alignment:
	size_t _alnScore;
	int _contigStart; // the first nucleotide to participate in the alignment times the strand
	int _contigEnd; // the last nucleotide to participate in the alignment times the strand
	int _strand;
    size_t _proteinMatchStart;
	size_t _proteinMatchEnd;
};


int collectoptimalset(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, 2, true, true);

    DBReader<unsigned int> resultReader(par.db1.c_str(), par.db1Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // TO DO: write results
    //DBWriter resultWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    //resultWriter.open();

// analyze each entry of the result DB, this is a swapped DB
// so the original targets play the role of queries...
// i.e., the results are from protein --> potentialExon
#pragma omp for schedule(dynamic, 100)
    for (size_t id = 0; id < resultReader.getSize(); id++) {
        Debug::printProgress(id);

        // unsigned int proteinID = resultReader.getDbKey(id);

        char *results = resultReader.getData(id);
        
        std::vector<potentialExon> plusStrandPotentialExons;
        plusStrandPotentialExons.reserve(10000);
        std::vector<potentialExon> minusStrandPotentialExons;
        minusStrandPotentialExons.reserve(10000);
        int currContigId = -1;
        bool isFirstIteration = true;
        
        while (*results != '\0') {
            char potentialExonKeyStr[255 + 1];
            Util::parseKey(results, potentialExonKeyStr);
            const unsigned int potentialExonKey = (unsigned int) strtoul(potentialExonKeyStr, NULL, 10);
            char *entry[255];
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

            potentialExon currPotentialExon(potentialExonToProteinAlnScore, potentialExonContigStart, potentialExonContigEnd, potentialExonStrand, proteinMatchStart, proteinMatchEnd);

            if (isFirstIteration) {
                currContigId = potentialExonContigId;
                isFirstIteration = false;
            }

            if (potentialExonContigId != currContigId) {
                // sort and then DP on vectors!
                // ...
                size_t numOnPlus = plusStrandPotentialExons.size();
                size_t numOnMinus = minusStrandPotentialExons.size();

                // empty vectors:
                plusStrandPotentialExons.clear();
                minusStrandPotentialExons.clear();
                currContigId = potentialExonContigId;
            }
            
            // push current potentialExon struct to vector:
            if (potentialExonStrand == PLUS) {
                plusStrandPotentialExons.emplace_back(currPotentialExon);
            } else {
                minusStrandPotentialExons.emplace_back(currPotentialExon);
            }

            results = Util::skipLine(results);
        }

        // sort and then DP on vectors - required for the matches of the last contig against the protein
        // ...
        size_t numOnPlus = plusStrandPotentialExons.size();
        size_t numOnMinus = minusStrandPotentialExons.size();

        bool stam = false;

        // int stopMCount = 0;
        // int mCount = 0;
        // for (size_t seqIdx = 0; seqIdx < stopPositions.size(); seqIdx++){
        //     stopMCount += stopPositions[seqIdx].hasStopM;
        //     mCount += stopPositions[seqIdx].hasM;
        // }
        // if (stopPositions.size() > 1){
        //     const float frequency = static_cast<float>(stopMCount) / static_cast<float>(stopPositions.size());
        //     if (frequency >= threshold){
        //         for (size_t seqIdx = 0; seqIdx < stopPositions.size(); seqIdx++){
        //             int target;

        //             int curVal = stopPositions[seqIdx].mPos;
        //             __atomic_load(&addStopAtPosition[stopPositions[seqIdx].id], &target ,__ATOMIC_RELAXED);
        //             do {
        //                 if (target >= curVal) break;
        //             } while (!__atomic_compare_exchange(&addStopAtPosition[stopPositions[seqIdx].id],  &target,  &curVal , false,  __ATOMIC_RELAXED, __ATOMIC_RELAXED));
        //         }
        //     }
        // }
    }

    // cleanup
    //resultWriter.close(DBReader<unsigned int>::DBTYPE_AA);
    resultReader.close();
    
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}



