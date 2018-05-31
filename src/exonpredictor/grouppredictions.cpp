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
#include <string.h>
#include "itoa.h"

#ifdef OPENMP
#include <omp.h>
#endif

const size_t EXPECTED_NUM_PREDICTIONS = 100000;

struct prediction {
    // constructor
    prediction(unsigned int iProteinContigStrandId, unsigned int iProteinMMSeqs2Key, unsigned int iContigMMSeqs2Key, 
                int iStrand, int iCombinedNormalizedAlnBitScore, unsigned int iNumExons, unsigned int iLowContigCoord, 
                unsigned int iHighContigCoord, char * iExonCharptr) : proteinContigStrandId(iProteinContigStrandId), proteinMMSeqs2Key(iProteinMMSeqs2Key), 
                contigMMSeqs2Key(iContigMMSeqs2Key), strand(iStrand),
                combinedNormalizedAlnBitScore(iCombinedNormalizedAlnBitScore), numExons(iNumExons),
                lowContigCoord(iLowContigCoord), highContigCoord(iHighContigCoord) {
        // parsing exon info: 20317*20995*21701*21720*21753*21795\n
        std::vector<std::string> exonIdsStr = Util::split(iExonCharptr, "\n");
        std::vector<std::string> exonIdsSepStrs = Util::split(exonIdsStr[0].c_str(), "*");
        for (size_t i = 0; i < exonIdsSepStrs.size(); ++i) {
            exonIds.push_back(Util::fast_atoi<int>(exonIdsSepStrs[i].c_str()));
        }

        // initialize the cluster assignment as self (will change later):
        clusterProteinContigStrandId = proteinContigStrandId;
    }

    unsigned int proteinContigStrandId;
    unsigned int proteinMMSeqs2Key;
    unsigned int contigMMSeqs2Key;
    int strand;
    int combinedNormalizedAlnBitScore;
    unsigned int numExons;
    unsigned int lowContigCoord;
    unsigned int highContigCoord;
    std::vector<unsigned int> exonIds;
    unsigned int clusterProteinContigStrandId; // this will hold the cluster identifier
};


int grouppredictions(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, 2, true, true);

    // Terminology: CS = Contig & Strand combination, TCS = Target & Contig & Strand (i.e. -  a prediction identifier)
    // db1 = a map from CS to all its TCS predictions, ordered by their start on the contig, reverse subordered by # exons 
    DBReader<unsigned int> contigStrandSortedMap(par.db1.c_str(), par.db1Index.c_str());
    contigStrandSortedMap.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // db2 = output, cluster format
    DBWriter writerGroupedPredictions(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    writerGroupedPredictions.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        // per thread variables
        std::vector<prediction> predictionToCluster;
        predictionToCluster.reserve(EXPECTED_NUM_PREDICTIONS);

        char *entry[255];       
        
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < contigStrandSortedMap.getSize(); id++) {
            Debug::printProgress(id);

            // these will serve to verify sorted order:
            unsigned int prevLowCoord = 0;
            unsigned int prevNumExons = 0;

            char *contigStrandSortedRecord = contigStrandSortedMap.getData(id);
            while (*contigStrandSortedRecord != '\0') {
                const size_t contigStrandSortedColumns = Util::getWordsOfLine(contigStrandSortedRecord, entry, 255);

                if (contigStrandSortedColumns != 9) {
                    Debug(Debug::ERROR) << "ERROR: the sorted contig+strand map record should contain 9 columns: proteinContigStrandId, proteinMMSeqs2Key, contigMMSeqs2Key, strand, combinedBitScore, numExons, lowContigCoord, highContigCoord, exonIDsStr. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int proteinContigStrandId = Util::fast_atoi<int>(entry[0]);
                unsigned int proteinMMSeqs2Key = Util::fast_atoi<int>(entry[1]);
                unsigned int contigMMSeqs2Key = Util::fast_atoi<int>(entry[2]);
                int strand = Util::fast_atoi<int>(entry[3]);
                int combinedNormalizedAlnBitScore = Util::fast_atoi<int>(entry[4]);
                unsigned int numExons = Util::fast_atoi<int>(entry[5]);
                unsigned int lowContigCoord = Util::fast_atoi<int>(entry[6]);
                unsigned int highContigCoord = Util::fast_atoi<int>(entry[7]);
                char * exonCharptr = entry[8];

                // verify sorted order:
                if (prevNumExons == 0) {
                    // first iteration:
                    prevNumExons = numExons;
                    prevLowCoord = lowContigCoord;
                }
                else if (prevLowCoord > lowContigCoord) {
                    Debug(Debug::ERROR) << "ERROR: Predictions are assumed to be sorted by their start position. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }
                else if ((prevLowCoord == lowContigCoord) && (prevNumExons < numExons)) {
                    Debug(Debug::ERROR) << "ERROR: Predictions are assumed to be reverse sub-sorted by their number of exons. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }
                else {
                    prevNumExons = numExons;
                    prevLowCoord = lowContigCoord;
                }

                predictionToCluster.emplace_back(proteinContigStrandId, proteinMMSeqs2Key, contigMMSeqs2Key, strand, combinedNormalizedAlnBitScore, numExons, lowContigCoord, highContigCoord, exonCharptr);

                contigStrandSortedRecord = Util::skipLine(contigStrandSortedRecord);
            }

            // by this stage we have collected all TCS predictions into a vector
            char clusterBuffer[10000];
            size_t i = 0;
            size_t next_i = 0;
            while (i < predictionToCluster.size()) {
                next_i = 0;
                for (size_t j = (i + 1); j < predictionToCluster.size(); ++j) {
                    if (predictionToCluster[j].lowContigCoord >= predictionToCluster[i].highContigCoord) {
                        // overlap is over - no need to compare to any more j - moving to the next i
                        if (next_i == 0) {
                            next_i = j;
                        }
                        break;
                    }
                    

                    bool doIandJshareAnExon = false;
                    for (size_t exon_id_i = 0; exon_id_i < predictionToCluster[i].exonIds.size(); ++exon_id_i) {
                        for (size_t exon_id_j = 0; exon_id_j < predictionToCluster[j].exonIds.size(); ++exon_id_j) {
                            if (predictionToCluster[i].exonIds[exon_id_i] == predictionToCluster[j].exonIds[exon_id_j]) {
                                doIandJshareAnExon = true;
                                goto endExonComparison;
                            }
                        }
                    }

                    // assign j to the cluster of i if it has a mutual exon and if it wasn't already assigned:
                    endExonComparison:
                    if (doIandJshareAnExon && (predictionToCluster[j].proteinContigStrandId == predictionToCluster[j].clusterProteinContigStrandId)) {
                        predictionToCluster[j].clusterProteinContigStrandId = predictionToCluster[i].proteinContigStrandId;
                    }

                    if ((! doIandJshareAnExon) && (next_i == 0)) {
                        next_i = j;
                    }
                }
            }
        }
    }

    
    // cleanup
    writerGroupedPredictions.close();
    contigStrandSortedMap.close();
    
    Debug(Debug::INFO) << "\nDone.\n";


    return EXIT_SUCCESS;
}



