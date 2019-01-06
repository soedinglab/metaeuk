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
                int iStrand, int iCombinedNormalizedAlnBitScore, double iCombinedEvalue, unsigned int iNumExons, unsigned int iLowContigCoord, 
                unsigned int iHighContigCoord, char * iExonCharptr) : proteinContigStrandId(iProteinContigStrandId), proteinMMSeqs2Key(iProteinMMSeqs2Key), 
                contigMMSeqs2Key(iContigMMSeqs2Key), strand(iStrand),
                combinedNormalizedAlnBitScore(iCombinedNormalizedAlnBitScore), combinedEvalue(iCombinedEvalue),
                numExons(iNumExons), lowContigCoord(iLowContigCoord), highContigCoord(iHighContigCoord) {
        
        // parsing exon info: 20317*20995*21701*21720*21753*21795\n
        // first replace the "\n" with "\0" to avoid potentially super long strings in split
        char exonIdsCharArr[2001];
        size_t indExonIdsCharArr = 0;
        while (*iExonCharptr != '\n') {
            if (indExonIdsCharArr == 2000) {
                Debug(Debug::ERROR) << "ERROR: exonIdsCharArr is too small to hold the exons.\n";
                EXIT(EXIT_FAILURE);
            }
            exonIdsCharArr[indExonIdsCharArr] = *iExonCharptr;
            iExonCharptr++;
            indExonIdsCharArr++;
        }
        exonIdsCharArr[indExonIdsCharArr] = '\0';

        std::vector<std::string> exonIdsSepStrs = Util::split(exonIdsCharArr, "*");
        for (size_t i = 0; i < exonIdsSepStrs.size(); ++i) {
            exonIds.push_back(Util::fast_atoi<int>(exonIdsSepStrs[i].c_str()));
        }
        // initialize cluster assignment as self (will change later):
        clusterProteinContigStrandId = proteinContigStrandId;
        // initialize cluster no overlap assignment as self (will change later):
        noOverlapProteinContigStrandId = proteinContigStrandId;
    }

    // to allow sorting a vector of predictions by their E-values
    static bool comparePredictionsByEvalue (const prediction & aPrediction, const prediction & anotherPrediction) {
        if(aPrediction.combinedEvalue < anotherPrediction.combinedEvalue)
            return true;
        if(aPrediction.combinedEvalue > anotherPrediction.combinedEvalue)
            return false;
        // the following lines will break even cases in a consistent way
        if(aPrediction.lowContigCoord < anotherPrediction.lowContigCoord)
            return true;
        if(aPrediction.lowContigCoord > anotherPrediction.lowContigCoord)
            return false;
        // this line will not be reached...
        return false;
    }

    unsigned int proteinContigStrandId;
    unsigned int proteinMMSeqs2Key;
    unsigned int contigMMSeqs2Key;
    int strand;
    int combinedNormalizedAlnBitScore;
    double combinedEvalue;
    unsigned int numExons;
    unsigned int lowContigCoord;
    unsigned int highContigCoord;
    std::vector<unsigned int> exonIds;
    unsigned int clusterProteinContigStrandId; // this will hold the cluster identifier
    unsigned int noOverlapProteinContigStrandId; // this will hold the cluster identifier for the no overlap data structure
};


int grouppredictions(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, 3, true, true);

    // Terminology: CS = Contig & Strand combination, TCS = Target & Contig & Strand (i.e. -  a prediction identifier)
    // db1 = a map from CS to all its TCS predictions, ordered by their start on the contig, reverse subordered by # exons 
    DBReader<unsigned int> contigStrandSortedMap(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    contigStrandSortedMap.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // db2 = output, cluster format: TCS representative to all TCSs that share an exon with it
    DBWriter writerGroupedPredictions(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_CLUSTER_RES);
    writerGroupedPredictions.open();

    // db3 = output, cluster format: TCS representative to all TCS representatives that overlap it
    DBWriter writerGroupedPredictionsNoOverlap(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_CLUSTER_RES);
    writerGroupedPredictionsNoOverlap.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        // per thread variables
        std::vector<prediction> predictionToCluster;
        predictionToCluster.reserve(EXPECTED_NUM_PREDICTIONS);
        std::vector<prediction> representativePredictions;
        representativePredictions.reserve(EXPECTED_NUM_PREDICTIONS);
        char *entry[255];
        char TCSKeyBuff[128];
        std::string clusterBuffer;
        clusterBuffer.reserve(10000);  
        
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < contigStrandSortedMap.getSize(); id++) {
            Debug::printProgress(id);

            // these will serve to verify sorted order:
            unsigned int prevLowCoord = 0;
            unsigned int prevNumExons = 0;
            int prevBitScore = 0;

            char *contigStrandSortedRecord = contigStrandSortedMap.getData(id, thread_idx);
            while (*contigStrandSortedRecord != '\0') {
                const size_t contigStrandSortedColumns = Util::getWordsOfLine(contigStrandSortedRecord, entry, 255);

                if (contigStrandSortedColumns != 10) {
                    Debug(Debug::ERROR) << "ERROR: the sorted contig+strand map record should contain 10 columns: proteinContigStrandId, proteinMMSeqs2Key, contigMMSeqs2Key, strand, combinedBitScore, combinedEvalue, numExons, lowContigCoord, highContigCoord, exonIDsStr. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int proteinContigStrandId = Util::fast_atoi<int>(entry[0]);
                unsigned int proteinMMSeqs2Key = Util::fast_atoi<int>(entry[1]);
                unsigned int contigMMSeqs2Key = Util::fast_atoi<int>(entry[2]);
                int strand = Util::fast_atoi<int>(entry[3]);
                int combinedNormalizedAlnBitScore = Util::fast_atoi<int>(entry[4]);
                double combinedEvalue = atof(entry[5]);
                unsigned int numExons = Util::fast_atoi<int>(entry[6]);
                unsigned int lowContigCoord = Util::fast_atoi<int>(entry[7]);
                unsigned int highContigCoord = Util::fast_atoi<int>(entry[8]);
                char * exonCharptr = entry[9];

                // verify sorted order:
                if (prevNumExons == 0) {
                    // first iteration:
                    prevNumExons = numExons;
                    prevLowCoord = lowContigCoord;
                    prevBitScore = combinedNormalizedAlnBitScore;
                }
                else if (prevLowCoord > lowContigCoord) {
                    Debug(Debug::ERROR) << "ERROR: Predictions are assumed to be sorted by their start position. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }
                else if ((prevLowCoord == lowContigCoord) && (prevNumExons < numExons)) {
                    Debug(Debug::ERROR) << "ERROR: Predictions are assumed to be reverse sub-sorted by their number of exons. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                } 
                else if ((prevLowCoord == lowContigCoord) && (prevNumExons == numExons) && (prevBitScore < combinedNormalizedAlnBitScore)) {
                    Debug(Debug::ERROR) << "ERROR: Predictions are assumed to be reverse sub-sorted by their bitscore. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }
                else {
                    prevNumExons = numExons;
                    prevLowCoord = lowContigCoord;
                    prevBitScore = combinedNormalizedAlnBitScore;
                }

                predictionToCluster.emplace_back(proteinContigStrandId, proteinMMSeqs2Key, contigMMSeqs2Key, strand, combinedNormalizedAlnBitScore, combinedEvalue, numExons, lowContigCoord, highContigCoord, exonCharptr);
                contigStrandSortedRecord = Util::skipLine(contigStrandSortedRecord);
            }

            // by this stage we have collected all TCS predictions into a vector
            // the index i iterates over cluster representatives: 
            for (size_t i = 0; i < predictionToCluster.size(); ++i) {
                // if i is already assigned - skip it, it is not a cluster representative
                if (predictionToCluster[i].proteinContigStrandId != predictionToCluster[i].clusterProteinContigStrandId) {
                    continue;
                }
                
                // initialize the new cluster:
                char *tmpBuff = Itoa::i32toa_sse2(static_cast<uint32_t>(predictionToCluster[i].proteinContigStrandId), TCSKeyBuff);
                clusterBuffer.append(TCSKeyBuff, tmpBuff - TCSKeyBuff - 1);
                clusterBuffer.append(1, '\n');
                // push the representative into the representative vector - will be used to avoid overlaps
                representativePredictions.emplace_back(predictionToCluster[i]);

                // collect cluster members:
                for (size_t j = (i + 1); j < predictionToCluster.size(); ++j) {
                    if (predictionToCluster[j].lowContigCoord >= predictionToCluster[i].highContigCoord) {
                        // overlap is over - no need to compare to other j - write and move to the next i
                        goto writeClusterTCSi;
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

                        tmpBuff = Itoa::i32toa_sse2(static_cast<uint32_t>(predictionToCluster[j].proteinContigStrandId), TCSKeyBuff);
                        clusterBuffer.append(TCSKeyBuff, tmpBuff - TCSKeyBuff - 1);
                        clusterBuffer.append(1, '\n');
                    }
                }

                writeClusterTCSi:
                writerGroupedPredictions.writeData(clusterBuffer.c_str(), clusterBuffer.size(), predictionToCluster[i].proteinContigStrandId, thread_idx);
                clusterBuffer.clear();
            }

            // output representatives with no overlap in coords:
            // sort vector by E-value:
            std::stable_sort(representativePredictions.begin(), representativePredictions.end(), prediction::comparePredictionsByEvalue);
            // go over sorted representatives - unassigned representatives are non-overlap representatives
            for (size_t i = 0; i < representativePredictions.size(); ++i) {
                // if i is already assigned - skip it, it is not a non-overlap cluster representative
                if (representativePredictions[i].proteinContigStrandId != representativePredictions[i].noOverlapProteinContigStrandId) {
                    continue;
                }
                // initialize the new cluster:
                char *tmpBuff = Itoa::i32toa_sse2(static_cast<uint32_t>(representativePredictions[i].proteinContigStrandId), TCSKeyBuff);
                clusterBuffer.append(TCSKeyBuff, tmpBuff - TCSKeyBuff - 1);
                clusterBuffer.append(1, '\n');

                // collect cluster members:
                for (size_t j = (i + 1); j < representativePredictions.size(); ++j) {
                    if (
                        ((representativePredictions[j].lowContigCoord < representativePredictions[i].lowContigCoord) && 
                        (representativePredictions[j].highContigCoord > representativePredictions[i].lowContigCoord)) || 
                        ((representativePredictions[j].lowContigCoord < representativePredictions[i].highContigCoord) && 
                        (representativePredictions[j].highContigCoord > representativePredictions[i].highContigCoord))) {
                            // i begins in the middle of j or j begins in the middle of i ==> overlap
                            representativePredictions[j].noOverlapProteinContigStrandId = representativePredictions[i].proteinContigStrandId;
                            tmpBuff = Itoa::i32toa_sse2(static_cast<uint32_t>(representativePredictions[j].proteinContigStrandId), TCSKeyBuff);
                            clusterBuffer.append(TCSKeyBuff, tmpBuff - TCSKeyBuff - 1);
                            clusterBuffer.append(1, '\n');
                    }
                }

                writerGroupedPredictionsNoOverlap.writeData(clusterBuffer.c_str(), clusterBuffer.size(), representativePredictions[i].proteinContigStrandId, thread_idx);
                clusterBuffer.clear();
            }
            
            // empty predictions vectors so the same thread can process another CS combination:
            predictionToCluster.clear();
            representativePredictions.clear();
        }        
    }
   
    // cleanup
    writerGroupedPredictions.close();
    writerGroupedPredictionsNoOverlap.close();
    contigStrandSortedMap.close();
    
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}
