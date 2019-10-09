#ifndef PREDICTION_PARSER_H
#define PREDICTION_PARSER_H

#include "LocalParameters.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "MathUtil.h"
#include "itoa.h"

const int PLUS = 1;
const int MINUS = -1;

struct PotentialExon {
    void setByAln(const char ** exonData) {
        // assumption - exonData has 20 columns!

        // from T<-->O alignment:
        targetKey = Util::fast_atoi<int>(exonData[0]);
        bitScore = Util::fast_atoi<int>(exonData[1]);
        seqId = strtod(exonData[2],NULL);
        evalue = strtod(exonData[3],NULL);

        int orfProtStart = Util::fast_atoi<int>(exonData[4]);
        int orfProtEnd = Util::fast_atoi<int>(exonData[5]);

        targetMatchStart =  Util::fast_atoi<int>(exonData[7]);
        targetMatchEnd = Util::fast_atoi<int>(exonData[8]);
        targetLen = Util::fast_atoi<int>(exonData[9]);

        // from O<-->C alignment:
        exonKey = Util::fast_atoi<int>(exonData[10]);

        int potentialExonContigStartBeforeTrim = Util::fast_atoi<int>(exonData[17]);
        int potentialExonContigEndBeforeTrim = Util::fast_atoi<int>(exonData[18]);
        
        // adjust strand and contig positions based on match to T
        // plus strand:
        if (potentialExonContigStartBeforeTrim < potentialExonContigEndBeforeTrim) {
            contigStart = potentialExonContigStartBeforeTrim + (orfProtStart * 3);
            contigEnd = potentialExonContigStartBeforeTrim + (orfProtEnd * 3) + 2;
            strand = PLUS;
        }
        // minus strand:
        else {
            // multiplying by minus allows carrying out same logic as with plus strand:
            contigStart = -1 * (potentialExonContigStartBeforeTrim - (orfProtStart * 3));
            contigEnd = -1 * (potentialExonContigStartBeforeTrim - (orfProtEnd * 3) - 2);
            strand = MINUS;
        }

        // compute length in nucleotides
        nucleotideLen = contigEnd - contigStart + 1;
        if (nucleotideLen % 3 != 0) {
            Debug(Debug::ERROR) << "seems like the coordiantes do not dictate a legal length for a codon segment.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    void setByDPRes (const char ** exonData) {
        // assumption - exonData has 17 columns!

        // from predictions part (same for all exons in prediction):
        targetKey = Util::fast_atoi<int>(exonData[0]);
        strand = Util::fast_atoi<int>(exonData[1]);

        // from exon part:
        exonKey = Util::fast_atoi<int>(exonData[7]);
        bitScore = Util::fast_atoi<int>(exonData[8]);
        seqId = strtod(exonData[9],NULL);
        evalue = strtod(exonData[10],NULL);

        targetMatchStart = Util::fast_atoi<int>(exonData[11]);
        targetMatchEnd = Util::fast_atoi<int>(exonData[12]);
        targetLen = Util::fast_atoi<int>(exonData[13]);

        contigStart = Util::fast_atoi<int>(exonData[14]);
        contigEnd = Util::fast_atoi<int>(exonData[15]);
        nucleotideLen = Util::fast_atoi<int>(exonData[16]);
    }

    static size_t exonToBuffer (char * exonBuffer, const PotentialExon & exon) {
        // writes 10 columns format
        char * basePos = exonBuffer;
        
        char * tmpBuff = Itoa::u32toa_sse2((uint32_t) exon.exonKey, exonBuffer);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(exon.bitScore, tmpBuff);
        *(tmpBuff-1) = '\t';

        float seqIdFlt = exon.seqId;
        if (seqIdFlt == 1.0) {
            *(tmpBuff) = '1';
            tmpBuff++;
            *(tmpBuff) = '.';
            tmpBuff++;
            *(tmpBuff) = '0';
            tmpBuff++;
            *(tmpBuff) = '0';
            tmpBuff++;
            *(tmpBuff) = '0';
            tmpBuff++;
            *(tmpBuff) = '\t';
            tmpBuff++;
        } else {
            *(tmpBuff) = '0';
            tmpBuff++;
            *(tmpBuff) = '.';
            tmpBuff++;
            if (seqIdFlt < 0.10) {
                *(tmpBuff) = '0';
                tmpBuff++;
            }
            if (seqIdFlt < 0.01) {
                *(tmpBuff) = '0';
                tmpBuff++;
            }
            int seqId = seqIdFlt*1000;
            tmpBuff = Itoa::i32toa_sse2(seqId, tmpBuff);
            *(tmpBuff-1) = '\t';
        }

        tmpBuff += sprintf(tmpBuff,"%.3E",exon.evalue);
        tmpBuff++;
        *(tmpBuff-1) = '\t';

        tmpBuff = Itoa::i32toa_sse2(exon.targetMatchStart, tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(exon.targetMatchEnd, tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(exon.targetLen, tmpBuff);
        *(tmpBuff-1) = '\t';

        tmpBuff = Itoa::i32toa_sse2(exon.contigStart, tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(exon.contigEnd, tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(exon.nucleotideLen, tmpBuff);

        *(tmpBuff-1) = '\n';
        *(tmpBuff) = '\0';
        return (tmpBuff - basePos);
    }

    // allow comparing PotentialExons by their start on the contig
    static bool comparePotentialExons (const PotentialExon & aPotentialExon, const PotentialExon & anotherPotentialExon) {
        if(aPotentialExon.contigStart < anotherPotentialExon.contigStart)
            return true;
        if(aPotentialExon.contigStart > anotherPotentialExon.contigStart)
            return false;
        // the following lines will break even cases in a consistent way
        if(aPotentialExon.contigEnd < anotherPotentialExon.contigEnd)
            return true;
        if(aPotentialExon.contigEnd > anotherPotentialExon.contigEnd)
            return false;
        // if this line is reached, it is the same PotentialExon (same start & same end)...
        return false;
    }

    // contig start and end refer to the first (and last) nucleotides to participate in the alignment
    // the coordinates are with respect to the contig start (5', plus strand) and are negative
    // in case of the minus strand. This way, on both strands, start < end.
    unsigned int exonKey;
    unsigned int targetKey;
    int strand;

    unsigned int bitScore;
    double seqId;
    double evalue;

    int targetMatchStart;
    int targetMatchEnd;
    int targetLen;

    int contigStart;
    int contigEnd;
    int nucleotideLen;
};

class Prediction {
    public:
    Prediction() {};

    Prediction(const unsigned int itargetKey, const int istrand, const unsigned int itotalBitscore, 
                const double icombinedEvalue, std::vector<PotentialExon> & ioptimalExonSet) {
        targetKey = itargetKey;
        strand = istrand;
        totalBitscore = itotalBitscore;
        combinedEvalue = icombinedEvalue;
        numExons = ioptimalExonSet.size();
        optimalExonSet = ioptimalExonSet;

        PotentialExon firstExon = optimalExonSet[0];
        PotentialExon lastExon = optimalExonSet[numExons - 1];

        // since contigStart and contigEnd are negative on the MINUS strand, we multiply by (-1)
        // to assure ContigCoords are always positive
        lowContigCoord = (firstExon.strand == PLUS) ? firstExon.contigStart : (-1 * lastExon.contigEnd);
        highContigCoord = (firstExon.strand == PLUS) ? lastExon.contigEnd : (-1 * firstExon.contigStart);

        // initialize cluster assignment:
        isClustered = false;
        clusterId = 0;

        // initialize cluster no overlap assignment:
        isNoOverlapClustered = false;
        noOverlapClusterId = 0;
    }

    void setByDPRes (const char** entry) {
        targetKey = Util::fast_atoi<int>(entry[0]);
        strand = Util::fast_atoi<int>(entry[1]);
        totalBitscore = Util::fast_atoi<int>(entry[2]);
        combinedEvalue = atof(entry[3]);
        numExons = Util::fast_atoi<int>(entry[4]);
        lowContigCoord = Util::fast_atoi<int>(entry[5]);
        highContigCoord = Util::fast_atoi<int>(entry[6]);

        // initialize cluster assignment:
        isClustered = false;
        clusterId = 0;

        // initialize cluster no overlap assignment:
        isNoOverlapClustered = false;
        noOverlapClusterId = 0;
    }

    void clearPred () {
        targetKey = 0;
        strand = 0;
        totalBitscore = 0;
        combinedEvalue = 0;
        numExons = 0;
        lowContigCoord = 0;
        highContigCoord = 0;

        // clear exon vector
        optimalExonSet.clear();

        // initialize cluster assignment:
        isClustered = false;
        clusterId = 0;

        // initialize cluster no overlap assignment:
        isNoOverlapClustered = false;
        noOverlapClusterId = 0;
    }

    void addExon (const char ** exonData) {
        PotentialExon exon;
        exon.setByDPRes(exonData);
        optimalExonSet.emplace_back(exon);
    }

    static int getTargetKey (const char** entry) {
        return (Util::fast_atoi<int>(entry[0]));
    }

    static int getStrand (const char** entry) {
        return (Util::fast_atoi<int>(entry[1]));
    }

    // to allow sorting a vector of predictions by their start on the contig
    static bool comparePredictionsByContigStart (const Prediction & aPrediction, const Prediction & anotherPrediction) {
        if(aPrediction.lowContigCoord < anotherPrediction.lowContigCoord)
            return true;
        if(aPrediction.lowContigCoord > anotherPrediction.lowContigCoord)
            return false;
        // the following lines will break even cases in a consistent way (longer comes before shorter)
        if(aPrediction.highContigCoord > anotherPrediction.highContigCoord)
            return true;
        if(aPrediction.highContigCoord < anotherPrediction.highContigCoord)
            return false;
        // the following lines will break even cases in a consistent way (higher bitscore comes before lower)
        if(aPrediction.totalBitscore > anotherPrediction.totalBitscore)
            return true;
        if(aPrediction.totalBitscore < anotherPrediction.totalBitscore)
            return false;
        // the following lines will break even cases in a consistent way (higher bitscore comes before lower)
        if(aPrediction.targetKey < anotherPrediction.targetKey)
            return true;
        if(aPrediction.targetKey > anotherPrediction.targetKey)
            return false;
        // this line will not be reached...
        return false;
    }

    // to allow sorting a vector of predictions by their E-values
    static bool comparePredictionsByEvalue (const Prediction & aPrediction, const Prediction & anotherPrediction) {
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

    // to allow sorting a vector of contig predictions by their targetKey and strand
    static bool comparePredictionsByTarget (const Prediction & aPrediction, const Prediction & anotherPrediction) {
        if(aPrediction.targetKey < anotherPrediction.targetKey)
            return true;
        if(aPrediction.targetKey > anotherPrediction.targetKey)
            return false;
        // the following lines will break even cases in a consistent way (MINUS before PLUS)
        if(aPrediction.strand < anotherPrediction.strand)
            return true;
        if(aPrediction.strand > anotherPrediction.strand)
            return false;
        // this line will not be reached...
        return false;
    }

    static void predictionToBuffer (std::string& predictionBuffer, char* exonBuffer, const Prediction & prediction) {
        for (size_t i = 0; i < prediction.optimalExonSet.size(); ++i) {
            char* tmpBuff = exonBuffer;
            // add the columns that are joint for all exons
            tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(prediction.targetKey), tmpBuff);
            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::i32toa_sse2(static_cast<uint32_t>(prediction.strand), tmpBuff);
            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(prediction.totalBitscore), tmpBuff);
            *(tmpBuff-1) = '\t';
            tmpBuff += sprintf(tmpBuff, "%.3E", prediction.combinedEvalue);
            tmpBuff++;
            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(prediction.numExons), tmpBuff);
            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(prediction.lowContigCoord), tmpBuff);
            *(tmpBuff-1) = '\t';
            tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(prediction.highContigCoord), tmpBuff);
            *(tmpBuff-1) = '\t';

            // add exon information
            size_t len = PotentialExon::exonToBuffer(tmpBuff, prediction.optimalExonSet[i]); 
            tmpBuff += len;

            // add a new line after each exon
            *(tmpBuff-1) = '\n';
            predictionBuffer.append(exonBuffer, tmpBuff - exonBuffer);
        }
    }

    static size_t predictionClusterToBuffer (char * clusterBuffer, const Prediction & prediction) {
        // write: Representative(T,S) , Member(T,S)
        char * basePos = clusterBuffer;
        char * tmpBuff = basePos;

        // clusterId is the TargetKey of the representative. The representative is on the same strand
        tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(prediction.clusterId), tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(static_cast<uint32_t>(prediction.strand), tmpBuff);
        *(tmpBuff-1) = '\t';

        tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(prediction.targetKey), tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(static_cast<uint32_t>(prediction.strand), tmpBuff);
        *(tmpBuff-1) = '\n';

        // close buffer
        *(tmpBuff) = '\0';
        return (tmpBuff - basePos);
    }

    // members
    unsigned int targetKey;
    int strand;
    unsigned int totalBitscore;
    double combinedEvalue;
    unsigned int numExons;
    unsigned int lowContigCoord;
    unsigned int highContigCoord;
    std::vector<PotentialExon> optimalExonSet;

    // members for grouping
    bool isClustered;
    unsigned int clusterId; // will hold the target key of the final representative

    bool isNoOverlapClustered; // will allow excluding same-strand overlaps
    unsigned int noOverlapClusterId; // will hold the target key of the final representative after resolving overlaps
};

#endif // PREDICTION_PARSER_H