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

//note: N->N, S->S, W->W, U->A, T->A
static const char* iupacReverseComplementTable =
"................................................................"
".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
"................................................................"
"................................................................";

inline char complement(const char c)
{
    return iupacReverseComplementTable[static_cast<unsigned char>(c)];
}

void reverseComplement (const std::string & seq, std::string & revCompSeq) {
    size_t seqLength = revCompSeq.size();
    for(size_t i = 0; i < seqLength; ++i) {
        char currNuc = seq[seqLength - i - 1];
        revCompSeq[i] = complement(currNuc);
    }
}

int unitesetstosequencedb(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, 4, true, true);

    // db1 = contigsDB (data + header):
    std::string contigsDBFilename = par.db1;
    std::string contigsDBIndexFilename = par.db1Index;
    std::string contigsDBHeaderFilename(contigsDBFilename);
    contigsDBHeaderFilename.append("_h");
    std::string contigsDBHeaderIndexFilename(contigsDBFilename);
    contigsDBHeaderIndexFilename.append("_h.index");

    DBReader<unsigned int> contigsData(contigsDBFilename.c_str(), contigsDBIndexFilename.c_str());
    contigsData.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> contigsHeaders(contigsDBHeaderFilename.c_str(), contigsDBHeaderIndexFilename.c_str());
    contigsHeaders.open(DBReader<unsigned int>::NOSORT);

    // db2 = proteinsDB (only the header)
    std::string proteinsDBFilename = par.db2;
    std::string proteinsDBHeaderFilename(proteinsDBFilename);
    proteinsDBHeaderFilename.append("_h");
    std::string proteinsDBHeaderIndexFilename(proteinsDBFilename);
    proteinsDBHeaderIndexFilename.append("_h.index");

    DBReader<unsigned int> proteinsHeaders(proteinsDBHeaderFilename.c_str(), proteinsDBHeaderIndexFilename.c_str());
    proteinsHeaders.open(DBReader<unsigned int>::NOSORT);
    
    // db3 = optimalExonsResuls (exons + map)
    std::string combinedResultBasename = par.db3;
    std::string optimalSetsExonRecordsFilename(combinedResultBasename);
    optimalSetsExonRecordsFilename.append("_optimal_exon_sets");
    std::string optimalSetsExonRecordsIndexFilename(optimalSetsExonRecordsFilename);
    optimalSetsExonRecordsIndexFilename.append(".index");

    std::string setsMapFilename(combinedResultBasename);
    setsMapFilename.append("_protein_contig_strand_map");
    std::string setsMapIndexFilename(setsMapFilename);
    setsMapIndexFilename.append(".index");

    DBReader<unsigned int> optimalSetsExonRecords(optimalSetsExonRecordsFilename.c_str(), optimalSetsExonRecordsIndexFilename.c_str());
    optimalSetsExonRecords.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> setMap(setsMapFilename.c_str(), setsMapIndexFilename.c_str());
    setMap.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    
    // db4 = output
    std::string outDBFilename(par.db4);
    std::string outDBIndexFilename(par.db4);
    outDBIndexFilename.append(".index");

    std::string outDBHeaderFilename(par.db4);
    outDBHeaderFilename.append("_h");
    std::string outDBHeaderIndexFilename(par.db4);
    outDBHeaderIndexFilename.append("_h.index");

    DBWriter concatenatedSetsHeaders(outDBHeaderFilename.c_str(), outDBHeaderIndexFilename.c_str(), par.threads);
    concatenatedSetsHeaders.open();

    DBWriter concatenatedSetsData(outDBFilename.c_str(), outDBIndexFilename.c_str(), par.threads);
    concatenatedSetsData.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        // per thread variables
        std::ostringstream joinedHeaderStream;
        std::ostringstream joinedExonsStream;
        std::string littleStringBuffer;
        littleStringBuffer.reserve(1000);

        char *entry[255];      
        
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < optimalSetsExonRecords.getSize(); id++) {
            Debug::printProgress(id);

            unsigned int setCombinationKey = optimalSetsExonRecords.getDbKey(id);

            // take the record corresponding to the same id from both result files:
            char *optimalExonRecord = optimalSetsExonRecords.getData(id);
            char *setRecord = setMap.getData(id);
            
            // parse setRecord (a single line)
            const size_t setColumns = Util::getWordsOfLine(setRecord, entry, 255);
            if (setColumns != 4) {
                Debug(Debug::ERROR) << "ERROR: the map record should contain 4 columns: proteinMMSeqs2Key, contigMMSeqs2Key, strand, combinedEvalue. This doesn't seem to be the case.\n";
                EXIT(EXIT_FAILURE);
            }
            unsigned int proteinMMSeqs2Key = Util::fast_atoi<int>(entry[0]);
            unsigned int contigMMSeqs2Key = Util::fast_atoi<int>(entry[1]);
            int strand = Util::fast_atoi<int>(entry[2]);
            double combinedEvalue = atof(entry[3]);
            setRecord = Util::skipLine(setRecord);

            // get non-MMSeqs2 identifiers from header files:
            const char* contigHeader = contigsHeaders.getDataByDBKey(contigMMSeqs2Key);
            const char* proteinHeader = proteinsHeaders.getDataByDBKey(proteinMMSeqs2Key);
            
            // get contig data:
            const char* contigData = contigsData.getDataByDBKey(contigMMSeqs2Key);
            
            // initialize header:
            size_t numWhiteCharsToTrim = 2; // a single space and a \n
            littleStringBuffer.append(proteinHeader, std::strlen(proteinHeader) - numWhiteCharsToTrim);
            joinedHeaderStream << littleStringBuffer << "|";
            littleStringBuffer.clear();
            littleStringBuffer.append(contigHeader, std::strlen(contigHeader) - numWhiteCharsToTrim);
            joinedHeaderStream << littleStringBuffer << "|";
            littleStringBuffer.clear();
            (strand == PLUS) ? joinedHeaderStream << "plus" : joinedHeaderStream << "minus";
            joinedHeaderStream << "|" << combinedEvalue;
            
            int lastTargetPosMatched = 0;
            while (*optimalExonRecord != '\0') {
                const size_t exonRecordColumns = Util::getWordsOfLine(optimalExonRecord, entry, 255);

                if (exonRecordColumns != 10) {
                    Debug(Debug::ERROR) << "ERROR: the exon record should contain 10 columnn (alignment format). This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                // unsigned int exonMMSeqs2Key = Util::fast_atoi<int>(entry[0]);
                // int exonToProteinAlnScore = Util::fast_atoi<int>(entry[1]);
                // double exonSequenceIdentity = atof(entry[2]);
                // double exonEvalue = atof(entry[3]);
                int proteinMatchStart =  Util::fast_atoi<int>(entry[4]);
                int proteinMatchEnd = Util::fast_atoi<int>(entry[5]);
                // int proteinLen = Util::fast_atoi<int>(entry[6]);
                int exonContigStart = Util::fast_atoi<int>(entry[7]);
                int exonContigEnd = Util::fast_atoi<int>(entry[8]);
                int exonNucleotideLen = Util::fast_atoi<int>(entry[9]);

                if (strand == MINUS) {
                    if ((exonContigStart > 0) || (exonContigEnd > 0)) {
                        Debug(Debug::ERROR) << "ERROR: strand is MINUS but the contig coordinates are positive. Soemthing is wrong.\n";
                        EXIT(EXIT_FAILURE);
                    }
                }

                // in order to avoid target overlaps, we remove a few codons from the start of the current exon if needed:
                int exonAdjustedContigStart = exonContigStart;
                int exonAdjustedNucleotideLen = exonNucleotideLen;
                if (lastTargetPosMatched > proteinMatchStart) {
                    int diffInAAs = lastTargetPosMatched - proteinMatchStart + 1;
                    exonAdjustedContigStart += (3 * diffInAAs * strand);
                    exonAdjustedNucleotideLen -= (3 * diffInAAs);
                }
                int lowContigCoord = (strand == PLUS) ? exonAdjustedContigStart : (-1 * exonContigEnd);
                // extract the segment from the contig:
                std::string exonContigSeq(&contigData[lowContigCoord], (size_t)exonAdjustedNucleotideLen);
                // update the last AA of the protein that was matched:
                lastTargetPosMatched = proteinMatchEnd;

                // write the header and data to streams:
                joinedHeaderStream << "|" << abs(exonContigStart) << "[" << abs(exonAdjustedContigStart) << "]:";
                joinedHeaderStream << abs(exonContigEnd) << "[" << abs(exonContigEnd) << "]:";
                joinedHeaderStream << exonNucleotideLen << "[" << exonAdjustedNucleotideLen << "]";

                if (strand == PLUS) {
                    joinedExonsStream << exonContigSeq;
                } else {
                    std::string exonContigRevCompSeq(exonContigSeq);
                    reverseComplement (exonContigSeq, exonContigRevCompSeq);
                    joinedExonsStream << exonContigRevCompSeq;
                }

                optimalExonRecord = Util::skipLine(optimalExonRecord);
            }

            // finished collecting all exons for the current combination, write streams:
            joinedHeaderStream << "\n";
            std::string headerToWrite = joinedHeaderStream.str();
            concatenatedSetsHeaders.writeData(headerToWrite.c_str(), headerToWrite.length(), setCombinationKey);
            joinedHeaderStream.str("");

            joinedExonsStream << "\n";
            std::string dataToWrite = joinedExonsStream.str();
            concatenatedSetsData.writeData(dataToWrite.c_str(), dataToWrite.length(), setCombinationKey);
            joinedExonsStream.str("");

        }
    }

    // cleanup
    contigsData.close();
    contigsHeaders.close();
    proteinsHeaders.close();
    optimalSetsExonRecords.close();
    setMap.close();
    concatenatedSetsHeaders.close();
    concatenatedSetsData.close();
    
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}



