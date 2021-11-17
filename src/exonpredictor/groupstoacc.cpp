#include "LocalParameters.h"
#include "PredictionParser.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "MathUtil.h"
#include "itoa.h"


#ifdef OPENMP
#include <omp.h>
#endif

int groupstoacc(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, true, 0, 0);

    // db1 = input, contigsDB (only the header is used)
    DBReader<unsigned int> contigsHeaders(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    contigsHeaders.open(DBReader<unsigned int>::NOSORT);

    // db2 = input, targetsDB (only the header is used)
    DBReader<unsigned int> targetsHeaders(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    targetsHeaders.open(DBReader<unsigned int>::NOSORT);

    // db3 = input, grouping of predictions: T,S of representatives to T,S of prediction
    DBReader<unsigned int> readerRepToMembers(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    readerRepToMembers.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_OMIT_FILE);
    writer.open();

    Debug::Progress progress(readerRepToMembers.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        // per thread variables
        std::string strToWrite;
        const char *entry[255];
        
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < readerRepToMembers.getSize(); id++) {
            progress.updateProgress();

            unsigned int contigKey = readerRepToMembers.getDbKey(id);

            char *results = readerRepToMembers.getData(id, thread_idx);

            // get contig header:
            const char* contigHeader = contigsHeaders.getDataByDBKey(contigKey, thread_idx);
            std::string contigHeaderAcc = Util::parseFastaHeader(contigHeader);

            // process a specific contig
            while (*results != '\0') {
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // each line informs of a representative and of a member
                if (columns != 6) {
                    Debug(Debug::ERROR) << "There should be 6 columns in the input file. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int repTargetKey = Util::fast_atoi<int>(entry[0]);
                int repStrand = Util::fast_atoi<int>(entry[1]);
                int repLowCoord = Util::fast_atoi<int>(entry[2]);
                unsigned int memTargetKey = Util::fast_atoi<int>(entry[3]);
                int memStrand = Util::fast_atoi<int>(entry[4]);
                int memLowCoord = Util::fast_atoi<int>(entry[5]);

                if (repStrand != memStrand) {
                    Debug(Debug::ERROR) << "A representative should always be on the same strand as its member. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                std::string strandStr = "+";
                if (repStrand == MINUS) {
                    strandStr = "-";
                }

                const char* repTargetHeader = targetsHeaders.getDataByDBKey(repTargetKey, thread_idx);
                std::string repTargetHeaderAcc = Util::parseFastaHeader(repTargetHeader);

                const char* memTargetHeader = targetsHeaders.getDataByDBKey(memTargetKey, thread_idx);
                std::string memTargetHeaderAcc = Util::parseFastaHeader(memTargetHeader);
                
                strToWrite.clear();
                strToWrite = repTargetHeaderAcc + "|" + contigHeaderAcc + "|" + strandStr + "|" + std::to_string(repLowCoord) + "\t" +
                             memTargetHeaderAcc + "|" + contigHeaderAcc + "|" + strandStr + "|" + std::to_string(memLowCoord) + "\n";

                writer.writeData(strToWrite.c_str(), strToWrite.size(), 0, thread_idx, false, false);

                results = Util::skipLine(results);
            }
        }
    }
    writer.close(true);
    FileUtil::remove(par.db4Index.c_str());

    contigsHeaders.close();
    targetsHeaders.close();
    readerRepToMembers.close();

    return EXIT_SUCCESS;
}
