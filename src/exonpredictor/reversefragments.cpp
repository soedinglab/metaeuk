#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>

#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

int reversefragments(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, 2, true, true);

    DBReader<unsigned int> aaFragmentsReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    aaFragmentsReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::string aaFragmentsReversed = par.db2;
    std::string aaFragmentsReversedIndex = par.db2Index;
    DBWriter aaFragmentsReversedWriter(aaFragmentsReversed.c_str(), aaFragmentsReversedIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aaFragmentsReversedWriter.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char revAaFragmentBuffer[32000];

#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < aaFragmentsReader.getSize(); id++) {
            Debug::printProgress(id);

            unsigned int fragmentKey = aaFragmentsReader.getDbKey(id);
            char *aaFragment = aaFragmentsReader.getData(id, thread_idx);
            int lenFragment = strlen(aaFragment) - 1; // last char is '\n'

            if ((lenFragment - 1) > 32000) {
                Debug(Debug::ERROR) << "AA fragment is longer than what is handled...\n";
                EXIT(EXIT_FAILURE);
            }

            for (int i = 0; i < ((lenFragment / 2) + 1); ++i) {
                revAaFragmentBuffer[i] = aaFragment[lenFragment - i - 1];
                revAaFragmentBuffer[lenFragment - i - 1] = aaFragment[i];
            }
            revAaFragmentBuffer[lenFragment] = '\n';
            
            aaFragmentsReversedWriter.writeData(revAaFragmentBuffer, (lenFragment + 1), fragmentKey, thread_idx, true);
        }
    }

    // cleanup
    aaFragmentsReversedWriter.close();
    aaFragmentsReader.close();

    FileUtil::symlinkAbs(par.hdr1, par.hdr2);
    FileUtil::symlinkAbs(par.hdr1Index, par.hdr2Index);
    
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}



