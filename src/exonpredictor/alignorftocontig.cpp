#include "LocalParameters.h"
#include "Matcher.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Orf.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int alignorftocontig(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, 3, true, true);

    // contig length is needed for computation:
    DBReader<unsigned int> contigsReader(par.db1.c_str(), par.db1Index.c_str());
    contigsReader.open(DBReader<unsigned int>::NOSORT);

    // info will be obtained from orf headers:
    DBReader<unsigned int> orfHeadersReader(par.hdr2.c_str(), par.hdr2Index.c_str());
    orfHeadersReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // writing in alignment format:
    DBWriter alignmentFormatWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    alignmentFormatWriter.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif 
        char orfToContigBuffer[LINE_MAX];
        
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < orfHeadersReader.getSize(); ++id) {
            Debug::printProgress(id);

            unsigned int orfKey = orfHeadersReader.getDbKey(id);
            char * orfHeader = orfHeadersReader.getData(id);
            Orf::SequenceLocation orfLocOnContigParsed;
            orfLocOnContigParsed = Orf::parseOrfHeader(orfHeader);

            // get contig key and its length in nucleotides
            int contigKey = orfLocOnContigParsed.id;
            unsigned int contigId = contigsReader.getId(contigKey);

            size_t contigLenWithEndings = contigsReader.getSeqLens(contigId);
            if (contigLenWithEndings < 2) {
                Debug(Debug::ERROR) << "Invalid contig record has less than two bytes\n";
                EXIT(EXIT_FAILURE);
            }

            size_t contigLen = contigLenWithEndings - 2; // remove \n\0

            // orf search returns the position right after the orf, this keeps consitency with alignemnt format
            // if strand == -1 (reverse), we need to recompute the coordinates with respect to the positive strand and swap them
            size_t contigFromWithStrand = (orfLocOnContigParsed.strand > 0) ? orfLocOnContigParsed.from : (contigLen - orfLocOnContigParsed.from - 1);
            size_t contigToWithStrand = (orfLocOnContigParsed.strand > 0) ? (orfLocOnContigParsed.to - 1) : (contigLen - orfLocOnContigParsed.to);

            // compute orf length
            int orfLen = (orfLocOnContigParsed.strand > 0) ? (contigToWithStrand - contigFromWithStrand + 1) : (contigFromWithStrand - contigToWithStrand + 1);

            Matcher::result_t orfToContigResult(contigKey, 1, 1, 0, 1, 0, orfLen, 0, (orfLen - 1), orfLen, contigFromWithStrand, contigToWithStrand, contigLen, "");
            size_t len = Matcher::resultToBuffer(orfToContigBuffer, orfToContigResult, true);
            alignmentFormatWriter.writeData(orfToContigBuffer, len, orfKey, thread_idx);
        }
    }

    // cleanup
    alignmentFormatWriter.close();
    orfHeadersReader.close();
    contigsReader.close();
    
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}




