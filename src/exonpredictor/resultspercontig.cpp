#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Orf.h"
#include "AlignmentSymmetry.h"
#include "Timer.h"
#include "IndexReader.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

struct compareByTarget {
    bool operator() (const std::pair<Matcher::result_t, Matcher::result_t>& lhs, const std::pair<Matcher::result_t, Matcher::result_t>& rhs) const {
        // sort by target id
        if (lhs.first.dbKey < rhs.first.dbKey) {
            return true;
        }
        // if a contig hits the same target with two orfs - sort by orf id
        if (lhs.second.dbKey < rhs.second.dbKey) {
            return true;
        }
        return false;
    }
};

int resultspercontig(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool touch = par.preloadMode != Parameters::PRELOAD_MODE_MMAP;
    IndexReader qContigDbr(par.db1.c_str(), par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX) : 0, DBReader<unsigned int>::USE_INDEX);
    IndexReader qOrfDbr(par.db2.c_str(), par.threads, IndexReader::HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    // contig length is needed for computation:
    DBReader<unsigned int> contigsReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    contigsReader.open(DBReader<unsigned int>::NOSORT);

    // info will be obtained from orf headers:
    DBReader<unsigned int> orfHeadersReader(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    orfHeadersReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    //IndexReader tSourceDbr(par.db3.c_str(), par.threads, IndexReader::HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

#ifdef OPENMP
    unsigned int totalThreads = par.threads;
#else
    unsigned int totalThreads = 1;
#endif

    unsigned int localThreads = totalThreads;
    if (alnDbr.getSize() <= totalThreads) {
        localThreads = alnDbr.getSize();
    }

    // Compute mapping from contig -> orf[] from orf[]->contig in headers
    unsigned int *contigLookup = NULL;
    unsigned int *contigOffsets = NULL;
    char *contigExists = NULL;
    unsigned int maxContigKey = 0;
    Timer timer;
    Debug(Debug::INFO) << "Computing ORF lookup\n";
    unsigned int maxOrfKey = alnDbr.getLastKey();
    unsigned int *orfLookup = new unsigned int[maxOrfKey + 2]();
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i <= maxOrfKey; ++i) {
            size_t queryId = qOrfDbr.sequenceReader->getId(i);
            if (queryId == UINT_MAX) {
                orfLookup[i] = UINT_MAX;
                continue;
            }
            unsigned int queryKey = qOrfDbr.sequenceReader->getDbKey(queryId);
            char *header = qOrfDbr.sequenceReader->getData(queryId, thread_idx);
            Orf::SequenceLocation qloc = Orf::parseOrfHeader(header);
            unsigned int id = (qloc.id != UINT_MAX) ? qloc.id : queryKey;
            orfLookup[i] = id;
        }
    }
    Debug(Debug::INFO) << "Computing contig offsets\n";
    maxContigKey = qContigDbr.sequenceReader->getLastKey();
    unsigned int *contigSizes = new unsigned int[maxContigKey + 2]();
#pragma omp parallel for schedule(static) num_threads(localThreads)
    for (size_t i = 0; i <= maxOrfKey ; ++i) {
        if(orfLookup[i] == UINT_MAX){
            continue;
        }
        __sync_fetch_and_add(&(contigSizes[orfLookup[i]]), 1);
    }
    contigOffsets = contigSizes;

    AlignmentSymmetry::computeOffsetFromCounts(contigOffsets, maxContigKey + 1);

    contigExists = new char[maxContigKey + 1]();
#pragma omp parallel for schedule(static) num_threads(localThreads)
    for (size_t i = 0; i < qContigDbr.sequenceReader->getSize(); ++i) {
        contigExists[qContigDbr.sequenceReader->getDbKey(i)] = 1;
    }

    Debug(Debug::INFO) << "Computing contig lookup\n";
    contigLookup = new unsigned int[maxOrfKey + 2]();
#pragma omp parallel for schedule(static) num_threads(localThreads)
    for (size_t i = 0; i <= maxOrfKey; ++i) {
        if(orfLookup[i] == UINT_MAX){
            continue;
        }
        size_t offset = __sync_fetch_and_add(&(contigOffsets[orfLookup[i]]), 1);
        contigLookup[offset] = i;
    }
    delete[] orfLookup;

    for (unsigned int i = maxContigKey + 1; i > 0; --i) {
        contigOffsets[i] = contigOffsets[i - 1];
    }
    contigOffsets[0] = 0;
    Debug(Debug::INFO) << "Time for contig lookup: " << timer.lap() << "\n";

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

    size_t entryCount = maxContigKey + 1;
    Debug::Progress progress(entryCount);
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char * buffer = new char[65536];

        std::string ss;
        ss.reserve(1024);

        std::vector<std::pair<Matcher::result_t, Matcher::result_t>> results;
        results.reserve(300);

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < entryCount; ++i) {
            progress.updateProgress();

            if (contigExists[i] == 0) {
                continue;
            }

            unsigned int contigKey = i;
            size_t contigId = qContigDbr.sequenceReader->getId(contigKey);

            //unsigned int qLen = std::max(qContigDbr.sequenceReader->getSeqLens(contigId), static_cast<size_t>(2)) - 2;
            unsigned int *orfKeys = &contigLookup[contigOffsets[i]];
            size_t orfCount = contigOffsets[i + 1] - contigOffsets[i];
            for (unsigned int j = 0; j < orfCount; ++j) {
                unsigned int orfKey = orfKeys[j];
                size_t orfId = alnDbr.getId(orfKey);
                // this is needed when alnDbr does not contain all identifiers of the queryDB
                if (orfId == UINT_MAX) {
                    continue;
                }
                
                Matcher::result_t orfToContig = Orf::getFromDatabase(orfKey, contigsReader, orfHeadersReader, thread_idx);
                // hack orfToContig to retain the orf key and not the contig key (that we will have as the db key)
                orfToContig.dbKey = orfKey;

                char *data = alnDbr.getData(orfId, thread_idx);
                while (*data != '\0') {
                    results.emplace_back(std::make_pair(Matcher::parseAlignmentRecord(data, true), orfToContig));
                    data = Util::skipLine(data);
                }
            }
        
            std::stable_sort(results.begin(), results.end(), compareByTarget());

            for (size_t i = 0; i < results.size(); i++) {
                Matcher::result_t &orfToTarget = results[i].first;
                bool hasBacktrace = (orfToTarget.backtrace.size() > 0);
                size_t len = Matcher::resultToBuffer(buffer, orfToTarget, hasBacktrace, false);
                // -1 for newline
                ss.append(buffer, len - 1);
                // add "\t"
                ss.append("\t");
                Matcher::result_t &orfToContig = results[i].second;
                len = Matcher::resultToBuffer(buffer, orfToContig, false, false);
                ss.append(buffer, len);
            }
            resultWriter.writeData(ss.c_str(), ss.length(), contigKey, thread_idx);

            ss.clear();
            results.clear();
        }
        delete[] buffer;
    }
    resultWriter.close();

    if (contigLookup != NULL) {
        delete[] contigLookup;
    }

    if (contigOffsets != NULL) {
        delete[] contigOffsets;
    }

    if (contigExists != NULL) {
        delete[] contigExists;
    }
    alnDbr.close();

    contigsReader.close();
    orfHeadersReader.close();

    return EXIT_SUCCESS;
}
