#include "LocalParameters.h"
#include "PredictionParser.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "MathUtil.h"
#include "itoa.h"
#include "Orf.h"
#include "TranslateNucl.h"

#include <sstream>

#ifdef OPENMP
#include <omp.h>
#endif

void reverseComplement (const std::string & seq, std::string & revCompSeq) {
    size_t seqLength = revCompSeq.size();
    for(size_t i = 0; i < seqLength; ++i) {
        char currNuc = seq[seqLength - i - 1];
        revCompSeq[i] = Orf::complement(currNuc);
    }
}

void preparePredDataAndHeader (const Prediction & pred, const std::string & targetHeaderAcc, const std::string & contigHeaderAcc, 
                                    const char* contigData, std::ostringstream & joinedHeaderStream, std::ostringstream & joinedExonsStream) {
    
    // clear streams:
    joinedHeaderStream.str("");
    joinedExonsStream.str("");

    // initialize header:
    joinedHeaderStream << targetHeaderAcc << "|" << contigHeaderAcc << "|";
    if (pred.strand == PLUS) {
        joinedHeaderStream << "+|";
    } else {
        joinedHeaderStream << "-|";
    }
    joinedHeaderStream << pred.totalBitscore << "|" << pred.combinedEvalue << "|" << pred.numExons << "|";
    joinedHeaderStream << pred.lowContigCoord << "|" << pred.highContigCoord;

    // add all exons:
    int lastTargetPosMatched = -1;
    for (size_t i = 0; i < pred.optimalExonSet.size(); ++i) {
        int targetMatchStart = pred.optimalExonSet[i].targetMatchStart;
        int targetMatchEnd = pred.optimalExonSet[i].targetMatchEnd;
        
        int exonContigStart = pred.optimalExonSet[i].contigStart;
        int exonContigEnd = pred.optimalExonSet[i].contigEnd;
        int exonNucleotideLen = pred.optimalExonSet[i].nucleotideLen;

        if (pred.strand == MINUS) {
            if ((exonContigStart > 0) || (exonContigEnd > 0)) {
                Debug(Debug::ERROR) << "ERROR: strand is MINUS but the contig coordinates are positive. Soemthing is wrong.\n";
                EXIT(EXIT_FAILURE);
            }
        }

        // in order to avoid target overlaps, we remove a few codons from the start of the current exon if needed:
        int exonAdjustedContigStart = exonContigStart;
        int exonAdjustedNucleotideLen = exonNucleotideLen;
        if (lastTargetPosMatched >= targetMatchStart) {
            int diffInAAs = lastTargetPosMatched - targetMatchStart + 1;
            exonAdjustedContigStart += (3 * diffInAAs * pred.strand);
            exonAdjustedNucleotideLen -= (3 * diffInAAs);
        }
        int lowContigCoord = (pred.strand == PLUS) ? exonAdjustedContigStart : (-1 * exonContigEnd);

        // extract the segment from the contig:
        std::string exonContigSeq(&contigData[lowContigCoord], (size_t)exonAdjustedNucleotideLen);
        // update the last AA of the target that was matched:
        lastTargetPosMatched = targetMatchEnd;

        // write the header and data to streams:
        joinedHeaderStream << "|" << abs(exonContigStart) << "[" << abs(exonAdjustedContigStart) << "]:";
        joinedHeaderStream << abs(exonContigEnd) << "[" << abs(exonContigEnd) << "]:";
        joinedHeaderStream << exonNucleotideLen << "[" << exonAdjustedNucleotideLen << "]";

        if (pred.strand == PLUS) {
            joinedExonsStream << exonContigSeq;
        } else {
            std::string exonContigRevCompSeq(exonContigSeq);
            reverseComplement (exonContigSeq, exonContigRevCompSeq);
            joinedExonsStream << exonContigRevCompSeq;
        }
    }
    joinedHeaderStream << "\n";
    joinedExonsStream << "\n";
}

int unitesetstofasta(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, true, 0, 0);

    // db1 = contigsDB (data + header):
    DBReader<unsigned int> contigsData(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    contigsData.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> contigsHeaders(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    contigsHeaders.open(DBReader<unsigned int>::NOSORT);

    // db2 = targetsDB (only the header is used)
    DBReader<unsigned int> targetsHeaders(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    targetsHeaders.open(DBReader<unsigned int>::NOSORT);

    // db3 = predictions per contig
    DBReader<unsigned int> predsPerContig(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    predsPerContig.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_OMIT_FILE);
    writer.open();

    // in case a translated result is required
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));

    Debug::Progress progress(predsPerContig.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        // per thread variables
        std::ostringstream joinedHeaderStream;
        std::ostringstream joinedExonsStream;
        const char *entry[255];

        Prediction plusPred;
        Prediction minusPred;

        char translatedSeqBuff[Parameters::MAX_SEQ_LEN];
        char targetKeyBuff[40];
        
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < predsPerContig.getSize(); id++) {
            progress.updateProgress();

            unsigned int contigKey = predsPerContig.getDbKey(id);

            char *results = predsPerContig.getData(id, thread_idx);
            
            // if the contig has no predictions - move on
            if (*results == '\0') {
                plusPred.clearPred();
                minusPred.clearPred();
                continue;
            }

            unsigned int currTargetKey = 0;
            bool isFirstIteration = true;
            
            // get contig data and header:
            const char* contigData = contigsData.getDataByDBKey(contigKey, thread_idx);
            const char* contigHeader = contigsHeaders.getDataByDBKey(contigKey, thread_idx);
            std::string contigHeaderAcc = Util::parseFastaHeader(contigHeader);

            // process a specific contig
            while (*results != '\0') {
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // each line informs of a prediction and a single exon
                // the first 7 columns describe the entire prediction
                // the last 10 columns describe a single exon
                if (columns != 17) {
                    Debug(Debug::ERROR) << "There should be 17 columns in the input file. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int targetKey = Prediction::getTargetKey(entry);
                int strand = Prediction::getStrand(entry);

                if (isFirstIteration) {
                    currTargetKey = targetKey;
                    isFirstIteration = false;
                }

                // after collecting all the exons for the current target:
                if (targetKey != currTargetKey) {
                    if (targetKey < currTargetKey) {
                        Debug(Debug::ERROR) << "The targets are assumed to be sorted in increasing order. This doesn't seem to be the case.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    
                    const char* targetHeader = targetsHeaders.getDataByDBKey(currTargetKey, thread_idx);
                    std::string targetHeaderAcc;
                    if (par.writeTargetKey == true) {
                        Itoa::u32toa_sse2(static_cast<uint32_t>(currTargetKey), targetKeyBuff);
                        targetHeaderAcc = std::string(targetKeyBuff);
                    } else {
                        targetHeaderAcc = Util::parseFastaHeader(targetHeader);
                    }
                    
                    if (plusPred.optimalExonSet.size() > 0) {
                        preparePredDataAndHeader(plusPred, targetHeaderAcc, contigHeaderAcc, contigData, joinedHeaderStream, joinedExonsStream);
                        std::string result = ">" + joinedHeaderStream.str();
                        writer.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                        result = joinedExonsStream.str();
                        if (par.shouldTranslate == true) {
                            size_t nuclLen = result.size() - 1; // \n at the end of result...
                            if (nuclLen % 3 != 0) {
                                Debug(Debug::ERROR) << "coding sequence does not divide by 3.\n";
                                EXIT(EXIT_FAILURE);
                            }
                            size_t aaLen = nuclLen / 3;
                            translateNucl.translate(translatedSeqBuff, result.c_str(), nuclLen);
                            translatedSeqBuff[aaLen] = '\n';
                            writer.writeData(translatedSeqBuff, (aaLen + 1), 0, thread_idx, false, false);
                        } else {
                            writer.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                        }
                    }
                    if (minusPred.optimalExonSet.size() > 0) {
                        preparePredDataAndHeader(minusPred, targetHeaderAcc, contigHeaderAcc, contigData, joinedHeaderStream, joinedExonsStream);
                        std::string result = ">" + joinedHeaderStream.str();
                        writer.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                        result = joinedExonsStream.str();
                        if (par.shouldTranslate == true) {
                            size_t nuclLen = result.size() - 1; // \n at the end of result...
                            if (nuclLen % 3 != 0) {
                                Debug(Debug::ERROR) << "coding sequence does not divide by 3.\n";
                                EXIT(EXIT_FAILURE);
                            }
                            size_t aaLen = nuclLen / 3;
                            translateNucl.translate(translatedSeqBuff, result.c_str(), nuclLen);
                            translatedSeqBuff[aaLen] = '\n';
                            writer.writeData(translatedSeqBuff, (aaLen + 1), 0, thread_idx, false, false);
                        } else {
                            writer.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                        }
                    }
                    
                    // move to another target:
                    plusPred.clearPred();
                    minusPred.clearPred();
                    currTargetKey = targetKey;
                }
                
                // add exon from same target:
                if (strand == PLUS) {
                    // these fields remain constant between exons of the same prediction
                    plusPred.setByDPRes(entry);
                    plusPred.addExon(entry);
                } else {
                    // these fields remain constant between exons of the same prediction
                    minusPred.setByDPRes(entry);
                    minusPred.addExon(entry);
                }
                results = Util::skipLine(results);
            }

            // handle the last target for the current contig:
            const char* targetHeader = targetsHeaders.getDataByDBKey(currTargetKey, thread_idx);
            std::string targetHeaderAcc;
            if (par.writeTargetKey == true) {
                Itoa::u32toa_sse2(static_cast<uint32_t>(currTargetKey), targetKeyBuff);
                targetHeaderAcc = std::string(targetKeyBuff);
            } else {
                targetHeaderAcc = Util::parseFastaHeader(targetHeader);
            }
            
            if (plusPred.optimalExonSet.size() > 0) {
                preparePredDataAndHeader(plusPred, targetHeaderAcc, contigHeaderAcc, contigData, joinedHeaderStream, joinedExonsStream);
                std::string result = ">" + joinedHeaderStream.str();
                writer.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                result = joinedExonsStream.str();
                if (par.shouldTranslate == true) {
                    size_t nuclLen = result.size() - 1; // \n at the end of result...
                    if (nuclLen % 3 != 0) {
                        Debug(Debug::ERROR) << "coding sequence does not divide by 3.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    size_t aaLen = nuclLen / 3;
                    translateNucl.translate(translatedSeqBuff, result.c_str(), nuclLen);
                    translatedSeqBuff[aaLen] = '\n';
                    writer.writeData(translatedSeqBuff, (aaLen + 1), 0, thread_idx, false, false);
                } else {
                    writer.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                }
            }
            if (minusPred.optimalExonSet.size() > 0) {
                preparePredDataAndHeader(minusPred, targetHeaderAcc, contigHeaderAcc, contigData, joinedHeaderStream, joinedExonsStream);
                std::string result = ">" + joinedHeaderStream.str();
                writer.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                result = joinedExonsStream.str();
                if (par.shouldTranslate == true) {
                    size_t nuclLen = result.size() - 1; // \n at the end of result...
                    if (nuclLen % 3 != 0) {
                        Debug(Debug::ERROR) << "coding sequence does not divide by 3.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    size_t aaLen = nuclLen / 3;
                    translateNucl.translate(translatedSeqBuff, result.c_str(), nuclLen);
                    translatedSeqBuff[aaLen] = '\n';
                    writer.writeData(translatedSeqBuff, (aaLen + 1), 0, thread_idx, false, false);
                } else {
                    writer.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                }
            }
            // move to another contig:
            plusPred.clearPred();
            minusPred.clearPred();
        }
    }
    writer.close(true);
    FileUtil::remove(par.db4Index.c_str());
 
    contigsData.close();
    contigsHeaders.close();
    targetsHeaders.close();
    predsPerContig.close();

    return EXIT_SUCCESS;
}



