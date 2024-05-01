#include "LocalParameters.h"
#include "PredictionParser.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
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

int findStartInString (const std::string & seq) {
    int lastPosOfClosestStart = 0;
    // In a case like this: ATGcccATGccc
    // The returned index will be 8 (the last char of the last ATG in the string)
    for (size_t i = 0; i < seq.size(); i += 3) {
        if ((seq.substr(i,3).compare("ATG") == 0) or (seq.substr(i,3).compare("atg") == 0)) {
            lastPosOfClosestStart = i + 2;
        }
    }
    return(lastPosOfClosestStart);
}

int searchForStartBeforeFirstExon(Prediction & pred, const char* contigData, std::ostringstream & joinedExonsStream, const int searchLen) {
    // check if first exon's first codon is start and if so return 0 - already a start codon
    if (pred.strand == PLUS) {
        std::string firstCodon(&contigData[pred.lowContigCoord], 3);
        int lastPosOfStart = findStartInString(firstCodon);
        if (lastPosOfStart > 0) {
            return 0;
        }
    }
    if (pred.strand == MINUS) {
        std::string firstCodon(&contigData[pred.highContigCoord - 2], 3);
        std::string firstCodonRevCompSeq(firstCodon);
        reverseComplement (firstCodon, firstCodonRevCompSeq);
        int lastPosOfStart = findStartInString(firstCodonRevCompSeq);
        if (lastPosOfStart > 0) {
            return 0;
        }
    }

    int seachLenLegal = searchLen - (searchLen % 3);

    // case PLUS                                        
    int coordToStartSearch = pred.lowContigCoord - seachLenLegal;
    int posAfterStopCodon = pred.optimalExonSet[0].potentialExonContigStartBeforeTrim;
    // case MINUS
    if (pred.strand == MINUS) {
        coordToStartSearch = pred.highContigCoord + 1;
        posAfterStopCodon = pred.optimalExonSet[0].potentialExonContigEndBeforeTrim;
    }

    // be careful at the edge of the first exon - don't go over stop
    // case PLUS
    if ((pred.strand == PLUS) && (coordToStartSearch < posAfterStopCodon)) {
        coordToStartSearch = posAfterStopCodon;
        seachLenLegal = pred.lowContigCoord - coordToStartSearch;
    }
    // case MINUS
    if ((pred.strand == MINUS) && ((posAfterStopCodon - pred.highContigCoord) < (size_t)seachLenLegal)) {
        seachLenLegal = (posAfterStopCodon - pred.highContigCoord);
        if (seachLenLegal % 3 != 0) {
            Debug(Debug::ERROR) << "ERROR: seachLenLegal mod 3 is not 0.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // extract the segment from the contig:
    std::string beforeFirstExonSeq(&contigData[coordToStartSearch], (size_t)seachLenLegal);
    std::string beforeFirstExonSeqRevCompSeq(beforeFirstExonSeq);
    if (pred.strand == MINUS) {
        reverseComplement (beforeFirstExonSeq, beforeFirstExonSeqRevCompSeq);
    }

    size_t posLastCharStart = 0;
    // if not found, findStartInString will return 0, not changing anything
    if (pred.strand == PLUS) {
        posLastCharStart = findStartInString(beforeFirstExonSeq);
        if (posLastCharStart > 0) {
            std::string seqToAdd(beforeFirstExonSeq.substr(posLastCharStart - 2));
            int numNucsAdded = seqToAdd.size();
            pred.lowContigCoord = pred.lowContigCoord - numNucsAdded;
            joinedExonsStream << seqToAdd;
            return numNucsAdded;
        }
    } else {
        posLastCharStart = findStartInString(beforeFirstExonSeqRevCompSeq);
        if (posLastCharStart > 0) {
            std::string seqToAdd(beforeFirstExonSeqRevCompSeq.substr(posLastCharStart - 2));
            int numNucsAdded = seqToAdd.size();
            pred.highContigCoord = pred.highContigCoord + numNucsAdded;
            joinedExonsStream << seqToAdd;
            return numNucsAdded;
        }
    }

    // start codon not found - not changing anything, returning 0
    return 0;
}

void preparePredDataAndHeader (Prediction & pred, const std::string & targetHeaderAcc, const std::string & contigHeaderAcc, 
                                    const char* contigData, std::ostringstream & joinedHeaderStream, std::ostringstream & joinedExonsStream, 
                                    const int writeFragCoords, const int lenSearchStart, const size_t contigLen) {
    
    // clear streams:
    joinedHeaderStream.str("");
    joinedExonsStream.str("");

    // if list is empty - return:
    size_t numExonsInPred = pred.optimalExonSet.size();
    if (numExonsInPred == 0) {
        return;
    }

    // initialize header:
    joinedHeaderStream << targetHeaderAcc << "|" << contigHeaderAcc << "|";
    if (pred.strand == PLUS) {
        joinedHeaderStream << "+|";
    } else {
        joinedHeaderStream << "-|";
    }
    joinedHeaderStream << pred.totalBitscore << "|" << pred.combinedEvalue << "|" << pred.numExons << "|";
    if (lenSearchStart == 0) {
        // default case: user doesn't want to search for start codon before the first exon:
        joinedHeaderStream << pred.lowContigCoord << "|" << pred.highContigCoord;
    } else {
        // special case: user wants to search for start codon before the first exon:
        // if found - padd it up/downstream and indicate the padding != 0 in []
        int numNucsAdded = searchForStartBeforeFirstExon(pred, contigData, joinedExonsStream, lenSearchStart);
        if (pred.strand == PLUS) {
            joinedHeaderStream << pred.lowContigCoord << "[" << numNucsAdded << "]" << "|" << pred.highContigCoord;
        } else {
            joinedHeaderStream << pred.lowContigCoord << "|" << pred.highContigCoord << "[" << numNucsAdded << "]";
        }
    }
    

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
            exonAdjustedContigStart += (3 * diffInAAs);
            exonAdjustedNucleotideLen -= (3 * diffInAAs);
        }
        int exonLowContigCoord = (pred.strand == PLUS) ? exonAdjustedContigStart : (-1 * exonContigEnd);
        
        // update the field:
        pred.optimalExonSet[i].adjustedContigStart = exonAdjustedContigStart;
        
        // extract the segment from the contig:
        std::string exonContigSeq(&contigData[exonLowContigCoord], (size_t)exonAdjustedNucleotideLen);
        // update the last AA of the target that was matched:
        lastTargetPosMatched = targetMatchEnd;

        // write the header and data to streams:
        joinedHeaderStream << "|";
        if (writeFragCoords == true) {
            joinedHeaderStream << "[" << pred.optimalExonSet[i].potentialExonContigStartBeforeTrim << "]";
        }
        joinedHeaderStream << abs(exonContigStart) << "[" << abs(exonAdjustedContigStart) << "]:";
        if (writeFragCoords == true) {
            joinedHeaderStream << "[" << pred.optimalExonSet[i].potentialExonContigEndBeforeTrim << "]";
        }
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

    // if flag is on, add the stop codon after the last exon (if exists)
    if ((writeFragCoords == true) && 
        (pred.optimalExonSet[numExonsInPred - 1].potentialExonContigEndBeforeTrim == abs(pred.optimalExonSet[numExonsInPred - 1].contigEnd))) {
        
        int lastCodingPosition = pred.optimalExonSet[numExonsInPred - 1].potentialExonContigEndBeforeTrim;
        int strand = pred.optimalExonSet[numExonsInPred - 1].strand;
        int stopCodonPosition = 0;
        if (strand == PLUS) {
            stopCodonPosition = lastCodingPosition + 1;
        } else {
            stopCodonPosition = lastCodingPosition - 3;
        }
        
        // handle edge case of last codon on the edge of the contig. 
        // don't touch memory that shouldn't be touched.
        if ((stopCodonPosition <= (int)(contigLen - 2)) && (stopCodonPosition >= 0)) {
            std::string stopCodonSeq(&contigData[stopCodonPosition], 3);
            if (strand == PLUS) {
                joinedExonsStream << stopCodonSeq;
            } else {
                std::string stopCodonSeqRevCompSeq(stopCodonSeq);
                reverseComplement (stopCodonSeq, stopCodonSeqRevCompSeq);
                joinedExonsStream << stopCodonSeqRevCompSeq;
            }
        }
    }

    joinedHeaderStream << "\n";
    joinedExonsStream << "\n";
}

void preparePredHeaderToInfo (const unsigned int contigKey, const Prediction & pred, const std::string & joinedHeaderStr, 
                                std::ostringstream & joinedPredHeadToInfoStream) {
    // clear stream:
    joinedPredHeadToInfoStream.str("");

    // structure mimics the headers produced by extractorfs (Orf::writeOrfHeader)
    // the first columns are therefore:
    // contigkey, contigStartPosition+contigLenIncludingIntrons, 0
    // followed by MetaEuk columns:
    // targetKey, strand, predHeader

    size_t contigLenIncludingIntrons = pred.highContigCoord - pred.lowContigCoord + 1;

    joinedPredHeadToInfoStream << contigKey << "\t";
    if (pred.strand == PLUS) {
        joinedPredHeadToInfoStream << pred.lowContigCoord << "+" << contigLenIncludingIntrons << "\t";
    } else {
        joinedPredHeadToInfoStream << pred.highContigCoord << "-" << contigLenIncludingIntrons << "\t";
    }
    joinedPredHeadToInfoStream << "0\t" << pred.targetKey << "\t" << pred.strand << "\t";
    // no need for \n because the header already has one!
    joinedPredHeadToInfoStream << joinedHeaderStr;
}

void preparePredHeaderToGFF (const std::string & contigHeaderAcc, const Prediction & pred, const std::string & targetHeaderAcc, 
                                std::ostringstream & joinedPredHeadToGFFStream) {
    // clear stream:
    joinedPredHeadToGFFStream.str("");

    // structure follows GFF format
    // Enteties "gene" and "mRNA" are the same
    // "exon" reflects the first reported coords
    // "CDS" reflects the second reported coords, might be trimmed a bit due to
    // the allowed short overlap on the target
    // GFF: 1-base offset so we add 1 to MetaEuk coords
    std::string sourceStr = "MetaEuk";
    std::string strandStr = "+";
    if (pred.strand != PLUS) {
        strandStr[strandStr.length() - 1] = '-';
    }

    // the combination of T+C+S and lowContigStart should be unique
    std::string TCS = targetHeaderAcc + "|" + contigHeaderAcc + "|" + strandStr + "|" + std::to_string(pred.lowContigCoord);
    
    // start with writing the gene and mrna:
    const char *levels[4] = { "gene", "mRNA", "exon", "CDS" };  
    for (size_t i = 0; i < 2; i++) {
        joinedPredHeadToGFFStream << contigHeaderAcc << "\t" << sourceStr << "\t" ;
        joinedPredHeadToGFFStream << levels[i] << "\t";
        joinedPredHeadToGFFStream << (pred.lowContigCoord + 1) << "\t" << (pred.highContigCoord + 1);
        joinedPredHeadToGFFStream << "\t" << pred.totalBitscore << "\t";
        joinedPredHeadToGFFStream << strandStr << "\t";
        joinedPredHeadToGFFStream << "." << "\t";
        joinedPredHeadToGFFStream << "Target_ID=" << targetHeaderAcc << ";";
        if (strcmp(levels[i],"gene") == 0) {
            joinedPredHeadToGFFStream << "TCS_ID=" << TCS << "\n";
        } else {
            joinedPredHeadToGFFStream << "TCS_ID=" << TCS << "_mRNA;";
            joinedPredHeadToGFFStream << "Parent=" << TCS << "\n";
        }
    }

    // write exons and cds
    for (size_t j = 0; j < pred.numExons; ++j) {
        
        size_t exonContigStart = abs(pred.optimalExonSet[j].contigStart) + 1;
        size_t exonAdjustedContigStart = abs(pred.optimalExonSet[j].adjustedContigStart) + 1;
        size_t exonContigEnd = abs(pred.optimalExonSet[j].contigEnd) + 1;
        
        for (size_t i = 2; i < 4; i++) {
            joinedPredHeadToGFFStream << contigHeaderAcc << "\t" << sourceStr << "\t" ;
            joinedPredHeadToGFFStream << levels[i] << "\t";

            if (strcmp(levels[i],"exon") == 0) {
                if (pred.strand == PLUS) {
                    joinedPredHeadToGFFStream << exonContigStart << "\t" << exonContigEnd;
                } else {
                    joinedPredHeadToGFFStream << exonContigEnd << "\t" << exonContigStart;
                }
                joinedPredHeadToGFFStream << "\t" << pred.optimalExonSet[j].bitScore << "\t";
                joinedPredHeadToGFFStream << strandStr << "\t";
                joinedPredHeadToGFFStream << "." << "\t";
                joinedPredHeadToGFFStream << "Target_ID=" << targetHeaderAcc << ";";
                joinedPredHeadToGFFStream << "TCS_ID=" << TCS << "_exon_" << j << ";";
                joinedPredHeadToGFFStream << "Parent=" << TCS << "_mRNA\n";
            } else {
                if (pred.strand == PLUS) {
                    joinedPredHeadToGFFStream << exonAdjustedContigStart << "\t" << exonContigEnd;
                } else {
                    joinedPredHeadToGFFStream << exonContigEnd << "\t" << exonAdjustedContigStart;
                }
                joinedPredHeadToGFFStream << "\t" << pred.optimalExonSet[j].bitScore << "\t";
                joinedPredHeadToGFFStream << strandStr << "\t";
                joinedPredHeadToGFFStream << "." << "\t";
                joinedPredHeadToGFFStream << "Target_ID=" << targetHeaderAcc << ";";
                joinedPredHeadToGFFStream << "TCS_ID=" << TCS << "_CDS_" << j << ";";
                joinedPredHeadToGFFStream << "Parent=" << TCS << "_exon_" << j << "\n";
            }
        }
    }
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
    
    std::string fastaAaFileName = par.db4 + ".fas";
    std::string fastaAaFileNameIndex = par.db4Index;

    // out AA fasta
    DBWriter fastaAaWriter(fastaAaFileName.c_str(), fastaAaFileNameIndex.c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    fastaAaWriter.open();

    std::string fastaCodonFileName = par.db4 + ".codon.fas";
    std::string fastaCodonFileNameIndex = par.db4 + ".codon.index";

    // out codon fasta
    DBWriter fastaCodonWriter(fastaCodonFileName.c_str(), fastaCodonFileNameIndex.c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    fastaCodonWriter.open();

    // out mapping - MetaEuk header to contig, target, etc. Mimicking the headers produced by extractorfs so this can later be plugged in easily
    std::string mapFileName = par.db4 + ".headersMap.tsv";
    std::string mapFileNameIndex = par.db4 + ".headersMap.tsv.index"; // not used
    DBWriter mapWriter(mapFileName.c_str(), mapFileNameIndex.c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    mapWriter.open();

    // out GFF
    std::string gffFileName = par.db4 + ".gff";
    std::string gffFileNameIndex = par.db4 + ".gff.index"; // not used
    DBWriter gffWriter(gffFileName.c_str(), gffFileNameIndex.c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    gffWriter.open();

    // for the translated result
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
    if ((par.translationTable != 1) && (par.lenSearchStart > 0)) {
        Debug(Debug::WARNING) << "Selected translation table is not canonical and search for start is turned on. Please note that only ATG/atg is considered a start codon for the search!\n";
    }
    // for now, opting to use only ATG/atg as start codon because:
    // In Euks ATG is super dominant, unlike proks
    // The other two alternative codons in the standard table: TTG, CTG are NOT translated as M inside the gene
    // std::vector<std::string> startCodonsInTable = translateNucl.getStartCodons();
    // for (std::vector<std::string>::iterator t=startCodonsInTable.begin(); t!=startCodonsInTable.end(); ++t) {
        // std::cout << *t << std::endl;
    // }

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
        std::ostringstream predHeaderToInfoStream;
        std::ostringstream predHeaderToGffStream;
        const char *entry[255];

        Prediction plusPred;
        Prediction minusPred;

        size_t translatedSeqBuffSize = par.maxSeqLen * sizeof(char);
        char* translatedSeqBuff = (char*)malloc(translatedSeqBuffSize);

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
            unsigned int currLowContigCoord = 0;
            
            // get contig data and header:
            size_t contigId = contigsData.getId(contigKey);
            if (contigId == UINT_MAX) {
                Debug(Debug::ERROR) << "Sequence " << contigKey << " does not exist in the sequence database\n";
                EXIT(EXIT_FAILURE);
            }
            const char* contigData = contigsData.getData(contigId, thread_idx);
            size_t contigLen = contigsData.getSeqLen(contigId);
            const char* contigHeader = contigsHeaders.getDataByDBKey(contigKey, thread_idx);
            std::string contigHeaderAcc = Util::parseFastaHeader(contigHeader);

            // process a specific contig
            while (*results != '\0') {
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // each line informs of a prediction and a single exon
                // the first 7 columns describe the entire prediction
                // the last 12 columns describe a single exon
                if (columns != 19) {
                    Debug(Debug::ERROR) << "There should be 19 columns in the input file. This doesn't seem to be the case.\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int targetKey = Prediction::getTargetKey(entry);
                int strand = Prediction::getStrand(entry);
                unsigned int lowContigCoord = Prediction::getLowContigCoord(entry);
                
                if (isFirstIteration) {
                    currTargetKey = targetKey;
                    isFirstIteration = false;
                    currLowContigCoord = lowContigCoord;
                }

                // after collecting all the exons for the current prediction:
                if ((targetKey != currTargetKey) || (lowContigCoord != currLowContigCoord)) {
                    if (targetKey < currTargetKey) {
                        Debug(Debug::ERROR) << "The targets are assumed to be sorted in increasing order. This doesn't seem to be the case.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    
                    const char* targetHeader = targetsHeaders.getDataByDBKey(currTargetKey, thread_idx);
                    std::string targetHeaderAcc;
                    if (par.writeTargetKey == true) {
                        targetHeaderAcc = SSTR(currTargetKey);
                    } else {
                        targetHeaderAcc = Util::parseFastaHeader(targetHeader);
                    }
                    
                    if (plusPred.optimalExonSet.size() > 0) {
                        preparePredDataAndHeader(plusPred, targetHeaderAcc, contigHeaderAcc, contigData, joinedHeaderStream, joinedExonsStream, par.writeFragCoords, par.lenSearchStart, contigLen);
                        std::string result = ">" + joinedHeaderStream.str();
                        fastaAaWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                        fastaCodonWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                        
                        preparePredHeaderToInfo(contigKey, plusPred, joinedHeaderStream.str(), predHeaderToInfoStream);
                        std::string headerInfo = predHeaderToInfoStream.str();
                        mapWriter.writeData(headerInfo.c_str(), headerInfo.size(), 0, thread_idx, false, false);

                        // write GFF:
                        preparePredHeaderToGFF(contigHeaderAcc, plusPred, targetHeaderAcc, predHeaderToGffStream);
                        std::string gffInfo = predHeaderToGffStream.str();
                        gffWriter.writeData(gffInfo.c_str(), gffInfo.size(), 0, thread_idx, false, false);

                        result = joinedExonsStream.str();
                        size_t nuclLen = result.size() - 1; // \n at the end of result...
                        if (nuclLen % 3 != 0) {
                            Debug(Debug::ERROR) << "coding sequence does not divide by 3.\n";
                            EXIT(EXIT_FAILURE);
                        }
                        size_t aaLen = nuclLen / 3;
                        if ((aaLen + 1) > translatedSeqBuffSize) {
                            translatedSeqBuffSize = (aaLen + 1) * 1.5 * sizeof(char);
                            translatedSeqBuff = (char*)realloc(translatedSeqBuff, translatedSeqBuffSize);
                            Util::checkAllocation(translatedSeqBuff, "Cannot reallocate translatedSeqBuff");
                        }
                        translateNucl.translate(translatedSeqBuff, result.c_str(), nuclLen);
                        translatedSeqBuff[aaLen] = '\n';
                        fastaAaWriter.writeData(translatedSeqBuff, (aaLen + 1), 0, thread_idx, false, false);
                        fastaCodonWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                    }
                    if (minusPred.optimalExonSet.size() > 0) {
                        preparePredDataAndHeader(minusPred, targetHeaderAcc, contigHeaderAcc, contigData, joinedHeaderStream, joinedExonsStream, par.writeFragCoords, par.lenSearchStart, contigLen);
                        std::string result = ">" + joinedHeaderStream.str();
                        fastaAaWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                        fastaCodonWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);

                        preparePredHeaderToInfo(contigKey, minusPred, joinedHeaderStream.str(), predHeaderToInfoStream);
                        std::string headerInfo = predHeaderToInfoStream.str();
                        mapWriter.writeData(headerInfo.c_str(), headerInfo.size(), 0, thread_idx, false, false);

                        // write GFF:
                        preparePredHeaderToGFF(contigHeaderAcc, minusPred, targetHeaderAcc, predHeaderToGffStream);
                        std::string gffInfo = predHeaderToGffStream.str();
                        gffWriter.writeData(gffInfo.c_str(), gffInfo.size(), 0, thread_idx, false, false);

                        result = joinedExonsStream.str();
                        size_t nuclLen = result.size() - 1; // \n at the end of result...
                        if (nuclLen % 3 != 0) {
                            Debug(Debug::ERROR) << "coding sequence does not divide by 3.\n";
                            EXIT(EXIT_FAILURE);
                        }
                        size_t aaLen = nuclLen / 3;
                        if ((aaLen + 1) > translatedSeqBuffSize) {
                            translatedSeqBuffSize = (aaLen + 1) * 1.5 * sizeof(char);
                            translatedSeqBuff = (char*)realloc(translatedSeqBuff, translatedSeqBuffSize);
                            Util::checkAllocation(translatedSeqBuff, "Cannot reallocate translatedSeqBuff");
                        }
                        translateNucl.translate(translatedSeqBuff, result.c_str(), nuclLen);
                        translatedSeqBuff[aaLen] = '\n';
                        fastaAaWriter.writeData(translatedSeqBuff, (aaLen + 1), 0, thread_idx, false, false);
                        fastaCodonWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                    }
                    
                    // move to another prediction:
                    plusPred.clearPred();
                    minusPred.clearPred();
                    currTargetKey = targetKey;
                    currLowContigCoord = lowContigCoord;
                }
                
                // add exon from same prediction:
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

            // handle the last prediction for the current contig:
            const char* targetHeader = targetsHeaders.getDataByDBKey(currTargetKey, thread_idx);
            std::string targetHeaderAcc;
            if (par.writeTargetKey == true) {
                targetHeaderAcc = SSTR(currTargetKey);
            } else {
                targetHeaderAcc = Util::parseFastaHeader(targetHeader);
            }
            
            if (plusPred.optimalExonSet.size() > 0) {
                preparePredDataAndHeader(plusPred, targetHeaderAcc, contigHeaderAcc, contigData, joinedHeaderStream, joinedExonsStream, par.writeFragCoords, par.lenSearchStart, contigLen);
                std::string result = ">" + joinedHeaderStream.str();
                fastaAaWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                fastaCodonWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                
                preparePredHeaderToInfo(contigKey, plusPred, joinedHeaderStream.str(), predHeaderToInfoStream);
                std::string headerInfo = predHeaderToInfoStream.str();
                mapWriter.writeData(headerInfo.c_str(), headerInfo.size(), 0, thread_idx, false, false);

                // write GFF:
                preparePredHeaderToGFF(contigHeaderAcc, plusPred, targetHeaderAcc, predHeaderToGffStream);
                std::string gffInfo = predHeaderToGffStream.str();
                gffWriter.writeData(gffInfo.c_str(), gffInfo.size(), 0, thread_idx, false, false);
                
                result = joinedExonsStream.str();
                size_t nuclLen = result.size() - 1; // \n at the end of result...
                if (nuclLen % 3 != 0) {
                    Debug(Debug::ERROR) << "coding sequence does not divide by 3.\n";
                    EXIT(EXIT_FAILURE);
                }
                size_t aaLen = nuclLen / 3;
                if ((aaLen + 1) > translatedSeqBuffSize) {
                    translatedSeqBuffSize = (aaLen + 1) * 1.5 * sizeof(char);
                    translatedSeqBuff = (char*)realloc(translatedSeqBuff, translatedSeqBuffSize);
                    Util::checkAllocation(translatedSeqBuff, "Cannot reallocate translatedSeqBuff");
                }
                translateNucl.translate(translatedSeqBuff, result.c_str(), nuclLen);
                translatedSeqBuff[aaLen] = '\n';
                fastaAaWriter.writeData(translatedSeqBuff, (aaLen + 1), 0, thread_idx, false, false);
                fastaCodonWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
            }
            if (minusPred.optimalExonSet.size() > 0) {
                preparePredDataAndHeader(minusPred, targetHeaderAcc, contigHeaderAcc, contigData, joinedHeaderStream, joinedExonsStream, par.writeFragCoords, par.lenSearchStart, contigLen);
                std::string result = ">" + joinedHeaderStream.str();
                fastaAaWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                fastaCodonWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
                
                preparePredHeaderToInfo(contigKey, minusPred, joinedHeaderStream.str(), predHeaderToInfoStream);
                std::string headerInfo = predHeaderToInfoStream.str();
                mapWriter.writeData(headerInfo.c_str(), headerInfo.size(), 0, thread_idx, false, false);
                
                // write GFF:
                preparePredHeaderToGFF(contigHeaderAcc, minusPred, targetHeaderAcc, predHeaderToGffStream);
                std::string gffInfo = predHeaderToGffStream.str();
                gffWriter.writeData(gffInfo.c_str(), gffInfo.size(), 0, thread_idx, false, false);

                result = joinedExonsStream.str();
                size_t nuclLen = result.size() - 1; // \n at the end of result...
                if (nuclLen % 3 != 0) {
                    Debug(Debug::ERROR) << "coding sequence does not divide by 3.\n";
                    EXIT(EXIT_FAILURE);
                }
                size_t aaLen = nuclLen / 3;
                if ((aaLen + 1) > translatedSeqBuffSize) {
                    translatedSeqBuffSize = (aaLen + 1) * 1.5 * sizeof(char);
                    translatedSeqBuff = (char*)realloc(translatedSeqBuff, translatedSeqBuffSize);
                    Util::checkAllocation(translatedSeqBuff, "Cannot reallocate translatedSeqBuff");
                }
                translateNucl.translate(translatedSeqBuff, result.c_str(), nuclLen);
                translatedSeqBuff[aaLen] = '\n';
                fastaAaWriter.writeData(translatedSeqBuff, (aaLen + 1), 0, thread_idx, false, false);
                fastaCodonWriter.writeData(result.c_str(), result.size(), 0, thread_idx, false, false);
            }
            // move to another contig:
            plusPred.clearPred();
            minusPred.clearPred();
        }
        free(translatedSeqBuff);
    }
    fastaAaWriter.close(true);
    fastaCodonWriter.close(true);
    FileUtil::remove(fastaCodonFileNameIndex.c_str());
    FileUtil::remove(fastaAaFileNameIndex.c_str());

    mapWriter.close(true);
    FileUtil::remove(mapFileNameIndex.c_str());

    gffWriter.close(true);
    FileUtil::remove(gffFileNameIndex.c_str());
 
    contigsData.close();
    contigsHeaders.close();
    targetsHeaders.close();
    predsPerContig.close();

    return EXIT_SUCCESS;
}



