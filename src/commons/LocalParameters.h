#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

const int CITATION_METAEUK = CITATION_END << 1;

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> taxpercontigworkflow;
    std::vector<MMseqsParameter*> easypredictworkflow;
    std::vector<MMseqsParameter*> predictexonsworkflow;
    std::vector<MMseqsParameter*> collectoptimalset;
    std::vector<MMseqsParameter*> reduceredundancy;
    std::vector<MMseqsParameter*> unitesetstofasta;

    PARAMETER(PARAM_REVERSE_FRAGMENTS)
    int reverseFragments;

    PARAMETER(PARAM_METAEUK_EVAL_THR)
    float metaeukEvalueThr;

    PARAMETER(PARAM_METAEUK_TARGET_COV_THR)
    float metaeukTargetCovThr;

    PARAMETER(PARAM_MAX_INTRON_LENGTH)
    size_t maxIntronLength;

    PARAMETER(PARAM_MIN_INTRON_LENGTH)
    size_t minIntronLength;

    PARAMETER(PARAM_MIN_EXON_AA_LENGTH)
    size_t minExonAaLength;

    PARAMETER(PARAM_MAX_AA_OVERLAP)
    size_t maxAaOverlap;

    PARAMETER(PARAM_MAX_EXON_SETS)
    size_t maxExonSets;

    PARAMETER(PARAM_GAP_OPEN_PENALTY)
    int setGapOpenPenalty;

    PARAMETER(PARAM_GAP_EXTEND_PENALTY)
    int setGapExtendPenalty;

    PARAMETER(PARAM_SHOULD_TRANSLATE)
    int shouldTranslate;

    PARAMETER(PARAM_ALLOW_OVERLAP)
    int overlapAllowed;

    PARAMETER(PARAM_WRITE_TKEY)
    int writeTargetKey;

    PARAMETER(PARAM_WRITE_FRAG_COORDS)
    int writeFragCoords;

    PARAMETER(PARAM_LEN_SCAN_FOR_START)
    int lenScanForStart;

    LocalParameters() : 
        Parameters(),
        PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
        PARAM_METAEUK_EVAL_THR(PARAM_METAEUK_EVAL_THR_ID,"--metaeuk-eval", "maximal combined evalue of an optimal set", "maximal combined evalue of an optimal set [0.0, inf]", typeid(float), (void *) &metaeukEvalueThr, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$"),
        PARAM_METAEUK_TARGET_COV_THR(PARAM_METAEUK_TARGET_COV_THR_ID,"--metaeuk-tcov", "minimal length ratio between combined optimal set and target", "minimal length ratio of combined set to target [0.0, 1.0]", typeid(float), (void *) &metaeukTargetCovThr, "^0(\\.[0-9]+)?|^1(\\.0+)?$"),
        PARAM_MAX_INTRON_LENGTH(PARAM_MAX_INTRON_LENGTH_ID,"--max-intron", "Maximal intron length", "Maximal allowed intron length", typeid(int), (void *) &maxIntronLength, "^[0-9]+$"),
        PARAM_MIN_INTRON_LENGTH(PARAM_MIN_INTRON_LENGTH_ID,"--min-intron", "Minimal intron length", "Minimal allowed intron length", typeid(int), (void *) &minIntronLength, "^[0-9]+$"),
        PARAM_MIN_EXON_AA_LENGTH(PARAM_MIN_EXON_AA_LENGTH_ID,"--min-exon-aa", "Minimal exon length aa", "Minimal allowed exon length in amino acids", typeid(int), (void *) &minExonAaLength, "^[0-9]+$"),
        PARAM_MAX_AA_OVERLAP(PARAM_MAX_AA_OVERLAP_ID,"--max-overlap", "Maximal overlap of exons", "Maximal allowed overlap of consecutive exons in amino acids", typeid(int), (void *) &maxAaOverlap, "^[0-9]+$"),
        PARAM_MAX_EXON_SETS(PARAM_MAX_EXON_SETS_ID,"--max-exon-sets", "Maximal number of exon sets", "Maximal number of exons sets that match a given target. If >1 suboptimal solutions will be reported", typeid(int), (void *) &maxExonSets, "^[0-9]+$"),
        PARAM_GAP_OPEN_PENALTY(PARAM_GAP_OPEN_PENALTY_ID,"--set-gap-open", "Gap open penalty", "Gap open penalty (negative) for missed target amino acids between exons", typeid(int), (void *) &setGapOpenPenalty, "^-[0-9]+$"),
        PARAM_GAP_EXTEND_PENALTY(PARAM_GAP_EXTEND_PENALTY_ID,"--set-gap-extend", "Gap extend penalty", "Gap extend penalty (negative) for missed target amino acids between exons", typeid(int), (void *) &setGapExtendPenalty, "^-[0-9]+$"),
        PARAM_SHOULD_TRANSLATE(PARAM_SHOULD_TRANSLATE_ID,"--protein", "translate codons to AAs", "translate the joint exons coding sequence to amino acids [0,1]", typeid(int), (void *) &shouldTranslate, "^[0-1]{1}$"),
        PARAM_ALLOW_OVERLAP(PARAM_ALLOW_OVERLAP_ID,"--overlap", "allow same-strand overlaps", "allow predictions to overlap another on the same strand. when not allowed (default), only the prediction with better E-value will be retained [0,1]", typeid(int), (void *) &overlapAllowed, "^[0-1]{1}$"),
        PARAM_WRITE_TKEY(PARAM_WRITE_TKEY_ID,"--target-key", "write target key instead of accession", "write the target key (internal DB identifier) instead of its accession. By default (0) target accession will be written [0,1]", typeid(int), (void *) &writeTargetKey, "^[0-1]{1}$"),
        PARAM_WRITE_FRAG_COORDS(PARAM_WRITE_FRAG_COORDS_ID,"--write-frag-coords", "write fragment contig coords", "write the contig coords of the stop-to-stop fragment in which putative exon lies. By default (0) only putative exon coords will be written [0,1]", typeid(int), (void *) &writeFragCoords, "^[0-1]{1}$"),
        PARAM_LEN_SCAN_FOR_START(PARAM_LEN_SCAN_FOR_START_ID,"--len-scan-for-start", "length to scan for start codon", "length to scan for a start codon before the first exon and in the same frame. By default (0) no scan", typeid(int), (void *) &lenScanForStart, "^[0-9]+$")
        
    {
        collectoptimalset.push_back(&PARAM_METAEUK_EVAL_THR);
        collectoptimalset.push_back(&PARAM_METAEUK_TARGET_COV_THR);
        collectoptimalset.push_back(&PARAM_MAX_INTRON_LENGTH);
        collectoptimalset.push_back(&PARAM_MIN_INTRON_LENGTH);
        collectoptimalset.push_back(&PARAM_MIN_EXON_AA_LENGTH);
        collectoptimalset.push_back(&PARAM_MAX_AA_OVERLAP);
        collectoptimalset.push_back(&PARAM_MAX_EXON_SETS);
        collectoptimalset.push_back(&PARAM_GAP_OPEN_PENALTY);
        collectoptimalset.push_back(&PARAM_GAP_EXTEND_PENALTY);
        collectoptimalset.push_back(&PARAM_SCORE_BIAS);
        collectoptimalset.push_back(&PARAM_THREADS);
        collectoptimalset.push_back(&PARAM_COMPRESSED);
        collectoptimalset.push_back(&PARAM_V);

        // predictexonsworkflow = combineList(extractorfs, translatenucs); // available through searchworkflow
        predictexonsworkflow = combineList(searchworkflow, collectoptimalset);
        predictexonsworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);
        predictexonsworkflow = removeParameter(predictexonsworkflow, PARAM_NUM_ITERATIONS);
        predictexonsworkflow = removeParameter(predictexonsworkflow, PARAM_ADD_ORF_STOP);

        reduceredundancy.push_back(&PARAM_ALLOW_OVERLAP);
        reduceredundancy.push_back(&PARAM_THREADS);
        reduceredundancy.push_back(&PARAM_COMPRESSED);
        reduceredundancy.push_back(&PARAM_V);

        unitesetstofasta.push_back(&PARAM_SHOULD_TRANSLATE);
        unitesetstofasta.push_back(&PARAM_TRANSLATION_TABLE);
        unitesetstofasta.push_back(&PARAM_WRITE_TKEY);
        unitesetstofasta.push_back(&PARAM_WRITE_FRAG_COORDS);
        unitesetstofasta.push_back(&PARAM_LEN_SCAN_FOR_START);
        unitesetstofasta.push_back(&PARAM_MAX_SEQ_LEN);
        unitesetstofasta.push_back(&PARAM_THREADS);
        unitesetstofasta.push_back(&PARAM_V);

        easypredictworkflow = combineList(predictexonsworkflow, reduceredundancy);
        easypredictworkflow = combineList(easypredictworkflow, unitesetstofasta);

        taxpercontigworkflow = combineList(taxonomy, aggregatetax);
        
        // default value 0 means no reverse of AA fragments
        reverseFragments = 0;

        // default value 0 means no same-strand overlap of representative predictions is allowed
        overlapAllowed = 0;

        // default values for an optimal set:
        metaeukEvalueThr = 0.001;
        metaeukTargetCovThr = 0.5;
        maxIntronLength = 10000;
        minIntronLength = 15;
        minExonAaLength = 11;
        maxAaOverlap = 10; // should be smaller than minExonAaLength
        maxExonSets = 1; // by default, find just the optimal
        setGapOpenPenalty = -1;
        setGapExtendPenalty = -1;

        // default value 0 means no transaltion to AAs
        shouldTranslate = 0;

        // default value 0 means the accession is written
        writeTargetKey = 0;

        // default value 0 means only coords of putative exon are written
        writeFragCoords = 0;

        // default value 0 means no scanning before the first exon
        lenScanForStart = 0;

        citations.emplace(CITATION_METAEUK, "Levy Karin E, Mirdita M, Soeding J: MetaEuk – sensitive, high-throughput gene discovery and annotation for large-scale eukaryotic metagenomics. biorxiv, 851964 (2019).");
    }
private:
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
