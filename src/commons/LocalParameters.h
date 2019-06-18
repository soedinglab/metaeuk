#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

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

    std::vector<MMseqsParameter*> predictexonsworkflow;
    std::vector<MMseqsParameter*> collectoptimalset;
    std::vector<MMseqsParameter*> reduceredundancyworkflow;
    std::vector<MMseqsParameter*> assigntaxonomyworkflow;
    std::vector<MMseqsParameter*> unitetoseqdbsworkflow;

    PARAMETER(PARAM_REVERSE_FRAGMENTS)
    int reverseFragments;

    PARAMETER(PARAM_METAEUK_EVAL_THR)
    float metaeukEvalueThr;

    PARAMETER(PARAM_MAX_INTRON_LENGTH)
    size_t maxIntronLength;

    PARAMETER(PARAM_MIN_INTRON_LENGTH)
    size_t minIntronLength;

    PARAMETER(PARAM_MIN_EXON_AA_LENGTH)
    size_t minExonAaLength;

    PARAMETER(PARAM_MAX_AA_OVERLAP)
    size_t maxAaOverlap;

    PARAMETER(PARAM_GAP_OPEN_PENALTY)
    int setGapOpenPenalty;

    PARAMETER(PARAM_GAP_EXTEND_PENALTY)
    int setGapExtendPenalty;

private:
    LocalParameters() : 
        Parameters(),
        PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
        PARAM_METAEUK_EVAL_THR(PARAM_METAEUK_EVAL_THR_ID,"--metaeuk-eval", "maximal combined evalue of an optimal set", "maximal combined evalue of an optimal set [0.0, inf]", typeid(float), (void *) &metaeukEvalueThr, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$"),
        PARAM_MAX_INTRON_LENGTH(PARAM_MAX_INTRON_LENGTH_ID,"--max-intron", "Maximal intron length", "Maximal allowed intron length", typeid(int), (void *) &maxIntronLength, "^[0-9]+$"),
        PARAM_MIN_INTRON_LENGTH(PARAM_MIN_INTRON_LENGTH_ID,"--min-intron", "Minimal intron length", "Minimal allowed intron length", typeid(int), (void *) &minIntronLength, "^[0-9]+$"),
        PARAM_MIN_EXON_AA_LENGTH(PARAM_MIN_EXON_AA_LENGTH_ID,"--min-exon-aa", "Minimal exon length aa", "Minimal allowed exon length in amino acids", typeid(int), (void *) &minExonAaLength, "^[0-9]+$"),
        PARAM_MAX_AA_OVERLAP(PARAM_MAX_AA_OVERLAP_ID,"--max-overlap", "Maximal overlap of exons", "Maximal allowed overlap of consecutive exons in amino acids", typeid(int), (void *) &maxAaOverlap, "^[0-9]+$"),
        PARAM_GAP_OPEN_PENALTY(PARAM_GAP_OPEN_PENALTY_ID,"--set-gap-open", "Gap open penalty", "Gap open penalty (negative) for missed target amino acids between exons", typeid(int), (void *) &setGapOpenPenalty, "^-[0-9]+$"),
        PARAM_GAP_EXTEND_PENALTY(PARAM_GAP_EXTEND_PENALTY_ID,"--set-gap-extend", "Gap extend penalty", "Gap extend penalty (negative) for missed target amino acids between exons", typeid(int), (void *) &setGapExtendPenalty, "^-[0-9]+$")
    {
        collectoptimalset.push_back(&PARAM_METAEUK_EVAL_THR);
        collectoptimalset.push_back(&PARAM_MAX_INTRON_LENGTH);
        collectoptimalset.push_back(&PARAM_MIN_INTRON_LENGTH);
        collectoptimalset.push_back(&PARAM_MIN_EXON_AA_LENGTH);
        collectoptimalset.push_back(&PARAM_MAX_AA_OVERLAP);
        collectoptimalset.push_back(&PARAM_GAP_OPEN_PENALTY);
        collectoptimalset.push_back(&PARAM_GAP_EXTEND_PENALTY);
        collectoptimalset.push_back(&PARAM_SCORE_BIAS);
        collectoptimalset.push_back(&PARAM_THREADS);
        collectoptimalset.push_back(&PARAM_V);

        // predictexonsworkflow = combineList(extractorfs, translatenucs); // available through searchworkflow
        predictexonsworkflow = combineList(predictexonsworkflow, searchworkflow);
        predictexonsworkflow = combineList(predictexonsworkflow, collectoptimalset);
        predictexonsworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);

        reduceredundancyworkflow = combineList(reduceredundancyworkflow, swapdb);
        reduceredundancyworkflow.push_back(&PARAM_THREADS);
        reduceredundancyworkflow.push_back(&PARAM_REMOVE_TMP_FILES);

        assigntaxonomyworkflow.push_back(&PARAM_THREADS);
        assigntaxonomyworkflow.push_back(&PARAM_REMOVE_TMP_FILES);

        unitetoseqdbsworkflow.push_back(&PARAM_THREADS);
        unitetoseqdbsworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        
        // default value 0 means no reverse of AA fragments
        reverseFragments = 0;

        // default values for an optimal set:
        metaeukEvalueThr = 0.001;
        maxIntronLength = 10000;
        minIntronLength = 15;
        minExonAaLength = 11;
        maxAaOverlap = 10; // should be smaller than minExonAaLength
        setGapOpenPenalty = -1;
        setGapExtendPenalty = -1;
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
