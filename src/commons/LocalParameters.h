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

private:
    LocalParameters() : 
        Parameters(),
        PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
        PARAM_METAEUK_EVAL_THR(PARAM_METAEUK_EVAL_THR_ID,"--metaeuk-eval", "maximal combined evalue of an optimal set", "maximal combined evalue of an optimal set [0.0, inf]", typeid(float), (void *) &metaeukEvalueThr, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$")
    {
        collectoptimalset.push_back(&PARAM_METAEUK_EVAL_THR);
        collectoptimalset.push_back(&PARAM_THREADS);
        collectoptimalset.push_back(&PARAM_V);

        // predictexonsworkflow = combineList(extractorfs, translatenucs); // available through searchworkflow
        predictexonsworkflow = combineList(predictexonsworkflow, searchworkflow);
        predictexonsworkflow = combineList(predictexonsworkflow, collectoptimalset);
        predictexonsworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);

        reduceredundancyworkflow.push_back(&PARAM_THREADS);
        reduceredundancyworkflow.push_back(&PARAM_REMOVE_TMP_FILES);

        assigntaxonomyworkflow.push_back(&PARAM_THREADS);
        assigntaxonomyworkflow.push_back(&PARAM_REMOVE_TMP_FILES);

        unitetoseqdbsworkflow.push_back(&PARAM_THREADS);
        unitetoseqdbsworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        
        // default value 0 means no reverse of AA fragments
        reverseFragments = 0;

        // default value for an optimal set is 0.001
        metaeukEvalueThr = 0.001;
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
