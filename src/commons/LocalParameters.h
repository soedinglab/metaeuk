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

std::vector<MMseqsParameter> predictexonsworkflow;
std::vector<MMseqsParameter> reduceredundancyworkflow;

private:
    LocalParameters() : Parameters() {
        predictexonsworkflow = combineList(extractorfs, translatenucs);
        predictexonsworkflow = combineList(predictexonsworkflow, searchworkflow);
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
