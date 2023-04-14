#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "easypredict.sh.h"

void setEasyPredictDefaults(Parameters *p) {
    p->orfStartMode = 1;
    // minimal exon length in codons:
    p->orfMinLength = 15;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    // evalue for search is high by default
    // The metaeuk Evalue Thr is lower
    p->evalThr = 100;
}

int easypredict(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    setEasyPredictDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("PREDICTEXONS_PAR", par.createParameterString(par.predictexonsworkflow).c_str());
    cmd.addVariable("REDUCEREDUNDANCY_PAR", par.createParameterString(par.reduceredundancy).c_str());
    cmd.addVariable("UNITESETSTOFASTA_PAR", par.createParameterString(par.unitesetstofasta).c_str());
    cmd.addVariable("VERBOSITY_COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    cmd.addVariable("THREAD_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());

    std::string program("/bin/sh " + tmpDir + "/easypredict.sh");
    FileUtil::writeFile(program, easypredict_sh, easypredict_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // should never get here
    return EXIT_FAILURE;
}
