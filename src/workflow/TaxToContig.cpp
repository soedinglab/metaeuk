#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "taxtocontig.sh.h"

void setTaxToContigDefaults(Parameters *p) {
    p->majorityThr = 0.5;
    p->taxonomyOutputMode = 2;
}

int taxtocontig(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    setTaxToContigDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.db6;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    
    // the call to aggregatetaxweights requires that we output both the taxonomy and the alignment
    par.taxonomyOutputMode = 2;
    cmd.addVariable("TAXONOMY_PAR", par.createParameterString(par.taxonomy).c_str());
    cmd.addVariable("AGGREGATETAX_PAR", par.createParameterString(par.aggregatetax).c_str());
    cmd.addVariable("CREATETSV_PAR", par.createParameterString(par.createtsv).c_str());
    cmd.addVariable("VERBOSITY_COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    cmd.addVariable("THREAD_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());

    std::string program(tmpDir + "/taxtocontig.sh");
    FileUtil::writeFile(program, taxtocontig_sh, taxtocontig_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // should never get here
    return EXIT_FAILURE;
}
