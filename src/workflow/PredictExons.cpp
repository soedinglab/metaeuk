#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "predictexons.sh.h"

void setPredictExonsDefaults(Parameters *p) {
    p->orfStartMode = 1;
    // minimal exon length in codons:
    p->orfMinLength = 15;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    // evalue for search is high by default
    // The metaeuk Evalue Thr is lower
    p->evalThr = 100;
}

int predictexons(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    setPredictExonsDefaults(&par);
    par.parseParameters(argc, argv, command, false, 0, 0);
    int targetDbType = FileUtil::parseDbType(par.db2.c_str());
    if (Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::INFO) << "Enforcing exhaustive profile search mode due to profile target database\n";
        par.lcaSearch = true;
    }
    par.printParameters(command.cmd, argc, argv, *command.params);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("REVERSE_FRAGMENTS", par.reverseFragments == 1 ? "TRUE" : NULL);
    cmd.addVariable("EXTRACTORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    // align module should return alignments of at least a minimal exon length
    par.alnLenThr = par.minExonAaLength;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    cmd.addVariable("THREAD_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("COLLECTOPTIMALSET_PAR", par.createParameterString(par.collectoptimalset).c_str());

    std::string program(tmpDir + "/predictexons.sh");
    FileUtil::writeFile(program, predictexons_sh, predictexons_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // should never get here
    return EXIT_FAILURE;
}
