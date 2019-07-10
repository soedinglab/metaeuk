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
    par.parseParameters(argc, argv, command, true, 0, 0);

    // check if temp dir exists and if not, try to create it:
    if (FileUtil::directoryExists(par.db4.c_str()) == false){
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false){
            Debug(Debug::ERROR) << "Could not crate tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.predictexonsworkflow);
    std::string tmpDir = par.db4 + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str())==false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("REVERSE_FRAGMENTS", par.reverseFragments == 1 ? "TRUE" : NULL);
    cmd.addVariable("EXTRACTORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    // align module should return alignments of at least a minimal exon length
    par.alnLenThr = par.minExonAaLength;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    float originalEvalThr = par.evalThr;
    par.evalThr = std::numeric_limits<float>::max();
    cmd.addVariable("SWAPRESULT_PAR", par.createParameterString(par.swapresult).c_str());
    par.evalThr = originalEvalThr;
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    cmd.addVariable("COLLECTOPTIMALSET_PAR", par.createParameterString(par.collectoptimalset).c_str());

    FileUtil::writeFile(par.db4 + "/predictexons.sh", predictexons_sh, predictexons_sh_len);
    std::string program(par.db4 + "/predictexons.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
