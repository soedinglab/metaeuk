#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "unitetoseqdbs.sh.h"

int unitetoseqdbs(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, 6);

    // check if temp dir exists and if not, try to create it:
    if (FileUtil::directoryExists(par.db6.c_str()) == false){
        Debug(Debug::INFO) << "Temporary folder " << par.db6 << " does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db6.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << par.db6 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created directory " << par.db6 << "\n";
        }
    }

    CommandCaller cmd;
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    FileUtil::writeFile(par.db6 + "/unitetoseqdbs.sh", unitetoseqdbs_sh, unitetoseqdbs_sh_len);
    std::string program(par.db6 + "/unitetoseqdbs.sh");
    cmd.execProgram(program.c_str(), par.filenames);


    return EXIT_SUCCESS;
}
