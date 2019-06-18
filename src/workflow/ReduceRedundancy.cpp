#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "reduceredundancy.sh.h"

int reduceredundancy(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, 4);

    // check if temp dir exists and if not, try to create it:
    if (FileUtil::directoryExists(par.db4.c_str()) == false){
        Debug(Debug::INFO) << "Temporary folder " << par.db4 << " does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created directory " << par.db4 << "\n";
        }
    }

    CommandCaller cmd;
    cmd.addVariable("SWAPDB_PAR", par.createParameterString(par.swapdb).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    FileUtil::writeFile(par.db4 + "/reduceredundancy.sh", reduceredundancy_sh, reduceredundancy_sh_len);
    std::string program(par.db4 + "/reduceredundancy.sh");
    cmd.execProgram(program.c_str(), par.filenames);


    return EXIT_SUCCESS;
}
