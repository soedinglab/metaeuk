#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "predictexons.sh.h"

int predictexons(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, 3);

    CommandCaller cmd;
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }

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

    FileUtil::writeFile(par.db4 + "/predictexons.sh", predictexons_sh, predictexons_sh_len);
    std::string program(par.db4 + "/predictexons.sh");
    cmd.execProgram(program.c_str(), par.filenames);


    return EXIT_SUCCESS;
}
