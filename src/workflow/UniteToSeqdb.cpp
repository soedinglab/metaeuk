#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "unitetoseqdb.sh.h"

int unitetoseqdb(int argc, const char **argv, const Command& command) {
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
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("UNITE_EXONS_REP", par.uniteExonsRep == 1 ? "TRUE" : NULL);
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    FileUtil::writeFile(par.db4 + "/unitetoseqdb.sh", unitetoseqdb_sh, unitetoseqdb_sh_len);
    std::string program(par.db4 + "/unitetoseqdb.sh");
    cmd.execProgram(program.c_str(), par.filenames);


    return EXIT_SUCCESS;
}
