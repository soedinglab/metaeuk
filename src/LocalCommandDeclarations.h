#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int predictexons(int argc, const char **argv, const Command& command);
extern int collectoptimalset(int argn, const char **argv, const Command& command);
extern int unitesetstofasta(int argn, const char **argv, const Command& command);
extern int reduceredundancy(int argc, const char **argv, const Command& command);
extern int grouppredictions(int argn, const char **argv, const Command& command);
extern int assigntaxonomy(int argc, const char **argv, const Command& command);
extern int resultspercontig(int argc, const char **argv, const Command& command);
extern int swapdb(int argc, const char **argv, const Command& command);

#endif
