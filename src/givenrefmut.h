#ifndef GIVENREFMUT_H
#define GIVENREFMUT_H

#include "bamalignfunc.h"
#include "mutmap.h"

/* Loading all given mutations into memory by a given chromosome name
*/
int loadGivenMutations(map<string,MutMap> & givenmut, string filename, string givenchr);


/* Loading all possible chromosome names in given mutations into memory by a given chromosome name
*/
int loadGivenMutationChrNames(map<string,MutMap> & givenmut, string filename );

#endif
