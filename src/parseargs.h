/* Arguments passing in this program */
#ifndef PARSEARGS_H
#define PARSEARGS_H

#include <string>

struct CallingArgs{
  /* BAM file name */
  std::string bamfilename;

  /* For calling mutations */
  int mut_span;
  /* minimum read count to be output */
  int min_read;
  /* The filename of a given mutation list */
  std::string mutfile;
  bool mutation_given;
};

/* parsing arguments */
int parseArguments(int argc, char* argv[],CallingArgs& args);

#endif
