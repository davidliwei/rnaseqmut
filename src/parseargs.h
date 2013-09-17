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

  /* To output 'N's as alternative option? */
  bool printn;

  /* The reference fasta file */
  std::string ref_fasta;
  bool has_fasta; //whether the fasta file is given
  
  /* Whether to skip indels */
  bool skipindel;
  /* Whether to skip reads with indels */
  bool skipindelread;

  /* Whether to use MD tag to call mutation */
  bool usemdtag;
};

/* parsing arguments */
int parseArguments(int argc, char* argv[],CallingArgs& args);

#endif
