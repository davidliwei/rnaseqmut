#ifndef BAMALIGNFUNC_H
#define BAMALIGNFUNC_H

#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;
using namespace BamTools;

/*
Defining a mismatch, including insertion, deletion, and substitution
*/
struct NMStruct{
  char type;  // can be either I, D, S, or U (undefined)
  
  /* the original and substituted basepair (s)
     for substitutions, they are strings of a single chacter;
     for insertions and deletions, their length may be greater than 1.
     For insertion and deletions, strings may extend from both ends to include information of repeats.
  */
  string origin; //the original basepair (i.e., in the reference)
  string sub; // the substituted basepair (i.e., in the reads)
  int chr_id;

  /*the position of the events; for insertion and deletion, this position is where the first letter of origin begins (as in the VCF file)
    the real position of the event is located at the real_pos 
  */
  long pos; 
  long real_pos;
  /* the length of the event.
  for substitution, it is always 1;
  for indels, these are the length of the insertion/deletion
  */
  int len; 
  // the following are internal information
  int relativepos; // the relative position of this event in the read, this is 1-base inclusive
  
  /* For insertion and deletion, fj and bj records the number of basepairs that are extended in front of the indel and at the end of the endel.
  */
  int fj;
  int bj;
  
  /* reset all the members of the struct */
  void clear();

  NMStruct(){
    clear();
  }

};

/* retrieve the mismatch info from a bam alignment
Parameter:
  al:	the BamAlignment structure
  vnms: the vector of NMStruct to store the results
  printdbginfo: the parameter to print debug information
Return value: 
  0 if success, -1 if error occurs
*/
int getMismatchInfo(BamAlignment & al, vector<NMStruct>& vnms, bool printdbginfo=false);
 


/* 
Check whether the given position is  within the range of a read
Input:
    read0: the beginning of the read
    cgo:   the parsed CIGAR option
Return value:
 -1 if it is outside the range;
 otherwise, return the relative position in a read.
NOTE: insertions are not counted in the MDTag, but are corrected in the posafterins .
*/
int posInRead(long read0,  vector<CigarOp>& cgo, long pos );

#endif
