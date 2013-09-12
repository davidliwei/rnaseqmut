
#ifndef MUTMAP_H
#define MUTMAP_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <list>
#include "bamalignfunc.h"

using namespace std;


/* a mutation information */
class MutInfo{
/* class members */
public:
  /* the basic information: position (stored somewhere else), origin seq and alternative seq */
  string ref;
  string alt;

  /* read counts in ref, alt, separating forward and backward */
  int refF;
  int refB;
  int altF;
  int altB;


/* class functions */
public:
  void clean(){
    refF=refB=altF=altB=0;
    ref=alt="";
  }
  MutInfo(){ clean(); }
};

/* define a ref record, used to calibrate the reference read counts */
struct REFREC{
  long start_pos;
  long end_pos;
  bool direction; // mapped direction, true for forward, false for backward 
  vector<CigarOp> CigarData;
  set<long> blackout;
};



/* all possible mutations in a map */
typedef map<long, vector<MutInfo> > MAPVMUT;
typedef map<long, vector<MutInfo> >::iterator  MAPVMUT_ITR;

class MutMap{
/* class members */
public:
  MAPVMUT mvm;
  // the following structures store the bam alignment to correct reference read counts
  list<REFREC> allrefrec;

/* class functions */
  void clear(){ mvm.clear();allrefrec.clear();}
  /* Add one mutation specified by pos:ref->alt 
  Paremeters:
    pos,ref,alt: the mutation information
    forward: whether this read is forward (true) or reverse (false)
    createifnoexist: whether to create a new record if this mutation does not exist
    increment: how many reads of this type should be added?
  */
  int addOneMut(long pos, string ref, string alt, bool forward, bool createifnoexist,int increment);

  /* remove all records, including mutations and ref sequences, before pos */
  int removeMutBeforePos(long pos);

  /* add the ref record for all mutations, except those specified in the blackout */ 
  int addOneRefAlign(REFREC& rfrc);
private:
  int addOneRefAlign(long read0, vector<CigarOp>& cgo,bool forward, set<long>& blackout);
public:
  
  /* update the REF record before the pos */
  int updateRefRecord(long pos);
};


#endif
