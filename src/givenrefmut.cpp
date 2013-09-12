/* Operations concerning a given list of mutations */
#include "givenrefmut.h"

/* Loading all given mutations into memory by a given chromosome name
*/
int loadGivenMutations(map<string,MutMap> & givenmut, string filename, string givenchr){
  ifstream ifs(filename.c_str());
  if(!ifs.is_open()){  
    cerr<<"Error opening mutation file "<<filename<<endl;
    return -1;
  }
  // clear all records in MutMap
  for(map<string,MutMap>::iterator mmit=givenmut.begin();mmit!=givenmut.end();mmit++)
    (mmit->second).clear();
  string oneline;
  int ncount=0;
  int ngoodcount=0;
  while(true){
    getline(ifs,oneline);
    ncount++;
    if( ifs.eof())break;
    if(oneline.size()==0)continue;
    if(oneline[0]=='#' || oneline[0]=='@')continue;
    // parsing 
    stringstream ss(oneline);
    string chrom;
    long position;
    string ref;
    string alt;
    ss>>chrom;
    ss>>position;
    ss>>ref;
    ss>>alt;
    if(ss.fail()){
      cerr<<"Error parsing line "<<ncount<<"; ignore this line.\n";
      continue;
    }
    if(chrom!=givenchr) continue;
    if(givenmut.count(chrom)==0)continue;
    ngoodcount++;
    MutMap& mm=givenmut[chrom];
    mm.addOneMut( position, ref, alt, true,true, 0);
  }
  cerr<<"Reading "<<ncount<<" lines, "<<givenchr<<": "<<ngoodcount<<" records."<<endl;
  ifs.close();
  return 0;
}


/* Loading all possible chromosome names in given mutations into memory by a given chromosome name
*/
int loadGivenMutationChrNames(map<string,MutMap> & givenmut, string filename ){
  ifstream ifs(filename.c_str());
  if(!ifs.is_open()){  
    cerr<<"Error opening mutation file "<<filename<<endl;
    return -1;
  }
  // clear all records in MutMap
  givenmut.clear();
  string oneline;
  int ncount=0;
  int ngoodcount=0;
  while(true){
    getline(ifs,oneline);
    ncount++;
    if( ifs.eof())break;
    if(oneline.size()==0)continue;
    if(oneline[0]=='#' || oneline[0]=='@')continue;
    // parsing 
    stringstream ss(oneline);
    string chrom;
    long position;
    string ref;
    string alt;
    ss>>chrom;
    ss>>position;
    ss>>ref;
    ss>>alt;
    if(ss.fail()){
      cerr<<"Error parsing line "<<ncount<<"; ignore this line.\n";
      continue;
    }
    if(givenmut.count(chrom)>0)continue;
    ngoodcount++;
    MutMap& mm=givenmut[chrom];
    mm.clear();
  }
  cerr<<"Reading "<<ncount<<" lines, "<<ngoodcount<<" chromosomes."<<endl;
  ifs.close();
  return 0;
}


