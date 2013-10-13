/* Call mutations from BAM file */
#include "bamalignfunc.h"
#include "parseargs.h"
#include "mutmap.h"
#include "givenrefmut.h"
#include "refio.h"

#include <fstream>
#include <sstream>
#include <climits>

using namespace std;



/*
Get the reference annotation
*/
RefVector RVREF;
int getReferenceInfo(BamReader & reader){

  RVREF=reader.GetReferenceData();
  return 0;
}

void printMutation(string refname,long pos,MutInfo& mti){
  string sep="\t";
  cout<<refname<<"\t"<<pos<<"\t"<<mti.ref<<"\t"<<mti.alt<<"\t"
      <<mti.refF<<sep<<mti.refB<<sep<<mti.altF<<sep<<mti.altB<<endl;
}

/*
Update the stored mutational position, exporting and removing all records before the mappos
*/
int updateMutPos(MutMap& mtmp, long mappos, string refseq, CallingArgs& args){
  mtmp.updateRefRecord(mappos);
  // print information for mut records before this mutation 
  for(MAPVMUT_ITR itr=mtmp.mvm.begin(); itr!=mtmp.mvm.end();itr++){
    if(itr->first >= mappos) break;
    //cout<<"Processing "<<itr->first<<endl;
    vector<MutInfo> & vmpt= itr->second;
    for(int j=0;j<vmpt.size();j++){
      int totalalt=vmpt[j].altF+vmpt[j].altB;
      bool outputmut=false;
      if(args.mutation_given) outputmut=true;
      else{
        if(totalalt>=args.min_read){
          if(args.printn || vmpt[j].alt != "N") outputmut=true;
        }
        if(args.skipindel && (vmpt[j].ref.size()!=1 || vmpt[j].alt.size()!=1) ) outputmut=false;
      }
      if(outputmut)
        printMutation(refseq,itr->first,vmpt[j]);
    }
  }
  // update the map information
  //cout<<"Removing before "<<mappos<<endl;
  mtmp.removeMutBeforePos(mappos);
}

map<string,MutMap> GIVENMUT;

int main(int argc, char* argv[]){

  //parse command line information
  CallingArgs args;
  if(parseArguments(argc,argv,args)==-1) return -1;
  
  // check if a given list of mutations is given; load chromosome names
  bool mutation_given=args.mutation_given;
  if(mutation_given){
    if(loadGivenMutationChrNames(GIVENMUT,args.mutfile)==-1){
      return -1;
    }
  }

  //BAM readers
  BamReader reader;
  string filename=string(args.bamfilename);
  if (!reader.Open(filename)){
    cerr<<"Could not open input BAM files."<<endl;
    return -1;
  }
  if(! reader.LocateIndex() || !reader.HasIndex()){
    cerr<<"Error loading index.\n";
    return -1;
  }
  
  // test for reading references
  getReferenceInfo(reader);
    
  // test for reading alignments
  BamAlignment al;
  long counter=0;
  long prevpos=-1;
  int prevrefid=-1;


  MutMap mtmp_0;
  MutMap* mtmp_itr=& mtmp_0;
  vector<NMStruct> vnms; //store the mismatches in one alignment
  set<long> blackout;
  REFREC thisrefrec;

  while(reader.GetNextAlignment(al)){
    counter++;
    if(counter%1000000 ==1)cerr<<counter<<"...\n";
    long mappos=(long) al.Position;
    int refid=(int)al.RefID;
    
    if(refid<0 || refid >=RVREF.size())continue;
    // update stored mutations before the current position
    if( prevrefid!=refid){
      string refidstr=RVREF[refid].RefName;
      if(prevrefid!=-1){
        updateMutPos(*mtmp_itr,LONG_MAX,RVREF[prevrefid].RefName,args);
      }
      // update reference
      if(mutation_given && GIVENMUT.count(refidstr)>0){
        // if we have unused reference before?
        for(map<string,MutMap>::iterator givenmut_itr=GIVENMUT.begin();
           givenmut_itr!=GIVENMUT.end();givenmut_itr++){
           if(givenmut_itr->first >= refidstr) break;
           if(givenmut_itr->first <= RVREF[prevrefid].RefName) continue;
           //For the unused reference, load the corresponding mutation list first
           loadGivenMutations(GIVENMUT, args.mutfile, givenmut_itr->first);
           if( (givenmut_itr->second).mvm.size()>0){
             updateMutPos(givenmut_itr->second, LONG_MAX, givenmut_itr->first,args);
           }
        }
        //load 
        loadGivenMutations(GIVENMUT, args.mutfile, refidstr);
        mtmp_itr=&(GIVENMUT[refidstr]);
      }else{
        mtmp_itr=&mtmp_0;
        mtmp_0.clear();
      }
      prevrefid=refid;
    }else{
      if(mappos>prevpos){ 
        updateMutPos(*mtmp_itr,mappos,RVREF[refid].RefName,args);
      }
    }
    prevpos=mappos;
    vnms.clear();
    blackout.clear();
    
    // check if we need to skip this read
    if(!al.IsPrimaryAlignment() || (al.IsPaired() && ! al.IsProperPair() ))continue;
    vector<CigarOp> & cgo= al.CigarData;
    // check if the read contains indels?
    if(args.skipindelread){
      bool hasindel=false;
      for(int i=0;i<cgo.size();i++){
        if(cgo[i].Type=='I' || cgo[i].Type=='D') {hasindel=true;break;}
      }
      if(hasindel) continue;
    }
  
    int getmismatchret=0;
    if(!args.usemdtag)
      getmismatchret=getMismatchInfoWithRefSeq(al,vnms,RVREF[refid].RefName,false);
    else
      getmismatchret=getMismatchInfo(al,vnms,false);
    if(getmismatchret!=0){
      cerr<<"Error in line "<<counter<<endl;
    }
    //string mdtag; al.GetTag("MD",mdtag);
    //cout<<"SEQ:"<<al.Name<<", MD tag:"<<mdtag<<",length"<<al.Length<<", CIGAR:"; for(int i=0;i<cgo.size();i++) cout<<cgo[i].Length<<cgo[i].Type<<" "; cout<<endl; 
    //cout<<"dup:"<<al.IsDuplicate()<<",failed:"<<al.IsFailedQC()<<",Mapped:"<<al.IsMapped()<<",second:"<<al.IsSecondMate()<<",MateMapped:"<<al.IsMateMapped()
    //  <<",MateRS:"<<al.IsMateReverseStrand()<<",paired:"<<al.IsPaired()<<",primary:"<<al.IsPrimaryAlignment()<<",proper:"<<al.IsProperPair()<<endl;
    
    for(int i=0;i<vnms.size();i++){
      vnms[i].chr_id=refid;
      blackout.insert(vnms[i].real_pos);
      if(vnms[i].relativepos < args.mut_span || al.Length- vnms[i].relativepos < args.mut_span ) continue;
      bool isrev=! al.IsReverseStrand();
      // compare with given reference sequence, check for valid
      if( args.has_fasta && vnms[i].origin.size()==1 && vnms[i].sub.size()==1){
        string mutinrefseq;
        int rfsqrs=refseq_getseq(RVREF[refid].RefName,vnms[i].pos-1,vnms[i].origin.size(),mutinrefseq); // in ref chr, it is 0-base, and in vcf it is 1-base
        if( rfsqrs==0 && vnms[i].origin != mutinrefseq){
          cerr<<"Error: incorrect REF sequence in "<<RVREF[refid].RefName<<":"<<vnms[i].pos<<" (real:"<<vnms[i].real_pos<<"), given "<<vnms[i].origin<<", should be "<<mutinrefseq<<". Check the reference genome or MD tag of your BAM file."<<endl;
        }
      }
      mtmp_itr->addOneMut( vnms[i].real_pos, vnms[i].origin, vnms[i].sub, isrev,!mutation_given,1);
    }
    
    // update ref alignment
    //cout<<"Blackout:"; for(set<long>::iterator bit=blackout.begin();bit!=blackout.end();bit++) cout<<*bit<<" ";cout<<endl;
    thisrefrec.start_pos=mappos+1; // 0-base to 1-base  
    thisrefrec.end_pos=al.GetEndPosition();
    thisrefrec.CigarData=cgo;
    thisrefrec.blackout=blackout;
    thisrefrec.direction=!al.IsReverseStrand();
    mtmp_itr->addOneRefAlign(thisrefrec);

  }
  if(prevrefid>=0 && prevrefid <RVREF.size()){
    if(!mutation_given){
      updateMutPos(*mtmp_itr,LONG_MAX,RVREF[prevrefid].RefName,args);
    }else{
        // if we have unused reference before?
        for(map<string,MutMap>::iterator givenmut_itr=GIVENMUT.begin();
           givenmut_itr!=GIVENMUT.end();givenmut_itr++){
           if(givenmut_itr->first < RVREF[prevrefid].RefName) continue;
           //For the unused reference, load the corresponding mutation list first
           if(givenmut_itr->first > RVREF[prevrefid].RefName)
             loadGivenMutations(GIVENMUT, args.mutfile, givenmut_itr->first);
           if( (givenmut_itr->second).mvm.size()>0){
             updateMutPos(givenmut_itr->second, LONG_MAX, givenmut_itr->first,args);
           }
        }//end for
    }//end if
  }//end if
  reader.Close();
  if(args.has_fasta){
    refseq_destroy();
  }
  return 0;
}


