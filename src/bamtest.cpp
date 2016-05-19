
#include "bamalignfunc.h"

RefVector RVREF;

/*
Get the reference annotation
*/
int getReferenceInfo(BamReader & reader){
  RVREF=reader.GetReferenceData();
  return 0;
}


int main(int argc, char* argv[]){
  BamReader reader;
  string filename=string(argv[1]);
  if (!reader.Open(filename)){
    cerr<<"Could not open input BAM files."<<endl;
    return -1;
  }
  if(reader.LocateIndex()){
    cout<<"Successfully located index.\n";
  }
  if(reader.HasIndex()){
    cout<<"Index opened.\n";
  }
  // test for reading references
  getReferenceInfo(reader);
    
  // test for reading alignments
  BamAlignment al;
  long counter=0;
  while(reader.GetNextAlignment(al)){
    counter++;
    if(counter%10000==1)   cout<<"Reading "<<counter<<" records.\n";
    
    long mappos=(long) al.Position;
    int refid=(int)al.RefID;
    
    uint32_t nmtag=0;
    if(al.HasTag("NM")){
      if(!al.GetTag("NM",nmtag)){
        continue;
      }   
    }else{
      if(al.HasTag("nM")){
        if(!al.GetTag("nM",nmtag)){
          continue;
        }   
      }else{
          continue;
      }
    }
    if(nmtag<1) continue;
    //MD type
    string mdtag;
    if(al.HasTag("MD")){
      if(al.GetTag("MD",mdtag)){
        cout<<RVREF[refid].RefName<<":"<<mappos<<", MDtag:"<<mdtag<<", NM:"<<nmtag<<endl;
      }
    }
    // CIGAR
    vector<CigarOp> & co=al.CigarData;
    for(int i=0;i<co.size();i++) cout<<co[i].Type<<" "<<co[i].Length<<";"; cout<<endl;
    
    cout<<"Q:	"<<al.QueryBases<<endl;
    cout<<"A:	"<<al.AlignedBases<<endl;
    vector<NMStruct> vnms;
    getMismatchInfo(al,vnms,true);
  
  }
  reader.Close();
  return 0;
}


