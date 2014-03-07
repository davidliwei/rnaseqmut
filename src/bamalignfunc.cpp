/* 
bam  alignment IOs 
NOTE: all positions are 0-based. VCF files are 1-based.
*/

#include "bamalignfunc.h"
#include "refio.h"

void NMStruct::clear(){

    type='U';
    origin=""; sub="";
    pos=0;real_pos=0;
    len=0; relativepos=0;
    fj=0;bj=0;
    chr_id=-1;
}


/* 
Get the relative and absolute position of the insertion, together with the length of the insertion from the CIGAR
Input:
    read0: the beginning of the read
    cgo:   the parsed CIGAR option
Return value:
    abspos: the absolute position,located exactly on the base before insertion begins
    relativepos: the relative position; (1-base)
    inslength: the length of the insertion
NOTE: insertions are not counted in the MDTag, but are corrected in the posafterins .
*/
int getInsertPos(long read0,  vector<CigarOp>& cgo, 
  vector<long>& abspos, vector<int>& relativepos, vector<int>& inslength  ){
  long read1=read0;
  int relpos=0;//record the current positions in the CIGAR to process so far
  int tojump=0;//record the positions to jump over
  for(int i=0;i<cgo.size();i++){
    switch(cgo[i].Type){
      case 'M': // match
        relpos+=cgo[i].Length;
        break;
      case 'N': // splicing junctions
        tojump+=cgo[i].Length;
        break;
      case 'I': // insertion relative to the reference; do not count towards the relative positions
        read1=read0+tojump+relpos;
        abspos.push_back(read1);
        relativepos.push_back(relpos+1);
        inslength.push_back(cgo[i].Length);
        relpos+=cgo[i].Length;
        break;
      case 'D': // deletion relative to the reference; jump 1 base
        tojump+=cgo[i].Length;
        break;
      case 'H':
      case 'S':
      case 'P':
      case '=':
      case 'X':
        //cerr<<"Error: CIGAR operator H, S, P, = and X are not supported currently.\n";
        return -1;
      default:
        //cerr<<"Error: unrecognized CIGAR character "<<cgo[i].Type<<endl;
        return -1;
    }
  }
  return abspos.size();
}


/* 
Convert the relative position in a read (pos, 1-based) to an absolute position, given the CIGAR character 
Return value:
    The real position;
NOTE: insertions are not counted in the MDTag, but are corrected in the posafterins .
*/
long getChrRealPos(long read0, int pos, vector<CigarOp>& cgo, int& posafterins){
  long read1=read0;
  int relpos=0;//record the current positions in the CIGAR to process so far
  int tojump=0;//record the positions to jump over
  bool quitfor=false;
  posafterins=pos;
  for(int i=0;i<cgo.size();i++){
    switch(cgo[i].Type){
      case 'M': // match
        if(relpos+cgo[i].Length >= pos) quitfor=true;
        relpos+=cgo[i].Length;
        break;
      case 'N': // splicing junctions
        tojump+=cgo[i].Length;
        break;
      case 'I': // insertion relative to the reference; do not count towards the relative positions
        relpos+=cgo[i].Length;
        posafterins+=cgo[i].Length; // correct for pos considering insertions
        break;
      case 'D': // deletion relative to the reference; jump 1 base
        tojump+=cgo[i].Length;
        break;
      case 'H':
      case 'S':
      case 'P':
      case '=':
      case 'X':
        //cerr<<"Error: CIGAR operator H, S, P, = and X are not supported currently.\n";
        return -1;
      default:
        cerr<<"Error: unrecognized CIGAR character "<<cgo[i].Type<<endl;
        return -1;
    }
    if(quitfor)break;
  }
  return read0+pos+tojump;
}

/* functions to parse substitutions and deletions from MD Tag 
Parameters:
  al:	  the BamAlignment record
  tag:	  the MD tag
Return value:
  subpos: the position of the substitution (relative to the read start)
  subchar:the character of the sub (i.e., the ref character). If the beginning of the character is ^, then this is an deletion event
*/
int parseSubFromMDTag(BamAlignment & al, string tag, vector<int>& subpos, vector<string>& subchar){
  stringstream ss(tag);
  int pos=-1;
  while(!ss.eof()){
    ss>>pos; 
    if(ss.eof())break;
    char c;
    string afterpos="";
    do{
      ss>>c;
      afterpos.push_back(c);
    }while(!ss.eof() && (c<'0' || c >'9'));
    ss.putback(c);
    //pos and afterpos[-1] are positions, sub chars
    if(afterpos.size()>0){
      subpos.push_back(pos);
      subchar.push_back(string(afterpos.substr(0,afterpos.size()-1)));
    }
  }
  return 0;
}

/* retrieve the mismatch info from a bam alignment
Parameter:
  al:	the BamAlignment structure
  nmsv: the vector of NMStruct to store the results
  printdbginfo: the parameter to print debug information
Return value: 
  0 if success, -1 if error occurs
*/
int getMismatchInfo(BamAlignment & al, vector<NMStruct>& nmsv, bool printdbginfo){
  long mappos=(long) al.Position;
  int refid=(int)al.RefID;
  //uint32_t nmtag=0;
  int32_t nmtag_i=0;
  int nmtag=0;
  if(!al.HasTag("NM")){
    cerr<<"Error: no NM tag.\n";
    return -1;
  }
  if( !al.GetTag("NM",nmtag_i)){
     // switch to uint32_t for some versions of BAM
     uint32_t nmtag_u=0;
     if( !al.GetTag("NM",nmtag_u)){
       cerr<<"Error reading NM tags in the read alignment; skipping the reads. NMtag="<<nmtag<<".\n";
       return -1;
     }else{
       nmtag=nmtag_u;
     }
  }else{
    nmtag=nmtag_i;
  }
  if(nmtag<1) return 0;
  //MD type
  string mdtag;
  if(!al.HasTag("MD") ||!al.GetTag("MD",mdtag) ){
    cerr<<"Error reading MD tags in the read alignment; skipping the reads..\n";
    return -1;
  }
  // CIGAR
  vector<CigarOp> & co=al.CigarData;
  
  // parse all substitutions
  vector<int> subpos;
  vector<string> subchar;
  parseSubFromMDTag(al, mdtag,subpos,subchar);
  // convert to abs positions
  int relpos=0;
  int relposafterins=0;
  // processing substitutions and deletions
  for(int i=0;i<subpos.size();i++){
    relpos+=subpos[i]; // move 1 bp forward to exactly where the substitution happens
    if(subchar[i][0]!='^'){
      // substitutions
      relpos++; // this is the exact position of the substitution
      long realpos=getChrRealPos(mappos,relpos,co,relposafterins); // the mapped real positions, right now matches exactly the VCF standards
      char originbase=al.QueryBases[relposafterins-1];

      NMStruct nms;
      nms.type='S'; nms.pos=realpos; nms.real_pos=realpos; //here, realpos should be 1-based 
      nms.origin=subchar[i]; nms.sub.assign(1,originbase);
      nms.relativepos=relposafterins;
      nms.len=1;
      nmsv.push_back(nms);

      //cout<<"SUB: "<<subpos[i]<<"("<<subchar[i]<<"->"<<originbase<<"), real position:"<<realpos<<"; relative pos:"<<relpos<<", after insertion:"<<relposafterins<<endl;
    }else{
      // deletions
    
      long realpos=getChrRealPos(mappos,relpos,co,relposafterins); // the mapped real positions, right now matches exactly the VCF standards
      string insertedseq=subchar[i].substr(1);
       // the deleted sequence is right after this position
      string disstr="";
      int deljump=0;
      for(int j=relposafterins-1;j>=0;j--){
        disstr.insert(disstr.begin(),al.QueryBases[j]);
        if(al.QueryBases[j] != insertedseq[0]) break;
        deljump++;
        realpos--;
      }
      int forwardextend=0;
      for(int j=relposafterins;j<al.QueryBases.size();j++){
        if(al.QueryBases[j]!=insertedseq[insertedseq.size()-1]) break; 
        disstr.append(1,al.QueryBases[j]);
        forwardextend++;
      }
      
      NMStruct nms;
      nms.type='D'; nms.pos=realpos; nms.real_pos=realpos+deljump;
      nms.origin=disstr+insertedseq;
      nms.sub=disstr;
      nms.relativepos=relposafterins;
      nms.len=insertedseq.size();
      nms.bj=deljump; nms.fj=forwardextend;
      nmsv.push_back(nms);
      
      //cout<<"DEL: "<<subpos[i]<<"("<<disstr+insertedseq<<"->"<<disstr<<"), real position:"<<realpos<<"; relative pos:"<<relpos<<", after insertion:"<<relposafterins<<", backward jump:"<<deljump<<endl;
    }
  }
  
  //processing insertions
  vector<long> insabspos;
  vector<int>  insrelativepos;
  vector<int> inslength;
  if(getInsertPos(mappos,co,insabspos,insrelativepos,inslength)>0){
    for(int i=0;i<insabspos.size();i++){
      long iabspos=insabspos[i];
      int irelpos=insrelativepos[i];
      int ilen=inslength[i];
      string insseq=al.QueryBases.substr(irelpos-1,ilen); //the position is exactly where the insertion occurs, and need 1-base to 0-base conversion
      //search left and right hand side to extend
      string disstr="";
      string insseq2=insseq;
      int insbackwardjump=0;
      for(int j=irelpos-1-1;j>=0;j--){
        disstr.insert(disstr.begin(),al.QueryBases[j]);
        insseq.insert(insseq.begin(),al.QueryBases[j]);
        if(al.QueryBases[j] != insseq2[0]) break;
        insbackwardjump++;
        iabspos--;
      }
      int insforwardjump=0;
      for(int j=irelpos+ilen-1;j<al.QueryBases.size();j++){
        if(al.QueryBases[j]!=insseq2[insseq2.size()-1]) break; 
        disstr.append(1,al.QueryBases[j]);
        insseq.append(1,al.QueryBases[j]);
        insforwardjump++;
      }

      NMStruct nms;
      nms.type='I'; nms.pos=iabspos; nms.real_pos=iabspos+insbackwardjump;
      nms.origin=disstr;
      nms.sub=insseq;
      nms.relativepos=irelpos;
      nms.len=ilen;
      nms.bj=insbackwardjump; nms.fj=insforwardjump;
      nmsv.push_back(nms);

      //cout<<"INS: "<<irelpos<<"("<<disstr<<"->"<<insseq<<")"<<", real position:"<<iabspos<<", length:"<<ilen<<", backward jump:"<<insbackwardjump<<", forward jump:"<<insforwardjump<<endl;
    }
  }

  // print debug information?
  if(printdbginfo){
    for(int i=0;i<nmsv.size();i++){
      cout<<"NMS Type:"<<nmsv[i].type<<" ("<<nmsv[i].origin<<"->"<<nmsv[i].sub<<"), pos:"<<nmsv[i].pos<<",real_pos="<<nmsv[i].real_pos<<", len="<<nmsv[i].len<<", relative: "<<nmsv[i].relativepos<<"; fj="<<nmsv[i].fj<<", bj="<<nmsv[i].bj<<endl; 
    }
  }
  return 0;
}

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
int posInRead(long read0,  vector<CigarOp>& cgo, long pos ){
  int difpos=pos-read0;
  if(difpos<0) return -1;
  long read1=read0;
  int relpos=0;//record the current positions in the CIGAR to process so far
  int tojump=0;//record the positions to jump over
  for(int i=0;i<cgo.size();i++){
    switch(cgo[i].Type){
      case 'M': // match
        relpos+=cgo[i].Length;
        if(difpos-tojump<relpos){
          return difpos-tojump;
        }
        break;
      case 'N': // splicing junctions
        tojump+=cgo[i].Length;
        if(difpos-tojump<relpos){
          return -1;
        }
        break;
      case 'I': // insertion relative to the reference; do not count towards the relative positions
        read1=read0+tojump+relpos;
        relpos+=cgo[i].Length;
        break;
      case 'D': // deletion relative to the reference; jump 1 base
        tojump+=cgo[i].Length;
        if(difpos-tojump<relpos){
          return -1;
        }
        break;
      case 'H':
      case 'S':
      case 'P':
      case '=':
      case 'X':
        //cerr<<"Error: CIGAR operator H, S, P, = and X are not supported currently.\n";
        return -1;
      default:
        cerr<<"Error: unrecognized CIGAR character "<<cgo[i].Type<<endl;
        return -1;
    }
  }
  return -1;
}


/* retrieve the mismatch info from a bam alignment using given reference sequence, without using a MD tag
Parameter:
  al:	the BamAlignment structure
  nmsv: the vector of NMStruct to store the results
  printdbginfo: the parameter to print debug information
Return value: 
  0 if success, -1 if error occurs
*/
int getMismatchInfoWithRefSeq(BamAlignment & al, vector<NMStruct>& nmsv, string refidstr, bool printdbginfo){
  long mappos=(long) al.Position;
  int refid=(int)al.RefID;
  //uint32_t nmtag=0;
  int32_t nmtag_i=0;
  int nmtag=0;
  if(!al.HasTag("NM")){
    cerr<<"Error: no NM tag.\n";
    return -1;
  }
  if( !al.GetTag("NM",nmtag_i)){
     // switch to uint32_t for some versions of BAM
     uint32_t nmtag_u=0;
     if( !al.GetTag("NM",nmtag_u)){
       cerr<<"Error reading NM tags in the read alignment; skipping the reads. NMtag="<<nmtag<<".\n";
       return -1;
     }else{
       nmtag=nmtag_u;
     }
  }else{
    nmtag=nmtag_i;
  }
  if(nmtag<1) return 0;

  // CIGAR
  vector<CigarOp> & co=al.CigarData;
  string albase=al.QueryBases;

  int relposafterins=0;
  
  for(int n=1;n<=albase.size();n++){ // n uses 1 base
    long realpos=getChrRealPos(mappos,n,co,relposafterins); // the mapped real positions, right now matches exactly the VCF standards
    if(realpos==-1) return 0;
    string refseq;
    if(refseq_getseq(refidstr,realpos-1,1,refseq)!=0) continue; // to access the sequence, use 0-base
    string oriseq=albase.substr(n-1,1); // to access the substr, use 0-base
    if(oriseq!=refseq){
      NMStruct nms;
      nms.type='S'; nms.pos=realpos; nms.real_pos=realpos; //should be 1-base
      nms.origin=refseq; nms.sub=oriseq;
      nms.relativepos=relposafterins; // 1-base
      nms.len=1;
      nmsv.push_back(nms);
    }
  }
   // print debug information?
  if(printdbginfo){
    for(int i=0;i<nmsv.size();i++){
      cout<<"NMS Type:"<<nmsv[i].type<<" ("<<nmsv[i].origin<<"->"<<nmsv[i].sub<<"), pos:"<<nmsv[i].pos<<",real_pos="<<nmsv[i].real_pos<<", len="<<nmsv[i].len<<", relative: "<<nmsv[i].relativepos<<"; fj="<<nmsv[i].fj<<", bj="<<nmsv[i].bj<<endl; 
    }
  }
  return 0;
}


