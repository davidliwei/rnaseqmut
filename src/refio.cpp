/* Rading and writing reference sequences */
#include "refio.h"
#include "kseq.h"
#include <stdio.h>
#include <zlib.h>
#include <map>


/* Initialize refseq sequences */

//KSEQ_INIT(FILE, read);
KSEQ_INIT(gzFile, gzread);

gzFile refFile;
kseq_t *refseq=NULL;
string reffilename="";
string currentchr="";

map<string,long> refgenome;

/* buffer operation */
string ref_buff;
long buff_a=-1;
long buff_b=-1;
/* the window size for buffer */
long buff_a_size=100000;
long buff_b_size=10000000;


/* Finish loading the reference sequence */
int refseq_destroy(){
  // destroy
  if(refseq != NULL){
    kseq_destroy(refseq);
  }
  refseq=NULL;
  return 0;
}


/* Initialize the reference sequence */
int refseq_init(string filename){
  reffilename=filename; 

  refFile=gzopen(reffilename.c_str(),"r");
  if(refFile == NULL){
    cerr<<"Error opening reference genome file: "<<reffilename<<endl;
    return -1;
  }
  refseq_destroy();
  refseq=kseq_init(refFile);
  int l=0;
  bool foundchr=false;
  while((l=kseq_read(refseq))>=0){
    string seqname(refseq->name.s);
    int chrlen=strlen(refseq->seq.s);
    refgenome[seqname]=chrlen;
    cerr<<"Checking "<<seqname<<", length:"<<chrlen<<endl;
  }
  gzclose(refFile);
  return 0;
}

/* load the chromosome */
int refseq_loadchr(string chrname){
  refFile=gzopen(reffilename.c_str(),"r");
  if(refFile == NULL){
    cerr<<"Error opening reference genome file: "<<reffilename<<endl;
    return -1;
  }
  refseq_destroy();
  refseq=kseq_init(refFile);
  int l=0;
  bool foundchr=false;
  while((l=kseq_read(refseq))>=0){
    string seqname(refseq->name.s);
    int chrlen=strlen(refseq->seq.s);
    if(seqname ==  chrname) {
      cerr<<"Switch to "<<seqname<<", length:"<<chrlen<<endl;
      foundchr=true;
      break;
    }
  }
  gzclose(refFile);
  if(foundchr==true)
    currentchr=chrname;
  else{
    cerr<<"Error: sequence "<<chrname<<" not found.\n";
    currentchr="";
    return -1;
  }
  return 0;
}

/* Get the sequence content */
int refseq_getseq(string chrname, long a, int len, string& seq){
  if(refgenome.count(chrname)==0) return -1;
  if( currentchr != chrname){ 
    if(refseq_loadchr(chrname)!=0){
      return -1;
    }
    buff_a = buff_b =-1;
  }
  
  // update string buffer
  long b=a+len;
  if(a<=buff_a || a >= buff_b || b <= buff_a || b >= buff_b){
    //get the reference sequence 
    buff_a=a-buff_a_size;
    if (buff_a <0) buff_a=0;
    buff_b=b+buff_b_size;
    if (buff_b > strlen(refseq->seq.s)-1) buff_b =strlen(refseq->seq.s)-1; 
    
    string nstr(refseq->seq.s,buff_a,buff_b-buff_a+1);
    ref_buff=nstr;
  }
  seq=ref_buff.substr(a-buff_a,len);
  // to upper case
  for(int k=0;k<seq.size();k++) 
    if(seq[k] >='a' && seq[k] <='z') 
      seq[k]-=32;
  return 0;
}




