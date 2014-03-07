/* store, search, and modify a set of mutations */

#include "mutmap.h"

/* Find the record pos:ref->alt, and return the corresponding index in vector<MutInfo>.
If no such record is found:
  If add=T: create a new MutInfo, set up the ref and alt (but not any else), and return this index;
  If add=F: return -1
*/
int findRecord(vector<MutInfo>& cvm, string ref, string alt, bool add){
  int recid=-1;
  for(int i=0;i<cvm.size();i++){
    if(cvm[i].ref==ref && cvm[i].alt==alt){
      recid=i;
      break;
    }
  }
  if(recid==-1){
    if(add){
      MutInfo mti;
      mti.ref=ref; mti.alt=alt;
      cvm.push_back(mti);
      recid=cvm.size()-1;
    }
  }
  return recid;
}
 
/* Add one mutation specified by pos:ref->alt */
int MutMap::addOneMut(long pos, string ref, string alt, bool forward,bool createifnoexist,int increment){
  
  vector<MutInfo>& cvm = mvm[pos];
  int recid=findRecord(cvm, ref, alt, createifnoexist);
  if(createifnoexist==false && recid==-1) return -1;
  if(forward){
    cvm[recid].altF+=increment;
  }else{
    cvm[recid].altB+=increment;
  }
  //cout<<"Adding "<<pos<<","<<cvm[recid].altF<<":"<<cvm[recid].altB<<endl;
  return 0;
}

/* remove all records and REFREC records before pos */
int MutMap::removeMutBeforePos(long pos){
  for(MAPVMUT_ITR itr=mvm.begin();itr!=mvm.end();){
    if(itr->first < pos){
      mvm.erase(itr++);
    }else{
      break;
    }
  }
  for(list<REFREC>::iterator lit=allrefrec.begin();lit!=allrefrec.end();){
    if(lit->end_pos<pos){
      allrefrec.erase(lit++);
    }else{
      ++lit;
    }
  }
  return 0;
}

int MutMap::addOneRefAlign(REFREC& rfrc){
  allrefrec.push_back(rfrc);
  return 0;
}

/* add the ref record for all mutations, except those specified in the blackout */ 
int MutMap::addOneRefAlign(long read0, vector<CigarOp>& cgo,bool forward, set<long>& blackout){
  for(MAPVMUT_ITR itr=mvm.begin();itr!=mvm.end();itr++){
    int bccnt=blackout.count(itr->first);
    if(bccnt>0) continue;
    int pir=posInRead(read0,cgo, itr->first);
    //cout<<"-- "<<itr->first<<" in read beginning with "<<read0<<": "<<pir<<endl;
    if(pir!=-1){
      vector<MutInfo>& cvm=itr->second;
      for(int i=0;i<cvm.size();i++){
        if(forward){
          cvm[i].refF++;
        }else{
          cvm[i].refB++;
        }
      }
    }//end if
  }//end for
  return 0;
}

/* update the REF record before the pos */
int MutMap::updateRefRecord(long pos){

  for(MAPVMUT_ITR itr=mvm.begin();itr!=mvm.end();itr++){
    if(itr->first >= pos) break;
    //update ref record
    for(list<REFREC>::iterator lit=allrefrec.begin();lit!=allrefrec.end();lit++){
      REFREC& rr=*lit;
      // addOneRefAlign(rr.start_pos, rr.CigarData,rr.direction,rr.blackout);
      int bccnt=rr.blackout.count(itr->first);
      if(bccnt>0) continue;
      int pir=posInRead(rr.start_pos,rr.CigarData, itr->first);
      //cout<<"-- "<<itr->first<<" in read beginning with "<<rr.start_pos<<": "<<pir<<", allrefrec.size()="<<allrefrec.size()<<endl;
      if(pir==-1) continue;
      vector<MutInfo>& cvm=itr->second;
      for(int i=0;i<cvm.size();i++){
        if(rr.direction){
          cvm[i].refF++;
        }else{
          cvm[i].refB++;
        }
      }// end for 
    }//end looping for all references
  }
  return 0;
}



