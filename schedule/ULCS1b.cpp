/*
  File: ULCS1b.h
  Brief: Record all the coalition structure information of SA
        1b: keep one more new structure
        allSupStru[totalNodes]

  Author: Steven
  Date: 2012/09/27.
  Latest Progress:

  Abbreviation:

  Quote:
  Index means the array index, so the data store in the vecHeadName is head name. However we need index to access it.
  Name means the number of each machine

*/
#include "ULCS1b.h"
#include<cstdio>
#include<iostream>
#include<cassert>


ULCS1b::ULCS1b() {}
ULCS1b::ULCS1b(int inputTotalNodes, int inputMaxChNum)
{
  systemon = false;
  maxChNum = inputMaxChNum;
  totalNodes = inputTotalNodes;
  clusterStru = new bool* [maxChNum];
  for (int i =0; i<maxChNum ; i++) clusterStru[i] = new bool [totalNodes];
  for(int i=0; i<maxChNum; i++)
  {
    for(int j=0; j<totalNodes; j++)
    {
      clusterStru[i][j] = false;
    }
  }
  allSupStru = new bool [totalNodes];
  for (int i=0; i<totalNodes; i++)
    allSupStru[i]=false;
  listCluMember = new list<list<int> > ;
  listUnSupport = new list<int>;
}
ULCS1b::~ULCS1b()
{
  for (int i =0; i<maxChNum; i++) delete [] clusterStru[i];
  delete [] clusterStru;
  delete [] allSupStru;
  delete listCluMember;
  delete listUnSupport;
}
void ULCS1b::addNewHeadCs(int inputHeadName)
{
  vector <bool> temp;
  int tempHeadIndex = -1;

  vecHeadName.push_back(inputHeadName);
  tempHeadIndex = vecHeadName.size()-1;
  clusterStru[tempHeadIndex][vecHeadName[tempHeadIndex]]=true;
  allSupStru[vecHeadName[tempHeadIndex]]=true;
  vecClusterSize.push_back(0);//Size of each Cluster = 0 because it hasn't known it's real head and member yet
  (*listCluMember).resize((*listCluMember).size()+1);
}

/*
    add new member in cluster system
    called by ULSA3b::addMemberSA_MachineCentric(int headIndex, int memberName)
*/
void ULCS1b::addMemberCs(int headIndex, int memberName,bool iniDone)
{
  clusterStru[headIndex][memberName] = true; //Please notice "memberName = memberIndex" here.
  allSupStru[memberName]=true;
  vecClusterSize[headIndex]++;
  list <list<int> >::iterator it = (*listCluMember).begin();
  for(int i=0; i<headIndex; i++) it++;
  (*it).push_back(memberName);
  //Only after the intial structure is constructed, we start to manipulate the listUnSupport

  //*****
  //the erase job needed to be done in decide , which will speed up the program

  if(iniDone)
  {
    list<int> ::iterator it2 =listUnSupport->begin();

    bool checkFoundUnsupportAdd =false;
    for(; it2!=listUnSupport->end(); it2++)
    {
      if(memberName==(*it2))
      {
        it2 =listUnSupport->erase(it2);
        checkFoundUnsupportAdd = true;
        break;
      }

    }
    if(!checkFoundUnsupportAdd){
        cout<<"Cannot find desired node:"<<memberName<<endl;
        it2 =listUnSupport->begin();
        cout<<"Unsupportset=";
        for(; it2!=listUnSupport->end(); it2++)
        {
            cout<<*it2<<", ";
        }
        cout<<endl;
        assert(checkFoundUnsupportAdd);
    }
  }
}
/*
     discard certain member form it's head
*/
void ULCS1b::discardMemberCs(int headIndex, int memberName )//headIndex := headName
{
  clusterStru[headIndex][memberName] = false; //Please note that accoring to my design memberName = memberIndex
  allSupStru[memberName]=false;
  vecClusterSize[headIndex]--;
  listUnSupport->push_back(memberName);
  list <list <int> > ::iterator it1=listCluMember->begin();
  for(int i=0; i<headIndex; i++)it1++;
  iteLastDiscard = it1 ;// record the last discarded iterator
  list <int>::iterator it2= it1->begin();
  bool veriFlag = false;
  for(; it2!=it1->end(); it2++)
  {
    if(memberName == (*it2))
    {
      it2 = it1->erase(it2);//the return value(it2) means the iterator after original it2.However, we don't need that here.
      veriFlag = true;
      break;
    }
  }
  it1=listCluMember->begin();
  for(int i=0; i<headIndex; i++)it1++;
  assert(it1==iteLastDiscard);

  assert(veriFlag);
}

void ULCS1b::join_FromHead(int JoiningHeadIndex,int targetHeadIndex){
     list<list <int> >::iterator it_LiInt=listCluMember->begin();
     list<list <int> >::iterator it_LiInt2=listCluMember->begin();

     for(int i=0;i<JoiningHeadIndex;i++)it_LiInt++;
     for(int i=0;i<targetHeadIndex;i++)it_LiInt2++;

     //add to another
     list<int>::iterator it_Int=it_LiInt->begin();
     for(unsigned int i=0;i<it_LiInt->size();i++,it_Int++){

        it_LiInt2->push_back(*it_Int);
        clusterStru[targetHeadIndex][*it_Int]=1;
     }
     vecClusterSize[targetHeadIndex]+=vecClusterSize[JoiningHeadIndex];


     it_Int=it_LiInt->begin();
     *it_Int=-1;
     int OriSiz=it_LiInt->size();
     for(unsigned int i=0;i<(OriSiz-1);i++)
         it_LiInt->pop_back();


     assert(it_LiInt->size()==1);
     for (int i=0;i< totalNodes;i++)
        clusterStru[JoiningHeadIndex][i]=0;
     vecHeadName[JoiningHeadIndex]=-1;
     vecClusterSize[JoiningHeadIndex]=0;
}

void ULCS1b::isolate_FromHead(int isoNodeName,int isolatedCluI,int targetH){
    //add;
    clusterStru[targetH][isoNodeName]=1;
    vecHeadName[targetH]=isoNodeName;
    list<list <int> >::iterator it_LiInt=listCluMember->begin();
    for(int i=0;i<targetH;i++)it_LiInt++;
    list<int>::iterator it_Int=it_LiInt->begin();
    *it_Int=isoNodeName;
    vecClusterSize[targetH]=1;

    //delete;
    clusterStru[isolatedCluI][isoNodeName]=0;
    it_LiInt=listCluMember->begin();
    for(int i=0;i<isolatedCluI;i++)it_LiInt++;
    bool found=false;
    it_Int=it_LiInt->begin();
    for(unsigned int i=0;i<it_LiInt->size();i++,it_Int++){
        if(*it_Int==isoNodeName)it_Int=it_LiInt->erase(it_Int);
    }
    vecClusterSize[isolatedCluI]--;
}
/*
    After "add" to the neighbor we and we decide we are not going to in that neighbor
*/
void ULCS1b::reverseAdd(int targetHeadIndex,int targetMember)
{
  //delete the just added node in the cSystem
  list <list <int> > ::iterator it1 = listCluMember->begin();
  for(int i=0; i<targetHeadIndex; i++)it1++;
  it1->pop_back();
  listUnSupport->push_back(targetMember);
  clusterStru[targetHeadIndex][targetMember] = false;
  allSupStru[targetMember]=false;
  vecClusterSize[targetHeadIndex]--;
}

/*
  After "Discard", but we want to reverse in the end

*/
void ULCS1b::reverseDiscard(int targetHeadIndex, int targetMember)
{
  if (iteLastDiscard != listCluMember->end())
    iteLastDiscard->push_back(targetMember);
  else
  {
    cout<<"error iteLastDiscard it shouldn't be NULL here"<<endl;
    assert(1);
  }
  listUnSupport->pop_back();
  clusterStru[targetHeadIndex][targetMember] = true;
  allSupStru[targetMember]=true;
  vecClusterSize[targetHeadIndex]++;
  iteLastDiscard = listCluMember->end();
}
//
//  @Clear all the existed cluster system structure
//

void ULCS1b::reverseJoin(int lastJoinHeadIndex,int lastJoinHead, int lastTargetIndex,vector<int>&lastjoinall){
    //add back
    vecHeadName[lastJoinHeadIndex]=lastJoinHead;
    vecClusterSize[lastJoinHeadIndex]=lastjoinall.size();
    for(unsigned int i=0;i<lastjoinall.size();i++){
        clusterStru[lastJoinHeadIndex][lastjoinall[i]]=1;
    }
    list<list <int> >::iterator it_LiInt=listCluMember->begin();

    for(int i=0;i<lastJoinHeadIndex;i++)it_LiInt++;
    list<int>::iterator it_Int=it_LiInt->begin();

    for(unsigned int i=0;i<lastjoinall.size();i++){
        if(i==0)*it_Int=lastjoinall[i];
        else
            it_LiInt->push_back(lastjoinall[i]);
    }

    //delete joined
    it_LiInt=listCluMember->begin();
    for(int i=0;i<lastTargetIndex;i++)it_LiInt++;

    for(unsigned int i=0;i<lastjoinall.size();i++){
        it_LiInt->pop_back();
        clusterStru[lastTargetIndex][lastjoinall[i]]=0;
    }
    vecClusterSize[lastTargetIndex]-=lastjoinall.size();
}
void ULCS1b::reverseisolate(int lastisoNodName,int lastisoCluI,int lastisoTargetH){
    //add back
    vecClusterSize[lastisoCluI]++;
    list<list <int> >::iterator it_LiInt=listCluMember->begin();
    for(int i=0;i<lastisoCluI;i++)it_LiInt++;
    it_LiInt->push_back(lastisoNodName);
    clusterStru[lastisoCluI][lastisoNodName]=1;

    //delete added
    vecClusterSize[lastisoTargetH]--;
    it_LiInt=listCluMember->begin();
    for(int i=0;i<lastisoTargetH;i++)it_LiInt++;
    list<int>::iterator it_Int=it_LiInt->begin();
    *it_Int=-1;
    vecHeadName[lastisoTargetH]=-1;
    clusterStru[lastisoTargetH][lastisoNodName]=0;




}

void ULCS1b::resetSystem()
{
  vecHeadName.clear();
  vecClusterSize.clear();
  vecHeadName.reserve(maxChNum);
  vecClusterSize.reserve(maxChNum);
  listCluMember->clear();
  listUnSupport->clear();
  for(int i=0;i<maxChNum;i++)
  {
    allSupStru[i]=false;

    for(int j=0;j<totalNodes;j++){
      clusterStru[i][j]=false;
    }
  }
}
bool ULCS1b::returnIfClusterSmall(int threshold, int &numOfClu){
    numOfClu=0;
    bool tmpFlag=false;
    for(int i=0;i<maxChNum;i++){
        if(vecClusterSize[i]>0&&vecClusterSize[i]<=threshold){
            tmpFlag=true;
            numOfClu++;
        }
    }
    return tmpFlag;
}

int ULCS1b::getCHIdxByCHName(const int& name){
  assert( vecHeadName.size() != 0 );
  int headIdx = 0;
  for ( headIdx = 0; headIdx < vecHeadName.size(); headIdx++ ) {
    if( vecHeadName[headIdx] == name ) break; 
  }
  return headIdx; 
}





//--------------------------------------------//
//Internal Tool                               //
//--------------------------------------------//

/*
    Initialization ONlY: return head iterator of certain headIndex
*/
int* ULCS1b::returnHeadPtr(int inputHeadIndex)
{
  int* ptrtemp = NULL;
  vector <int>::iterator  it;
  it = vecHeadName.begin();
  for(int i=0; i<inputHeadIndex; i++)it++;
  ptrtemp = &(*it);
  return ptrtemp;
}


/*
   Compute the number of support's node in the cSystem->listCluMember
   This is the payoff of Version ULSA2x
*/
int ULCS1b::calSupNodes()
{
  list <list<int> >:: iterator it1 = listCluMember->begin();
  int supNodes =0;
  for(; it1!=listCluMember->end(); it1++)supNodes+=it1->size();
  return supNodes;
}






















/*
    Overload a new operator of "=" for class ULCS1b
*/
/*ULCS1b& ULCS1b::operator=(const ULCS1b& rhs)
{
    totalNodes = rhs.totalNodes;
    maxChNum = rhs.maxChNum;
    systemon = rhs.systemon;
    vecHeadName = rhs.vecHeadName;
    for (int i=0;i<maxChNum;i++)
        for(int j=0;j<totalNodes;j++)
            clusterStru[i][j] = rhs.clusterStru[i][j];
    return (*this);
}*/
