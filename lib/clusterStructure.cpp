#include <cassert>
#include "clusterStructure.h"

ClusterStructure::ClusterStructure(const int numNodes, const int maxNumHeads):
  m_numNodes(numNodes), m_maxNumHeads(maxNumHeads), m_vecSupport(numNodes),
  m_vecPower(numNodes), m_matClusterStru(maxNumHeads,vector<int>(numNodes)),
  m_allSupStru(numNodes)
{
  fill(m_vecSupport.begin(), m_vecSupport.end(), 0);
  fill(m_vecPower.begin(), m_vecPower.end(), 0.0);
  fill(m_allSupStru.begin(), m_allSupStru.end(), 1);

  for(int i=0; i<m_maxNumHeads; i++) {
    for(int j=0; j<m_numNodes; j++) {
      m_matClusterStru[i][j] = 0;
    }
  }
}

ClusterStructure::~ClusterStructure()
{

}

void
ClusterStructure::SetRecord( 
    const vector<int>& inVecHeadName, 
    const list<list<int> >& inListCluMember,
    const vector<int>& vecAllSupport
    )
{
  m_allSupStru.assign(vecAllSupport.begin(), vecAllSupport.end());
  m_vecCHIdxForNodes.resize(m_numNodes);
  std::fill(m_vecCHIdxForNodes.begin(),m_vecCHIdxForNodes.end(),-1);
  m_vecCHNameForNodes.resize(m_numNodes);
  std::fill(m_vecCHNameForNodes.begin(),m_vecCHNameForNodes.end(),-1);
  m_vecHeadName = inVecHeadName;
  m_listCluMember = inListCluMember;
  list<list<int> >::const_iterator iterRow = m_listCluMember.begin();
  for (int i = 0; iterRow != m_listCluMember.end(); ++i, ++iterRow) {
    list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      m_vecCHIdxForNodes[*iterCol] = i; 
      m_vecCHNameForNodes[*iterCol] = m_vecHeadName[i];
    }
  }
}

void
ClusterStructure::Print() const
{
  list<list<int> >::const_iterator iterRows = m_listCluMember.begin();
  for (int i = 0;iterRows != m_listCluMember.end(); ++iterRows, ++i) {
    cout << "cluster: " << i <<"-th, Head:" ;
    list<int>::const_iterator iterCols = iterRows->begin();
    cout << m_vecCHNameForNodes.at(*iterCols) << ", ";
    for (;iterCols != iterRows->end(); ++iterCols) {
      cout << *iterCols << ' ';
    }
    cout << endl;
  }
}

int 
ClusterStructure::GetChIdxByName( const int& nodeName ) const 
{
  return m_vecCHIdxForNodes[nodeName];
}

int 
ClusterStructure::GetChNameByName( const int& nodeName) const
{
  return m_vecCHNameForNodes[nodeName];
}


void 
ClusterStructure::addNewHeadCs(int inputHeadName)
{
  vector <bool> temp;
  int tempHeadIndex = -1;

  m_vecHeadName.push_back(inputHeadName);
  tempHeadIndex = m_vecHeadName.size()-1;
  m_matClusterStru[tempHeadIndex][m_vecHeadName[tempHeadIndex]]=1;
  m_allSupStru[m_vecHeadName[tempHeadIndex]] = 1;
  m_vecClusterSize.push_back(0);//Size of each Cluster = 0 because it hasn't known it's real head and member yet
  m_listCluMember.resize(m_listCluMember.size()+1);
}

/*
    add new member in cluster system
    called by ULSA3b::addMemberSA_MachineCentric(int headIndex, int memberName)
*/
void 
ClusterStructure::addMemberCs(int headIndex, int memberName,bool iniDone)
{
  m_matClusterStru[headIndex][memberName] = 1; //Please notice "memberName = memberIndex" here.
  m_allSupStru[memberName] = 1;
  m_vecClusterSize[headIndex]++;
  list <list<int> >::iterator it = m_listCluMember.begin();
  for(int i=0; i<headIndex; i++) it++;
  (*it).push_back(memberName);
  //Only after the intial structure is constructed, we start to manipulate the listUnSupport

  //*****
  //the erase job needed to be done in decide , which will speed up the program

  if(iniDone)
  {
    list<int> ::iterator it2 =m_listUnSupport.begin();

    bool checkFoundUnsupportAdd =false;
    for(; it2!=m_listUnSupport.end(); it2++)
    {
      if(memberName==(*it2))
      {
        it2 =m_listUnSupport.erase(it2);
        checkFoundUnsupportAdd = true;
        break;
      }

    }
    if(!checkFoundUnsupportAdd){
        cout<<"Cannot find desired node:"<<memberName<<endl;
        it2 =m_listUnSupport.begin();
        cout<<"Unsupportset=";
        for(; it2!=m_listUnSupport.end(); it2++)
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
void 
ClusterStructure::discardMemberCs(int headIndex, int memberName )//headIndex := headName
{
  m_matClusterStru[headIndex][memberName] = 0; //Please note that accoring to my design memberName = memberIndex
  m_allSupStru[memberName] = 0;
  m_vecClusterSize[headIndex]--;
  m_listUnSupport.push_back(memberName);
  list <list <int> > ::iterator it1=m_listCluMember.begin();
  for(int i=0; i<headIndex; i++)it1++;
  m_iteLastDiscard = it1 ;// record the last discarded iterator
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
  it1=m_listCluMember.begin();
  for(int i=0; i<headIndex; i++)it1++;
  assert(it1==m_iteLastDiscard);

  assert(veriFlag);
}

void 
ClusterStructure::join_FromHead(int JoiningHeadIndex,int targetHeadIndex){
     list<list <int> >::iterator it_LiInt=m_listCluMember.begin();
     list<list <int> >::iterator it_LiInt2=m_listCluMember.begin();

     for(int i=0;i<JoiningHeadIndex;i++)it_LiInt++;
     for(int i=0;i<targetHeadIndex;i++)it_LiInt2++;

     //add to another
     list<int>::iterator it_Int=it_LiInt->begin();
     for(unsigned int i=0;i<it_LiInt->size();i++,it_Int++){

        it_LiInt2->push_back(*it_Int);
        m_matClusterStru[targetHeadIndex][*it_Int]=1;
     }
     m_vecClusterSize[targetHeadIndex]+=m_vecClusterSize[JoiningHeadIndex];


     it_Int=it_LiInt->begin();
     *it_Int=-1;
     int OriSiz=it_LiInt->size();
     for(unsigned int i=0;i<(OriSiz-1);i++)
         it_LiInt->pop_back();


     assert(it_LiInt->size()==1);
     for (int i=0;i< m_numNodes;i++)
        m_matClusterStru[JoiningHeadIndex][i]=0;
     m_vecHeadName[JoiningHeadIndex]=-1;
     m_vecClusterSize[JoiningHeadIndex]=0;
}

void 
ClusterStructure::isolate_FromHead(int isoNodeName,int isolatedCluI,int targetH){
    //add;
    m_matClusterStru[targetH][isoNodeName]=1;
    m_vecHeadName[targetH]=isoNodeName;
    list<list <int> >::iterator it_LiInt=m_listCluMember.begin();
    for(int i=0;i<targetH;i++)it_LiInt++;
    list<int>::iterator it_Int=it_LiInt->begin();
    *it_Int=isoNodeName;
    m_vecClusterSize[targetH]=1;

    //delete;
    m_matClusterStru[isolatedCluI][isoNodeName]=0;
    it_LiInt=m_listCluMember.begin();
    for(int i=0;i<isolatedCluI;i++)it_LiInt++;
    bool found=false;
    it_Int=it_LiInt->begin();
    for(unsigned int i=0;i<it_LiInt->size();i++,it_Int++){
        if(*it_Int==isoNodeName)it_Int=it_LiInt->erase(it_Int);
    }
    m_vecClusterSize[isolatedCluI]--;
}

/*
    After "add" to the neighbor we and we decide we are not going to in that neighbor
*/
void 
ClusterStructure::reverseAdd(int targetHeadIndex,int targetMember)
{
  //delete the just added node in the cSystem
  list <list <int> > ::iterator it1 = m_listCluMember.begin();
  for(int i=0; i<targetHeadIndex; i++)it1++;
  it1->pop_back();
  m_listUnSupport.push_back(targetMember);
  m_matClusterStru[targetHeadIndex][targetMember] = 0;
  m_allSupStru[targetMember] = 0;
  m_vecClusterSize[targetHeadIndex]--;
}

/*
  After "Discard", but we want to reverse in the end

*/
void 
ClusterStructure::reverseDiscard(int targetHeadIndex, int targetMember)
{
  if (m_iteLastDiscard != m_listCluMember.end())
    m_iteLastDiscard->push_back(targetMember);
  else
  {
    cout<<"error m_iteLastDiscard it shouldn't be NULL here"<<endl;
    assert(1);
  }
  m_listUnSupport.pop_back();
  m_matClusterStru[targetHeadIndex][targetMember] = 1;
  m_allSupStru[targetMember] = 1;
  m_vecClusterSize[targetHeadIndex]++;
  m_iteLastDiscard = m_listCluMember.end();
}
//
//  @Clear all the existed cluster system structure
//

void 
ClusterStructure::reverseJoin(int lastJoinHeadIndex,int lastJoinHead, int lastTargetIndex,vector<int>&lastjoinall){
    //add back
    m_vecHeadName[lastJoinHeadIndex]=lastJoinHead;
    m_vecClusterSize[lastJoinHeadIndex]=lastjoinall.size();
    for(unsigned int i=0;i<lastjoinall.size();i++){
        m_matClusterStru[lastJoinHeadIndex][lastjoinall[i]]=1;
    }
    list<list <int> >::iterator it_LiInt=m_listCluMember.begin();

    for(int i=0;i<lastJoinHeadIndex;i++)it_LiInt++;
    list<int>::iterator it_Int=it_LiInt->begin();

    for(unsigned int i=0;i<lastjoinall.size();i++){
        if(i==0)*it_Int=lastjoinall[i];
        else
            it_LiInt->push_back(lastjoinall[i]);
    }

    //delete joined
    it_LiInt=m_listCluMember.begin();
    for(int i=0;i<lastTargetIndex;i++)it_LiInt++;

    for(unsigned int i=0;i<lastjoinall.size();i++){
        it_LiInt->pop_back();
        m_matClusterStru[lastTargetIndex][lastjoinall[i]]=0;
    }
    m_vecClusterSize[lastTargetIndex]-=lastjoinall.size();
}

void 
ClusterStructure::reverseisolate(int lastisoNodName,int lastisoCluI,int lastisoTargetH){
    //add back
    m_vecClusterSize[lastisoCluI]++;
    list<list <int> >::iterator it_LiInt=m_listCluMember.begin();
    for(int i=0;i<lastisoCluI;i++)it_LiInt++;
    it_LiInt->push_back(lastisoNodName);
    m_matClusterStru[lastisoCluI][lastisoNodName]=1;

    //delete added
    m_vecClusterSize[lastisoTargetH]--;
    it_LiInt=m_listCluMember.begin();
    for(int i=0;i<lastisoTargetH;i++)it_LiInt++;
    list<int>::iterator it_Int=it_LiInt->begin();
    *it_Int=-1;
    m_vecHeadName[lastisoTargetH]=-1;
    m_matClusterStru[lastisoTargetH][lastisoNodName]=0;




}

void 
ClusterStructure::resetSystem()
{
  m_vecHeadName.clear();
  m_vecClusterSize.clear();
  m_vecHeadName.reserve(m_maxNumHeads);
  m_vecClusterSize.reserve(m_maxNumHeads);
  m_listCluMember.clear();
  m_listUnSupport.clear();
  for(int i=0;i<m_maxNumHeads;i++)
  {
    m_allSupStru[i] = 0;

    for(int j=0;j<m_numNodes;j++){
      m_matClusterStru[i][j] = 0;
    }
  }
}

bool 
ClusterStructure::returnIfClusterSmall(int threshold, int &numOfClu){
    numOfClu=0;
    bool tmpFlag=false;
    for(int i=0;i<m_maxNumHeads;i++){
        if(m_vecClusterSize[i]>0&&m_vecClusterSize[i]<=threshold){
            tmpFlag=true;
            numOfClu++;
        }
    }
    return tmpFlag;
}

/*
   Compute the number of support's node in the cSystem->listCluMember
   This is the payoff of Version ULSA2x
*/
int 
ClusterStructure::calSupNodes()
{
  list <list<int> >:: iterator it1 = m_listCluMember.begin();
  int supNodes =0;
  for(; it1 != m_listCluMember.end(); it1++)supNodes+=it1->size();
  return supNodes;
}

/*
    Initialization ONlY: return head iterator of certain headIndex
*/
int* 
ClusterStructure::returnHeadPtr(int inputHeadIndex)
{
  int* ptrtemp = NULL;
  vector <int>::iterator  it;
  it = m_vecHeadName.begin();
  for(int i=0; i<inputHeadIndex; ++i) ++it;
  ptrtemp = &(*it);
  return ptrtemp;
}

