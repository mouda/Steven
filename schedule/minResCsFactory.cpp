#include "minResCsFactory.h"

MinResCsFactory::MinResCsFactory( Map const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  CsFactory(myMap, myMatComputer) 
{

}

MinResCsFactory::~MinResCsFactory()
{

}

ClusterStructure*
MinResCsFactory::CreateClusterStructure()
{
  if (m_ptrCS == 0) {
    vector<int> myVecHeadNames;
    list<list<int> > myListCluMembers;
    if (SASearch(myVecHeadNames, myListCluMembers)) {
      m_ptrCS = new ClusterStructure(m_ptrMap->GetNumNodes(), 
          m_ptrMap->GetNumInitHeads() );
      m_ptrCS->SetRecord(myVecHeadNames, myListCluMembers);
      return m_ptrCS;
    }
    else{
      return NULL;
    }
  }
}

bool
MinResCsFactory::SASearch( vector<int>& vecHeadNames, list<list<int> >& listCluMembers )
{
  return true;
}
