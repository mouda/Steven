#include "ClusterStructure.h"
#include <algorithm>


ClusterStructure::~ClusterStructure()
{

}

void
ClusterStructure::SetRecord( 
    const vector<int>& inVecHeadName, 
    const list<list<int> >& inListCluMember
    )
{

  m_vecCHIdxForNodes.resize(m_numNode);
  std::fill(m_vecCHIdxForNodes.begin(),m_vecCHIdxForNodes.end(),-1);
  m_vecCHNameForNodes.resize(m_numNode);
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

int 
ClusterStructure::GetChIdxByName( const int& nodeName )
{
  return m_vecCHIdxForNodes[nodeName];
}

int 
ClusterStructure::GetChNameByName( const int& nodeName)
{
  return m_vecCHNameForNodes[nodeName];
}


