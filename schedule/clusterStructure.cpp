#include "clusterStructure.h"

ClusterStructure::ClusterStructure(const int numNodes, const int maxNumHeads):
  m_numNodes(numNodes), m_maxNumHeads(maxNumHeads)
{

}

ClusterStructure::~ClusterStructure()
{

}

void
ClusterStructure::SetRecord( 
    const vector<int>& inVecHeadName, 
    const list<list<int> >& inListCluMember
    )
{

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
    cout << "cluster: " << i <<"-th ";
    list<int>::const_iterator iterCols = iterRows->begin();
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


