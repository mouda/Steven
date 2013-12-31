#include "clusterStructure.h"

ClusterStructure::ClusterStructure(const int numNodes, const int maxNumHeads):
  m_numNodes(numNodes), m_maxNumHeads(maxNumHeads)
{

}

ClusterStructure::~ClusterStructure()
{

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


