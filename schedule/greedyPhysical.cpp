#include "greedyPhysical.h"

GreedyPhysical::GreedyPhysical( const double txTime, 
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("GreedyPhysical"),
    m_commGraph(ptrMap->GetNumNodes()),
    m_conflictGraph(ptrMap->GetNumNodes())
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads =  m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();

  /* initialize the fixed power */
  m_vecNodePower.resize(m_numNodes);
  fill(m_vecNodePower.begin(), m_vecNodePower.end(), m_maxPower);

  /* initialize the scheduled record */
  m_vecSched.resize(m_numNodes);
  fill(m_vecSched.begin(), m_vecSched.end(), 0);


  /* construct the conflict graph */
  list<list<int> >::const_iterator  itRow = m_ptrCS->GetListCluMemeber().begin();
  for (int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
    if (m_ptrCS->GetChNameByName(i) == i) continue; /* head doesn't transmit */ 
    for (int j = 0; j < m_ptrCS->GetNumHeads(); ++j) {
      int headName = m_ptrCS->GetVecHeadName()[j];
      if (m_ptrCS->GetChIdxByName(i) != j) {
        property<edge_weight_t, double> tmpProperty(GetRxPower(i,headName)); 
        add_edge(
            vertex(i, m_conflictGraph),
            vertex(headName, m_conflictGraph),
            tmpProperty,
            m_conflictGraph
            );
      }
    }
  }

  /* construct the communicate graph */
  itRow = m_ptrCS->GetListCluMemeber().begin(); 
  for (; itRow != m_ptrCS->GetListCluMemeber().end(); ++itRow) {
    list<int>::const_iterator itCol = itRow->begin();
    for (; itCol != itRow->end(); ++itCol) {
      if (m_ptrCS->GetChNameByName(*itCol) != *itCol) {
        property<edge_weight_t, double> tmpProperty(GetRxPower(*itCol,m_ptrCS->GetChNameByName(*itCol)));
        add_edge(
            vertex(*itCol, m_commGraph),  /* memeber */
            vertex(m_ptrCS->GetChNameByName(*itCol), m_commGraph), /* cluster head */
            tmpProperty,
            m_commGraph
            );
      } 
    }
  }

#ifdef DEBUG
  EdgeIter ep, ep_end;
  BglVertex u, v;
  BglVertexMap vertexMap = get( vertex_index, m_conflictGraph);
  BglEdgeMap   edgeMap = get( edge_weight, m_conflictGraph);
  for (tie(ep,ep_end) = edges(m_conflictGraph); ep != ep_end; ++ep) {
    u = source(*ep, m_conflictGraph);
    v = target(*ep, m_conflictGraph);
    cout << "u: " << vertexMap[u] <<" v: " << vertexMap[v] << " weight " << edgeMap[*ep] <<endl;
  }
#endif

}

GreedyPhysical::~GreedyPhysical()
{

}


double 
GreedyPhysical::ScheduleOneSlot(std::vector<int>& vecSupport )
{
  return 0.0;
}

bool 
GreedyPhysical::ScheduleOneSlot(std::vector<int>& vecSupport, const std::vector<double>& vecVariance)
{
  CheckAllScheduled();
  std::vector<BglEdge> vecEdge;
  int selectedNode = -1;
  BglVertex u, v;
  while( vecEdge.size() < m_ptrCS->GetNumHeads() ) {
    selectedNode = GreedySelectOneNode(vecEdge);
    if (selectedNode != -1) {
      u = vertex(selectedNode, m_commGraph);
      v = vertex(m_ptrCS->GetChNameByName(selectedNode), m_commGraph);
      BglEdge myEdge = edge( u, v, m_commGraph).first;
      vecEdge.push_back(myEdge);
    }
    else {
      break;
    }
  }

  BglVertexMap vertexMap;
  for (int i = 0; i < vecEdge.size(); ++i) {
    u = source(vecEdge.at(i), m_commGraph);
    vecSupport.at(vertexMap[u]) = 1;
    m_vecSched.at(vertexMap[u]) = 1;
  }
  return true;
}

int
GreedyPhysical::GetInterferenceNumber( const int source, const int target)
{
  std::vector<BglEdge> nullVector;
  return GetInterferenceNumber( source, target, nullVector);
}

int
GreedyPhysical::GetInterferenceNumber( const int source, const int target, const std::vector<BglEdge>& vecEdge)
{
  /* we only need to consider the rx condition of cluster head */
  list<list<int> >::const_iterator itRow = m_ptrCS->GetListCluMemeber().begin();
  int counter = 0.0;
  BglVertex u, v;
  BglEdgeMap edgeMap = get( edge_weight, m_conflictGraph);
  for (; itRow != m_ptrCS->GetListCluMemeber().end(); ++itRow) {
    list<int>::const_iterator itCol = itRow->begin();
    if (m_ptrCS->GetChNameByName(*itCol) == target) continue; /* we don't need to the self cluster */ 
    for (; itCol != itRow->end(); ++itCol) {
      if (*itCol == m_ptrCS->GetChNameByName(*itCol)) continue; /* cannot be head */

      u = vertex(source, m_conflictGraph);  /* the interferencing node */
      v = vertex(m_ptrCS->GetChNameByName(*itCol), m_conflictGraph); /* the interferenced node */
      double rxPower =  edgeMap[edge(u,v,m_conflictGraph).first];

      for (int i = 0; i < vecEdge.size(); ++i) rxPower+=edgeMap[vecEdge.at(i)]; 

      int headName = m_ptrCS->GetChNameByName(*itCol);
      if (m_ptrMap->GetIdtEntropy() > 
          m_bandwidthKhz * m_txTimePerSlot*log2(1.0+ m_ptrMap->GetMaxPower() * m_ptrMap->GetGijByPair(*itCol,headName) 
            / (m_ptrMap->GetNoise() + rxPower))) {
        ++counter;
      }
    }
  }
  return counter;
} 

int 
GreedyPhysical::GreedySelectOneNode( const std::vector<BglEdge>& vecEdge)
{
  int maxInterferingNode = -1;
  int maxInterferingNum = 0;
  for (int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
    int headName = m_ptrCS->GetChNameByName(i);
    if (!ClusterSelected(vecEdge,i) && headName != i && m_vecSched.at(i) == 0 ) {
      int num = GetInterferenceNumber(i, headName, vecEdge); 
      if (num > maxInterferingNum) {
        maxInterferingNum = num;
        maxInterferingNode = i;
      }
    }
  }
  return maxInterferingNode;
}

bool
GreedyPhysical::ClusterSelected( const std::vector<BglEdge>& vecEdge, const int nodeName )
{
  BglVertex v;
  BglVertexMap vertexMap;
  for (int i = 0; i < vecEdge.size(); ++i) {
    v = target(vecEdge.at(i), m_commGraph);
    if (m_ptrCS->GetChNameByName(nodeName) == m_ptrCS->GetChNameByName(vertexMap[v]) ) {
      return true;
    }
  }
  return false;
}

/* rx power */
double 
GreedyPhysical::GetRxPower( const int source, const int target)
{
  return m_ptrMap->GetGijByPair(source, target) * m_maxPower;
}

bool
GreedyPhysical::CheckGroupScheduled( const int idx)
{
  assert(idx >= 0 && idx < m_ptrMap->GetNumInitHeads());
  for (int i = 0; i < m_numNodes; ++i) {
    if ( m_ptrCS->GetChNameByName(i) == i ) continue;
    if ( m_ptrCS->GetChIdxByName(i) == idx && m_vecSched.at(i) == 0 ) {
      return false;
    }
  }
  return true;
}

bool
GreedyPhysical::CheckAllScheduled()
{
  for (int i = 0; i < m_ptrMap->GetNumInitHeads(); ++i) {
    if (CheckGroupScheduled(i)) {
      ResetGroupScheduled(i);
    }
  }
  return true;
}

bool
GreedyPhysical::ResetGroupScheduled( const int idx)
{
  assert(idx >= 0 && idx < m_ptrMap->GetNumInitHeads());
  for (int i = 0; i < m_vecSched.size(); ++i) {
    if ( m_ptrCS->GetChIdxByName(i) == idx ) {
      m_vecSched.at(i) = 0;
    }
  }
  return true;
}
