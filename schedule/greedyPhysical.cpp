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

  /* construct the communicate graph */
  list<list<int> >::const_iterator  itRow = m_ptrCS->GetListCluMemeber().begin(); 
  for (; itRow != m_ptrCS->GetListCluMemeber().end(); ++itRow) {
    list<int>::const_iterator itCol = itRow->begin();
    for (; itCol != itRow->end(); ++itCol) {
      if (m_ptrCS->GetChNameByName(*itCol) != *itCol) {
        add_edge(
            vertex(*itCol, m_commGraph),  /* memeber */
            vertex(m_ptrCS->GetChNameByName(*itCol), m_commGraph), /* cluster head */
            m_commGraph
            );
      } 
    }
  }

  /* construct the conflict graph */
  itRow = m_ptrCS->GetListCluMemeber().begin();
  for (int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
    if (m_ptrCS->GetChNameByName(i) == i) continue; /* head doesn't transmit */ 
    for (int j = 0; j < m_ptrMap->GetNumNodes(); ++j) {
      if (m_ptrCS->GetChNameByName(j) == j) continue; /* head doesn't transmit */ 
      if (m_ptrCS->GetChNameByName(i) != m_ptrCS->GetChNameByName(j)) {
        property<edge_weight_t, double> tmpProperty(GetConflictEdgeWeight(i,j)); 
        add_edge(
            vertex(i, m_conflictGraph),
            vertex(j, m_conflictGraph),
            tmpProperty,
            m_conflictGraph
            );
      }
    }
  }

  EdgeIter ep, ep_end;
  BglVertex u, v;
  BglVertexMap vertexMap = get( vertex_index, m_conflictGraph);
  BglEdgeMap   edgeMap = get( edge_weight, m_conflictGraph);
  for (tie(ep,ep_end) = edges(m_conflictGraph); ep != ep_end; ++ep) {
    u = source(*ep, m_conflictGraph);
    v = target(*ep, m_conflictGraph);
    cout << "u: " << vertexMap[u] <<" v: " << vertexMap[v] << " weight " << edgeMap[*ep] <<endl;
  }


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
  return true;
}

double
GreedyPhysical::GetInterferenceNumber( const int source, const int target)
{
  /* we only need to consider the rx condition of cluster head */
  list<list<int> >::const_iterator itRow = m_ptrCS->GetListCluMemeber().begin();
  double counter = 0.0;
  BglVertex u, v;
  for (; itRow != m_ptrCS->GetListCluMemeber().end(); ++itRow) {
    list<int>::const_iterator itCol = itRow->begin();
    if (m_ptrCS->GetChNameByName(*itCol) == target) continue; /* we don't need to the self cluster */ 
    for (; itCol != itRow->end(); ++itCol) {
      u = vertex(source, m_conflictGraph);
      v = vertex(*itCol, m_conflictGraph);
    }
  }

  return 0;
}

/* rx power */
double 
GreedyPhysical::GetConflictEdgeWeight( const int source, const int target)
{
  return m_ptrMap->GetGijByPair(source, target) * m_maxPower;
}
