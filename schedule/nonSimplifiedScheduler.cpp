#include "nonSimplifiedScheduler.h"

NonSimplifiedScheduler::NonSimplifiedScheduler( const double txTime, 
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("NonSimplifiedScheduler")
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
}

NonSimplifiedScheduler::~NonSimplifiedScheduler()
{

}

double 
NonSimplifiedScheduler::ScheduleOneSlot(std::vector<int>& vecSupport )
{
  return 0.0;
}

double
NonSimplifiedScheduler::ScheduleOneSlot(std::vector<int>& vecSupport, std::vector<double>& vecPower, const std::vector<double>& vecVariance)
{
  return 0.0;
}
