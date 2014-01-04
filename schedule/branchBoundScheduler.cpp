#include "branchBoundScheduler.h"

BranchBoundScheduler::BranchBoundScheduler( const double txTime, 
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE const * const ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("BranchBound")
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads = m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();

}

BranchBoundScheduler::~BranchBoundScheduler()
{

}

double
BranchBoundScheduler::ScheduleOneSlot( vector<int>& vecSupport )
{
  return 0.0;
}
