#include "minPowerSA.h"

MinPowerSA::MinPowerSA( const double txTime, 
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("MinPowerSA")
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads =  m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();
  m_vecNodePower.resize(m_numNodes);
  fill(m_vecNodePower.begin(), m_vecNodePower.end(), m_maxPower);
  m_vecSched.resize(m_numNodes);
  fill(m_vecSched.begin(), m_vecSched.end(), 0);

  m_ptrSASolver = new SASolver(ptrMap, ptrMatComputer, ptrCS, m_txTimePerSlot);
}

MinPowerSA::~MinPowerSA()
{
  delete m_ptrSASolver;

}

double
MinPowerSA::ScheduleOneSlot( vector<int>& vecSupport )
{
  CheckAllScheduled(); 
  return 0.0;

}

double
MinPowerSA::ScheduleOneSlot( std::vector<int>& vecSupport, const std::vector<double>& vecVariance)
{
  m_ptrSASolver->Solve(vecSupport);
  return 0.0;
}

bool 
MinPowerSA::CheckFeasible( const vector<int>& supStru, double txTime2nd)
{
  for (int i = 0; i < m_numMaxHeads; i++) {
    int headName = m_ptrCS->GetVecHeadName()[i];
    int member = -1;
    double interference = 0.0;
    for (int j = 0; j < m_numNodes; j++) {
      if (supStru[j] == 1 && headName != m_ptrCS->GetChNameByName(j) ) {
        interference += m_ptrMap->GetGijByPair(headName,j) * m_maxPower;
      }
      else if(supStru[j] == 1 && headName == m_ptrCS->GetChNameByName(j) ){
        member = j;
      }
    }
    if (member == -1) continue; 
    if (m_ptrMap->GetIdtEntropy() > m_bandwidthKhz*txTime2nd*log2(1.0+m_maxPower * m_ptrMap->GetGijByPair(headName,member)
          / (m_ptrMap->GetNoise() + interference))) {
      return false;
    }
  }
  return true;
}

bool
MinPowerSA::CheckGroupScheduled( const int idx)
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
MinPowerSA::CheckAllScheduled()
{
  for (int i = 0; i < m_ptrMap->GetNumInitHeads(); ++i) {
    if (CheckGroupScheduled(i)) {
      ResetGroupScheduled(i);
    }
  }
  return true;
}

bool
MinPowerSA::ResetGroupScheduled( const int idx)
{
  assert(idx >= 0 && idx < m_ptrMap->GetNumInitHeads());
  for (int i = 0; i < m_vecSched.size(); ++i) {
    if ( m_ptrCS->GetChIdxByName(i) == idx ) {
      m_vecSched.at(i) = 0;
    }
  }
  return true;
}
