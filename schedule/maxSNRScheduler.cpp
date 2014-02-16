#include "maxSNRScheduler.h"

MaxSNRScheduler::MaxSNRScheduler( const double txTime, 
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("MaxSNR")
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads =  m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();
  m_vecNodePower.resize(m_numNodes);
  fill(m_vecNodePower.begin(), m_vecNodePower.end(), m_maxPower);
  m_vecSched.resize(m_numNodes);
  fill(m_vecSched.begin(), m_vecSched.end(), 0);

}

MaxSNRScheduler::~MaxSNRScheduler()
{

}

double
MaxSNRScheduler::ScheduleOneSlot( vector<int>& vecSupport )
{
  for (int i = 1; i < m_numMaxHeads; i++) {
    double maxRxPower = 0.0;
    int headName = m_ptrCS->GetVecHeadName()[i];
    int maxRxPowerNode = -1;
    for (int j = 0; j < m_numNodes; j++) {
      if (m_ptrCS->GetChNameByName(j) == j ) continue; 
      if (m_ptrCS->GetChIdxByName(j) == i && m_vecSched[j] == 0 ) {
        if (m_ptrMap->GetGijByPair(headName,j) * m_maxPower > maxRxPower) {
          maxRxPower = m_ptrMap->GetGijByPair(headName,j) * m_maxPower; 
          maxRxPowerNode = j;
        }
      }
    }
    if (maxRxPowerNode >= 0) {
      m_vecSched.at(maxRxPowerNode) = 1;
      vecSupport.at(maxRxPowerNode) = 1;
    }
  }

  /* test the Interference */
  while(!CheckFeasible(vecSupport, m_txTimePerSlot)){
    for (int i = 0; i < m_numNodes; i++) {
      if (vecSupport[i] == 1) {
        vecSupport[i] = 0;
        m_vecSched[i] = 0;
        break;
      }
    }
  }
  int activeNodes = 0;
  for (int i = 0; i < m_numNodes; i++) {
    if (vecSupport[i] == 1) ++activeNodes;
  }

  double result = activeNodes * m_ptrMap->GetIdtEntropy() + m_ptrMatComputer->computeLog2Det(1.0, vecSupport);
#ifdef DEBUG
  cout << "activeNodes: " << activeNodes << endl;
  cout << "MaxSNR: " << result << " " <<m_ptrMatComputer->computeLog2Det(1.0, mySupStru) <<endl;
#endif
  return result;

}

double
MaxSNRScheduler::ScheduleOneSlot( std::vector<int>& vecSupport, std::vector<double>& vecVariance)
{

  return this->ScheduleOneSlot(vecSupport);
}

bool 
MaxSNRScheduler::CheckFeasible( bool const * const supStru, double txTime2nd)
{
  for (int i = 0; i < m_numMaxHeads; i++) {
    int headName = m_ptrCS->GetVecHeadName()[i];
    int member = -1;
    double interference = 0.0;
    for (int j = 0; j < m_numNodes; j++) {
      if (supStru[j] == true && headName != m_ptrCS->GetChNameByName(j) ) {
        interference += m_ptrMap->GetGijByPair(headName,j) * m_maxPower;
      }
      else if(supStru[j] == true && headName == m_ptrCS->GetChNameByName(j) ){
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
MaxSNRScheduler::CheckFeasible( const vector<int>& supStru, double txTime2nd)
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

