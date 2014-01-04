#include "maxSNRScheduler.h"

MaxSNRScheduler::MaxSNRScheduler( const double txTime, 
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE const * const ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("MaxSNR")
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads =  m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();

}

MaxSNRScheduler::~MaxSNRScheduler()
{

}

double
MaxSNRScheduler::ScheduleOneSlot( vector<int>& vecSupport )
{
//  bool * mySupStru;
//  mySupStru = new bool [m_numNodes];
//  for (int i = 0; i < m_numNodes; i++) {
//    mySupStru[i] = false;
//  }
  for (int i = 1; i < m_numMaxHeads; i++) {
    double maxRxPower = 0.0;
    int headName = m_ptrCS->GetVecHeadName()[i];
    int maxRxPowerNode = -1;
    for (int j = 0; j < m_numNodes; j++) {
      if (m_ptrCS->GetChNameByName(j) == j ) continue; 
      if (m_ptrCS->GetChIdxByName(j) == i ) {
        if (vecSupport[j] == 1) {
          vecSupport[j] = 0;
          continue;
        }
        if (m_ptrMap->GetGijByPair(headName,j) * m_maxPower > maxRxPower) {
          maxRxPower = m_ptrMap->GetGijByPair(headName,j) * m_maxPower; 
          maxRxPowerNode = j;
        }
      }
    }
//    mySupStru[maxRxPowerNode] = true; 
    vecSupport[maxRxPowerNode] = 1;

  }

  /* test the Interference */
  while(!CheckFeasible(vecSupport, m_txTimePerSlot)){
    for (int i = 0; i < m_numNodes; i++) {
//      if (mySupStru[i] == true) {
      if (vecSupport[i] == 1) {
//        mySupStru[i] = false;
        vecSupport[i] = 0;
        break;
      }
    }
  }
  int activeNodes = 0;
  for (int i = 0; i < m_numNodes; i++) {
    //if (mySupStru[i] == true) ++activeNodes;
    if (vecSupport[i] == 1) ++activeNodes;
  }

  double result = activeNodes * m_ptrMap->GetIdtEntropy() + m_ptrMatComputer->computeLog2Det(1.0, vecSupport);
#ifdef DEBUG
  cout << "activeNodes: " << activeNodes << endl;
  cout << "MaxSNR: " << result << " " <<m_ptrMatComputer->computeLog2Det(1.0, mySupStru) <<endl;
#endif
//  delete [] mySupStru;
  return result;

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

