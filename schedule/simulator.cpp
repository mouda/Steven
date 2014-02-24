#include "simulator.h"



Simulator::Simulator(Map* myMap, 
    ClusterStructure* myCS, 
    Scheduler* myScheduler, 
    CORRE_MA_OPE* myField,
    const string& outputFileName):
  m_ptrMap(myMap), 
  m_ptrCS(myCS), 
  m_ptrSched(myScheduler),
  m_ptrGaussianField(myField),
  m_fileHandler(outputFileName)

{
  m_vecSupport = new vector<int>(m_ptrMap->GetNumNodes());
}

Simulator::~Simulator()
{
  list<Slot*>::const_iterator it = m_listSlot.begin();
  vector<double> vecEntropy;
  for (; it != m_listSlot.end(); ++it) {
    vecEntropy.push_back((*it)->GetEntropy());
  }
  m_fileHandler.WriteString(VecToString(vecEntropy));
  this->SequentialFree();
}

void
Simulator::SetEvents(double t_ms)
{

}

void
Simulator::SequentialRun(double t_ms)
{
  std::vector<int> vecSupport(m_ptrMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 0);
  std::vector<double> currVecVariance(m_ptrMap->GetNumNodes());
  fill(currVecVariance.begin(), currVecVariance.end(), 1.0);
  std::vector<double> nextVecVariance(m_ptrMap->GetNumNodes());
  fill(nextVecVariance.begin(), nextVecVariance.end(), 1.0);

  m_ptrSched->ScheduleOneSlot(vecSupport, currVecVariance);
  double entropy = m_ptrGaussianField->GetJointEntropy(vecSupport, currVecVariance, 0.0, m_ptrMap->GetQBits());
  cout << "Entropy: " << entropy << ' ';  
  cout << "Solution: " << VecToString(vecSupport) << endl;
  cout << "Variance: " << VecToString(currVecVariance) << endl;
  m_ptrGaussianField->UpdateVariance(currVecVariance, nextVecVariance, vecSupport, m_ptrSched->GetTxTimePerSlot());

  Slot* ptrCurrSlot = new Slot(vecSupport, nextVecVariance, entropy);
  Slot* ptrNextSlot = 0;
  m_listSlot.push_back(ptrCurrSlot);

  for (double currTime = 0; currTime < t_ms; currTime+=m_ptrSched->GetTxTimePerSlot()) {
    ptrNextSlot = this->GetNextSlot(ptrCurrSlot);
    m_listSlot.push_back(ptrNextSlot);
    ptrCurrSlot = ptrNextSlot;
    ptrNextSlot = 0;
  }
}

Slot*
Simulator::GetNextSlot(Slot* mySlot)
{
  std::vector<int> vecSupport(m_ptrMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 0);
  std::vector<double> nextVecVariance(m_ptrMap->GetNumNodes());
  fill(nextVecVariance.begin(), nextVecVariance.end(), 1.0);
  
  m_ptrSched->ScheduleOneSlot(vecSupport, mySlot->GetVariance());
  double entropy = m_ptrGaussianField->GetJointEntropy(vecSupport, mySlot->GetVariance(), 0.0, m_ptrMap->GetQBits());
  cout << "Entropy: " << entropy << ' ';  
  cout << "Solution: " << VecToString(vecSupport) << endl;
  cout << "Variance: " << VecToString(mySlot->GetVariance()) << endl;
  m_ptrGaussianField->UpdateVariance(mySlot->GetVariance(), nextVecVariance, vecSupport, m_ptrSched->GetTxTimePerSlot());

  Slot* ptrSlot = new Slot(vecSupport, nextVecVariance, entropy);
  return ptrSlot;
}

void
Simulator::SequentialFree()
{
  std::list<Slot*>::iterator it = m_listSlot.begin();
  for (; it != m_listSlot.end(); ++it) {
    if (*it) {
      delete *it; 
    }
  }
}

bool
Simulator::SelfCheck()
{
  if (!m_ptrMap) {
    cerr << "Error: Uninitialized Map" << endl;
    return false;
  }
  if (!m_ptrCS) {
    cerr << "Error: Uninitialized Cluster Structure" << endl;
    return false;
  }
  if (!m_ptrSched) {
    cerr << "Error: Uninitialized Scheduler" << endl;
    return false;
  }
  /* print map */
  cout << "===================== Map ========================" << endl;
  cout << setw(20) << "MapId:" << setw(10) << m_ptrMap->GetMapId() << endl;
  cout << setw(20) << "Nodes:" << setw(10) << m_ptrMap->GetNumNodes() << endl;
  cout << setw(20) << "MaxHeads:" << setw(10) << m_ptrMap->GetNumInitHeads() << endl;
  cout << setw(20) << "Noise:" << setw(10) << m_ptrMap->GetNoise() << endl;
  cout << setw(20) << "IdtEntropy:" << setw(10) << m_ptrMap->GetIdtEntropy() << endl;
  cout << endl;
  cout << "================= Gaussian Field =================" << endl;
  cout << setw(20) << "Variance: " << setw(10) << m_ptrGaussianField->GetVariance() << endl;
  cout << setw(20) << "Spatial Correlation Factor: " << setw(10) << m_ptrGaussianField->GetSpatialCorrelationFactor() << endl;
  cout << setw(20) << "Temporal Correlation Factor: " << setw(10) << m_ptrGaussianField->GetTemporalCorrelationFactor() << endl;
  cout << endl;
  cout << "=============== Cluster Structure  ===============" << endl;
  /* print cluster structure */
  m_ptrCS->Print();
  cout << "==================================================" << endl;
  cout << endl;
  
  return true;
}

void
Simulator::Run()
{
  vector<int> vecSupport(m_ptrMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 0);
//  cout << m_ptrSched->PrintSelf() << endl;
//  cout << m_ptrSched->ScheduleOneSlot(vecSupport) << endl;
}
void 
Simulator::Run(const int numSlots)
{
  fill(m_vecSupport->begin(), m_vecSupport->end(), 0);
  for (int i = 0; i < numSlots; ++i) {
    fill(m_vecSupport->begin(), m_vecSupport->end(), 0);
    cout <<"Entropy: " << setw(8) << m_ptrSched->ScheduleOneSlot(*m_vecSupport)<< ' ';
    cout <<"Solution: " << VecToString(*m_vecSupport) << endl;
  }
}

vector<int>
Simulator::CheckConnection( const vector<int>& vecSupport )
{
  vector<int> returnVector = vecSupport;
  while(!CheckFeasible(vecSupport, m_ptrSched->GetTxTimePerSlot() )){
    for (int i = 0; i < m_ptrMap->GetNumNodes(); i++) {
      if (vecSupport[i] == 1) {
        returnVector[i] = 0;
        break;
      }
    }
  }
  return returnVector;
}

bool
Simulator::CheckFeasible(const vector<int>& supStru, double txTime2nd)
{
  for (int i = 0; i < m_ptrMap->GetNumInitHeads() ; i++) {
    int headName = m_ptrCS->GetVecHeadName()[i];
    int member = -1;
    double interference = 0.0;
    for (int j = 0; j < m_ptrMap->GetNumNodes(); j++) {
      if (supStru[j] == 1 && headName != m_ptrCS->GetChNameByName(j) ) {
        interference += m_ptrMap->GetGijByPair(headName,j) * m_ptrMap->GetMaxPower();
      }
      else if(supStru[j] == 1 && headName == m_ptrCS->GetChNameByName(j) ){
        member = j;
      }
    }
    if (member == -1) continue; 
    if (m_ptrMap->GetIdtEntropy() > m_ptrMap->GetBandwidth() * 
        txTime2nd*log2(1.0 + m_ptrMap->GetMaxPower() * m_ptrMap->GetGijByPair(headName,member)
          / (m_ptrMap->GetNoise() + interference))) {
      return false;
    }
  }
  return true;
}

void
Simulator::Print( const vector<int>& vec)
{
  for (int i = 0; i < vec.size(); ++i) {
    cout << vec[i] << ' ';
  }
}

template< class T>
string
Simulator::VecToString( const vector<T>& vec)
{
  stringstream ss;
  for (int i = 0; i < vec.size(); ++i) {
    ss << vec[i] << ' ';
  }
  return ss.str();
}


double
Simulator::GetTotalEntropy( const vector<int>& vecSupport) const
{
  return 0;
}

double
Simulator::Get1stSlotEntropy( const vector<int>& vecSupport) const
{
  return 0;
}

