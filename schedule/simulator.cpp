#include "simulator.h"

Simulator::Simulator()
{

}

Simulator::Simulator(Map* myMap, ClusterStructure* myCS, 
    Scheduler* myScheduler):
  m_ptrMap(myMap), m_ptrCS(myCS), m_ptrSched(myScheduler)
{

}

Simulator::Simulator(Map* myMap, ClusterStructure* myCS, 
    Scheduler* myScheduler, CORRE_MA_OPE* myField):
  m_ptrMap(myMap), m_ptrCS(myCS), m_ptrSched(myScheduler),
  m_ptrGaussianField(myField)
{
  m_vecSupport = new vector<int>(m_ptrMap->GetNumNodes());

}

Simulator::~Simulator()
{

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
  cout << setw(20) << "Correlation Factor: " << setw(10) << m_ptrGaussianField->GetCorrationFactor() << endl;
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
  cout << m_ptrSched->PrintSelf() << endl;
  cout << m_ptrSched->ScheduleOneSlot(vecSupport) << endl;
}
void 
Simulator::Run(const int numSlots)
{

  fill(m_vecSupport->begin(), m_vecSupport->end(), 0);
  for (int i = 0; i < numSlots; ++i) {
    fill(m_vecSupport->begin(), m_vecSupport->end(), 0);
    cout <<"Entropy: " << setw(8) <<m_ptrSched->ScheduleOneSlot(*m_vecSupport)<< ' ';
    cout <<"Solution: " << toString(*m_vecSupport) << endl;
  }
}

void
Simulator::Print( const vector<int>& vec)
{
  for (int i = 0; i < vec.size(); ++i) {
    cout << vec[i] << ' ';
  }
}

string
Simulator::toString( const vector<int>& vec)
{
  stringstream ss;
  for (int i = 0; i < vec.size(); ++i) {
    ss << vec[i] << ' ';
  }
  return ss.str();
}
