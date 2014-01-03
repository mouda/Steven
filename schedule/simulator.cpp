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

}
