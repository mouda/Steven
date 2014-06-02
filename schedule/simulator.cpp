#include "simulator.h"
#include <cstdio>



Simulator::Simulator(Map* myMap, 
    ClusterStructure* myCS, 
    Scheduler* myScheduler, 
    CORRE_MA_OPE* myField,
    const string& entropyFName,
    const string& MSEFName,
    const string& solutionFName,
    const string& supportNumFName,
    const string& totalEntropyFName ):
  m_ptrMap(myMap), 
  m_ptrCS(myCS), 
  m_ptrSched(myScheduler),
  m_ptrGField(myField),
  m_entropyFHandler(entropyFName),
  m_MSEFHandler(MSEFName),
  m_solutionFHandler(solutionFName),
  m_supportNumFHandler(supportNumFName),
  m_totalEntropyFHandler(totalEntropyFName)
{
  m_vecSupport = new vector<int>(m_ptrMap->GetNumNodes());
  m_vecTotal.resize(m_ptrMap->GetNumNodes());
  for (int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
    if (m_ptrCS->GetChNameByName(i) != i) {
      m_vecTotal.at(i) = 1;
    }
    else {
      m_vecTotal.at(i) = 0;
    }
  }
}

Simulator::~Simulator()
{
  /* free the memory */
  this->SequentialFree();
}

void
Simulator::SetEvents(double t_ms)
{

}

void
Simulator::SequentialRun(int tier2NumSlot)
{
  std::vector<int> vecSupport(m_ptrMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 0);
  std::vector<double> currVecVariance(m_ptrMap->GetNumNodes());
  fill(currVecVariance.begin(), currVecVariance.end(), 1.0);
  std::vector<double> nextVecVariance(m_ptrMap->GetNumNodes());
  fill(nextVecVariance.begin(), nextVecVariance.end(), 1.0);
  std::vector<double> vecPower(m_ptrMap->GetNumNodes());
  fill(vecPower.begin(), vecPower.end(), 0.0);
  std::vector<int> vecSlots(m_ptrMap->GetNumNodes());
  fill(vecSlots.begin(), vecSlots.end(), 0);

  m_ptrSched->ScheduleOneSlot(vecSupport, vecPower, currVecVariance);
  double entropy      = m_ptrGField->GetJointEntropy(vecSupport, currVecVariance, 0.0, m_ptrMap->GetQBits());
  double MSE          = m_ptrGField->GetRateDistortion(vecSupport, currVecVariance, 0.0, m_ptrMap->GetQBits());
  double totalEntropy = m_ptrGField->GetJointEntropy(m_vecTotal, currVecVariance, 0.0, m_ptrMap->GetQBits());
  cout << "Entropy: " << entropy << " MSE: " << MSE << " Total: " << totalEntropy << ' ';  
  cout << "Solution: " << VecToString(vecSupport) << endl;
//  cout << "Variance: " << VecToString(currVecVariance) << endl;
  m_ptrGField->UpdateVariance(currVecVariance, nextVecVariance, vecSupport, vecSlots ,m_ptrSched->GetTxTimePerSlot());

  Slot* ptrCurrSlot = new Slot(vecSupport, nextVecVariance, entropy, totalEntropy, MSE);
  Slot* ptrNextSlot = 0;
  m_listSlot.push_back(ptrCurrSlot);

  for (int currSlot = 0; currSlot < tier2NumSlot - 1; ++currSlot) {
    ptrNextSlot = this->GetNextSlot(ptrCurrSlot, vecSlots);
    m_listSlot.push_back(ptrNextSlot);
    ptrCurrSlot = ptrNextSlot;
    ptrNextSlot = 0;
  }
}

Slot*
Simulator::GetNextSlot(Slot* mySlot, std::vector<int>& vecSlots)
{
  std::vector<int> vecSupport(m_ptrMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 0);
  std::vector<double> nextVecVariance(m_ptrMap->GetNumNodes());
  fill(nextVecVariance.begin(), nextVecVariance.end(), 1.0);
  std::vector<double> vecPower(m_ptrMap->GetNumNodes());
  fill(vecPower.begin(), vecPower.end(), 0.0);
  
  m_ptrSched->ScheduleOneSlot(vecSupport, vecPower, mySlot->GetVariance());
  double entropy      = m_ptrGField->GetJointEntropy(vecSupport, mySlot->GetVariance(), 0.0, m_ptrMap->GetQBits());
  double MSE          = m_ptrGField->GetRateDistortion(vecSupport, mySlot->GetVariance(), 0.0, m_ptrMap->GetQBits());
  double totalEntropy = m_ptrGField->GetJointEntropy(m_vecTotal, mySlot->GetVariance(), 0.0, m_ptrMap->GetQBits());
  cout << "Entropy: " << entropy << " MSE: " << MSE << " Total: " << totalEntropy << ' ';  
  cout << "Solution: " << VecToString(vecSupport) << endl;
//  cout << "Variance: " << VecToString(mySlot->GetVariance()) << endl;
  m_ptrGField->UpdateVariance(mySlot->GetVariance(), nextVecVariance, vecSupport, vecSlots, m_ptrSched->GetTxTimePerSlot());

  Slot* ptrSlot = new Slot(vecSupport, nextVecVariance, entropy, totalEntropy, MSE);
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
  cout << setw(20) << "Variance: " << setw(10) << m_ptrGField->GetVariance() << endl;
  cout << setw(20) << "Spatial Correlation Factor: " << setw(10) << m_ptrGField->GetSpatialCorrelationFactor() << endl;
  cout << setw(20) << "Temporal Correlation Factor: " << setw(10) << m_ptrGField->GetTemporalCorrelationFactor() << endl;
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

// -------------------------------------------------------------------------- //
// @Description: for the output
// @Provides: 
// -------------------------------------------------------------------------- //


void
Simulator::WriteCS( const string& fileName )
{
    FILE  *fid; 
    fid = fopen(fileName.c_str(), "a");
    fprintf( fid,"%d\n", m_ptrMap->GetNumNodes()   );
    fprintf( fid,"%d\n", m_ptrCS->GetNumHeads()    );
    fprintf( fid,"%e\n", m_ptrMap->GetMaxPower()   );
    fprintf( fid,"%e\n", m_ptrMap->GetIdtEntropy() );
    fprintf( fid,"%d\n", 0                         );  // Payoff (objective)
    fprintf( fid,"%d\n", 0                         );  // SA iteration
    fprintf( fid,"%5e\n", 0.0                      );  // 1st tier tx time (ms)
    fprintf( fid,"%5e\n", 0.0                      );  // 2nd tier tx time (ms)
    fprintf( fid,"%5e\n", 0.0                      );  // 1st tier Joule
    fprintf( fid,"%5e\n", 0.0                      );  // 2nd tier Joule


    for ( int i = 0; i < m_ptrCS->GetNumHeads(); ++i) {
        for( int j = 0; j < m_ptrMap->GetNumNodes(); ++j) {
          if (m_ptrCS->GetChIdxByName(j) == i ) 
            fprintf(fid,"%d ", 1);
          else 
            fprintf(fid,"%d ", 0);
        }
        fprintf(fid,"\n");
    }

    for ( int i = 0; i < m_ptrCS->GetNumHeads(); ++i) {
        fprintf(fid,"%d ", m_ptrCS->GetVecHeadName().at(i)+1);
    }
    fprintf(fid,"\n");

    /* updated power for each node */
    for ( int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
        fprintf(fid,"%E ", 0);
    }
    fprintf(fid,"\n");

    list<list<int> >::const_iterator itRows = m_ptrCS->GetListCluMemeber().begin();
    for (; itRows != m_ptrCS->GetListCluMemeber().end(); ++itRows ) {
        list <int>::const_iterator itCols = itRows->begin();
        for (; itCols!=itRows->end(); itCols++) fprintf(fid,"%d ", (*itCols)+1 );
        fprintf(fid,"\n");
    }

    fprintf(fid,"\n");
    fclose(fid);
}

void
Simulator::WriteEntropy()
{
  list<Slot*>::const_iterator it = m_listSlot.begin();
  vector<double> vecEntropy;
  for (; it != m_listSlot.end(); ++it) {
    vecEntropy.push_back((*it)->GetEntropy());
  }
  m_entropyFHandler.WriteString(VecToString(vecEntropy));
}


void
Simulator::WriteMSE()
{
  vector<double> vecMSE;
  list<Slot*>::const_iterator it = m_listSlot.begin();
  for (; it != m_listSlot.end(); ++it) {
    vecMSE.push_back((*it)->GetMSE());
  }
  m_MSEFHandler.WriteString(VecToString(vecMSE));
}

void 
Simulator::WriteSolution()
{
  list<Slot*>::const_iterator it = m_listSlot.begin();
  for (; it != m_listSlot.end(); ++it) {
    m_solutionFHandler.WriteString(VecToString((*it)->GetSupport()));
  }
}

void
Simulator::WriteSupport()
{
  list<Slot*>::const_iterator it = m_listSlot.begin();
  vector<int> vecSupportNumber;
  for (; it != m_listSlot.end(); ++it) {
    vector<int> vecTmp = (*it)->GetSupport();
    int counter = 0;
    for (int i = 0; i < vecTmp.size(); ++i) {
      if ( vecTmp.at(i) == 1) {
        ++counter;
      }
    }
    vecSupportNumber.push_back(counter);
  }
  m_supportNumFHandler.WriteString(VecToString(vecSupportNumber));
}

void
Simulator::WriteTotalEntropy()
{
  list<Slot*>::const_iterator it = m_listSlot.begin();
  vector<double> vecEntropy;
  for (; it != m_listSlot.end(); ++it) {
    vecEntropy.push_back((*it)->GetTotalEntropy());
  }
  m_totalEntropyFHandler.WriteString(VecToString(vecEntropy));
}
