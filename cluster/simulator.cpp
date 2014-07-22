#include "simulator.h"
#include "imagePowerUpdater.h"
#include <cstdio>
#include <sstream>
#include <numeric>




Simulator::Simulator(Map* myMap, 
    ClusterStructure* myCS, 
    CORRE_MA_OPE* myField
    ):
  m_ptrMap(myMap), 
  m_ptrCS(myCS), 
  m_ptrGField(myField)
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
  delete m_vecSupport;
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


void
Simulator::WriteWorseCaseTier2Power( const string& fileName, 
    const double TxTimePerSlot, 
    const double maxPower) 
{
  std::list<list<int> >::const_iterator iterRow = m_ptrCS->GetListCluMemeber().begin();
  std::vector<int> vecSolution(m_ptrMap->GetNumNodes(), 0);
  for (; iterRow != m_ptrCS->GetListCluMemeber().end(); ++iterRow) {
    std::list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      vecSolution.at(*iterCol) = 1;
    }
  }

  std::vector<double> myVecPower(m_ptrMap->GetNumNodes());
  ImagePowerUpdater myImagePowerUpdater(m_ptrMap, m_ptrCS,  TxTimePerSlot);
  myImagePowerUpdater.Solve(myVecPower, vecSolution);
  double totalPower = std::accumulate(myVecPower.begin(), myVecPower.end(), 0.0);
  cout << totalPower << endl;
  return;
}



// -------------------------------------------------------------------------- //
// @Description: for the output
// @Provides: 
// -------------------------------------------------------------------------- //


void
Simulator::WriteCS( const string& fileName )
{
    FILE  *fid; 
    fid = fopen(fileName.c_str(), "w");
    fprintf( fid,"%d\n", m_ptrMap->GetNumNodes()   );
    fprintf( fid,"%d\n", m_ptrCS->GetNumHeads()    );
    fprintf( fid,"%e\n", m_ptrMap->GetMaxPower()   );
    fprintf( fid,"%e\n", m_ptrMap->GetIdtEntropy() );
    fprintf( fid,"%e\n", m_ptrCS->GetTier1TotalPower());  // Payoff (objective)
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

