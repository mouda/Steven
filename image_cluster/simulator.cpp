#include "simulator.h"
#include "newImagePowerUpdater.h"
#include <cstdio>
#include <sstream>
#include <numeric>




Simulator::Simulator(
    ImageMap* myMap, 
    ClusterStructure* myCS 
    ):
  m_ptrMap(myMap), 
  m_ptrCS(myCS) 
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
    cerr << "Error: Uninitialized ImageMap" << endl;
    return false;
  }
  if (!m_ptrCS) {
    cerr << "Error: Uninitialized Cluster Structure" << endl;
    return false;
  }
  /* print map */
  cout << "===================== ImageMap ========================" << endl;
  cout << setw(20) << "ImageMapId:" << setw(10) << m_ptrMap->GetMapId() << endl;
  cout << setw(20) << "Nodes:" << setw(10) << m_ptrMap->GetNumNodes() << endl;
  cout << setw(20) << "MaxHeads:" << setw(10) << m_ptrMap->GetNumInitHeads() << endl;
  cout << setw(20) << "Noise:" << setw(10) << m_ptrMap->GetNoise() << endl;
  cout << endl;
  cout << "=============== Cluster Structure  ===============" << endl;
  /* print cluster structure */
  m_ptrCS->Print();
  cout << "==================================================" << endl;
  cout << endl;
  
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
    const double TxTimeNumSlot,
    const double maxPower) 
{
  stringstream ss1;
  ss1 << TxTimePerSlot;
  FileHandler myFile(string("data/")+fileName+string("_")+ss1.str()+string(".out"), "app");
  std::list<list<int> >::const_iterator iterRow = m_ptrCS->GetListCluMemeber().begin();
  std::vector<int> vecSolution(m_ptrMap->GetNumNodes(), 0);
  for (; iterRow != m_ptrCS->GetListCluMemeber().end(); ++iterRow) {
    std::list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      vecSolution.at(*iterCol) = 1;
    }
  }

  std::vector<double> myVecPower(m_ptrMap->GetNumNodes());
  NewImagePowerUpdater myNewImagePowerUpdater(m_ptrMap, m_ptrCS,  TxTimePerSlot, TxTimeNumSlot);
  myNewImagePowerUpdater.Solve(myVecPower, vecSolution);
  double totalPower = std::accumulate(myVecPower.begin(), myVecPower.end(), 10e-30);
//  for (int i = 0; i < myVecPower.size(); ++i) {
//    cout << myVecPower.at(i) << endl;
//  }
  stringstream ss2;
  ss2 << totalPower;
  cout << totalPower << endl;
  myFile.WriteString(ss2.str());
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
    fprintf( fid,"%e\n", 0.0 );
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

