#include "simulator.h"
#include <cstdio>



Simulator::Simulator(ImageMap* myMap, 
    ClusterStructure* myCS, 
    Scheduler* myScheduler, 
    ImageSource* ptrImageSource,
    const string& entropyFName,
    const string& MSEFName,
    const string& solutionFName,
    const string& supportNumFName,
    const string& totalEntropyFName,
    const string& powerFName):
  m_ptrImageMap(myMap), 
  m_ptrCS(myCS), 
  m_ptrSched(myScheduler),
  m_ptrImageSource(ptrImageSource),
  m_entropyFHandler(entropyFName),
  m_MSEFHandler(MSEFName),
  m_solutionFHandler(solutionFName),
  m_supportNumFHandler(supportNumFName),
  m_totalEntropyFHandler(totalEntropyFName),
  m_vecPowerFHandler(powerFName)
{
  m_vecSupport = new vector<int>(m_ptrImageMap->GetNumNodes());
  m_vecTotal.resize(m_ptrImageMap->GetNumNodes());
  for (int i = 0; i < m_ptrImageMap->GetNumNodes(); ++i) {
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
Simulator::SequentialRun(int tier2NumSlot)
{
  std::vector<int> vecSupport(m_ptrImageMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 0);
  std::vector<double> currVecVariance(m_ptrImageMap->GetNumNodes());
  fill(currVecVariance.begin(), currVecVariance.end(), 1.0);
  std::vector<double> nextVecVariance(m_ptrImageMap->GetNumNodes());
  fill(nextVecVariance.begin(), nextVecVariance.end(), 1.0);
  std::vector<double> vecPower(m_ptrImageMap->GetNumNodes());
  fill(vecPower.begin(), vecPower.end(), 0.0);
  std::vector<int> vecSlots(m_ptrImageMap->GetNumNodes());
  fill(vecSlots.begin(), vecSlots.end(), 0);
  double currPower   = 0; 
  double entropy     = 0;
  double MSE         = 0;  
  double totalEntropy= 0;
  currPower    = m_ptrSched->ScheduleOneSlot(vecSupport, vecPower, currVecVariance);
  if (m_ptrImageSource) {
    entropy      = m_ptrImageSource->GetJointEntropy(vecSupport, currVecVariance, 0.0, 0.0);
    totalEntropy = m_ptrImageSource->GetJointEntropy(m_vecTotal, currVecVariance, 0.0, 0.0);
  }
  cout << "Entropy: " << entropy << " MSE: " << MSE << " Total: " << totalEntropy << ' ';  
  cout << "Solution: " << VecToString(vecSupport) << endl;
  cout << "Power: " << VecToString(vecPower) << endl;
//  cout << "Variance: " << VecToString(currVecVariance) << endl;

  Slot* ptrCurrSlot = new Slot(vecSupport, nextVecVariance, vecPower, entropy, totalEntropy, MSE, currPower);
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
  std::vector<int> vecSupport(m_ptrImageMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 0);
  std::vector<double> nextVecVariance(m_ptrImageMap->GetNumNodes());
  fill(nextVecVariance.begin(), nextVecVariance.end(), 1.0);
  std::vector<double> vecPower(m_ptrImageMap->GetNumNodes());
  fill(vecPower.begin(), vecPower.end(), 0.0);
  double currPower   = 0; 
  double entropy     = 0;
  double MSE         = 0;  
  double totalEntropy= 0;
  currPower    = m_ptrSched->ScheduleOneSlot(vecSupport, vecPower, mySlot->GetVariance());
  if (m_ptrImageSource) {
    entropy      = m_ptrImageSource->GetJointEntropy(vecSupport, mySlot->GetVariance(), 0.0, 0.0);
    totalEntropy = m_ptrImageSource->GetJointEntropy(m_vecTotal, mySlot->GetVariance(), 0.0, 0.0);
  }
  
  cout << "Entropy: " << entropy << " MSE: " << MSE << " Total: " << totalEntropy << ' ';  
  cout << "Solution: " << VecToString(vecSupport) << endl;
  cout << "Power: " << VecToString(vecPower) << endl;
//  cout << "Variance: " << VecToString(mySlot->GetVariance()) << endl;

  Slot* ptrSlot = new Slot(vecSupport, nextVecVariance, vecPower, entropy, totalEntropy, MSE, currPower);
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
  if (!m_ptrImageMap) {
    cerr << "Error: Uninitialized ImageMap" << endl;
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
  cout << "===================== ImageMap ========================" << endl;
  cout << setw(20) << "ImageMapId:" << setw(10) << m_ptrImageMap->GetMapId() << endl;
  cout << setw(20) << "Nodes:" << setw(10) << m_ptrImageMap->GetNumNodes() << endl;
  cout << setw(20) << "MaxHeads:" << setw(10) << m_ptrImageMap->GetNumInitHeads() << endl;
  cout << setw(20) << "Noise:" << setw(10) << m_ptrImageMap->GetNoise() << endl;
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
  vector<int> vecSupport(m_ptrImageMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 0);
//  cout << m_ptrSched->PrintSelf() << endl;
//  cout << m_ptrSched->ScheduleOneSlot(vecSupport) << endl;
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
    fprintf( fid,"%d\n", m_ptrImageMap->GetNumNodes()   );
    fprintf( fid,"%d\n", m_ptrCS->GetNumHeads()    );
    fprintf( fid,"%e\n", m_ptrImageMap->GetMaxPower()   );
    fprintf( fid,"%e\n", 0.0                       );
    fprintf( fid,"%d\n", 0                         );  // Payoff (objective)
    fprintf( fid,"%d\n", 0                         );  // SA iteration
    fprintf( fid,"%5e\n", 0.0                      );  // 1st tier tx time (ms)
    fprintf( fid,"%5e\n", 0.0                      );  // 2nd tier tx time (ms)
    fprintf( fid,"%5e\n", 0.0                      );  // 1st tier Joule
    fprintf( fid,"%5e\n", 0.0                      );  // 2nd tier Joule


    for ( int i = 0; i < m_ptrCS->GetNumHeads(); ++i) {
        for( int j = 0; j < m_ptrImageMap->GetNumNodes(); ++j) {
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
    for ( int i = 0; i < m_ptrImageMap->GetNumNodes(); ++i) {
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

void
Simulator::WriteVecPower()
{
  list<Slot*>::const_iterator it = m_listSlot.begin();
  for (; it != m_listSlot.end(); ++it) {
    m_vecPowerFHandler.WriteString(VecToString((*it)->GetVecPower()));
  }
}
