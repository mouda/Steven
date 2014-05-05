#include "minResCsFactory.h"

MinResCsFactory::MinResCsFactory( Map const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  CsFactory(myMap, myMatComputer),
  m_vecClusterHeadBits(myMap->GetNumInitHeads()),
  m_vecClusterHeadMS(myMap->GetNumInitHeads()),
  m_vecClusterHeadWatt(myMap->GetNumInitHeads()) 
{

}

MinResCsFactory::~MinResCsFactory()
{

}

ClusterStructure*
MinResCsFactory::CreateClusterStructure()
{
  if (m_ptrCS == 0) {
    vector<int> myVecHeadNames;
    list<list<int> > myListCluMembers;
    if (Kmedoid(myVecHeadNames, myListCluMembers)) {
      m_ptrCS = new ClusterStructure(m_ptrMap->GetNumNodes(), 
          m_ptrMap->GetNumInitHeads() );
      m_ptrCS->SetRecord(myVecHeadNames, myListCluMembers);
      SASearch();
       
      return m_ptrCS;
    }
    else{
      return NULL;
    }
  }
}

bool
MinResCsFactory::SASearch()
{
  bool inClu[m_ptrMap->GetNumInitHeads()];
  vector<int>     vecSupport(m_ptrMap->GetNumNodes());
  vector<double>  vecPower(m_ptrMap->GetNumNodes());
  fill(vecSupport.begin(), vecSupport.end(), 1);
  fill(vecPower.end(), vecPower.end(), 0.0);
  double sysRedundancy = m_ptrMatComputer->computeLog2Det(1.0, vecSupport);
  m_wholeSystemEntropy = m_ptrMap->GetNumInitHeads() * m_ptrMap->GetIdtEntropy() +sysRedundancy;
  CSPowerUpdater myCSPowerUpdater(m_ptrMap);
  bool flagAnsFound = false;
  int curSupNum = m_ptrMap->GetNumNodes(); 
  double cur2nd_ms = 0.0;
  
  if( curSupNum > m_ptrMap->GetNumInitHeads() ) 
    cur2nd_ms = myCSPowerUpdater.Solve_withT2Adj_BinerySearch_2(10.0, vecPower, vecSupport, m_ptrCS);
  else
    cur2nd_ms = 0;

  return true;
}

bool
MinResCsFactory::Kmedoid( vector<int>& vecHeadNames, list<list<int> >& listCluMembers )
{
  int retryTimes = 0;
  double* tempHeadX  = new double [m_numMaxHeads];
  double* tempHeadY  = new double [m_numMaxHeads];
  int* tempHeadList  = new int [m_numMaxHeads];
  bool convergedFlag = false;
  bool sameHeadFlag = true;
  vector <vector <int> > tempGroup;

  while(sameHeadFlag) {
    sameHeadFlag = false;
    convergedFlag = false;
    //Clear before added members
    if (retryTimes > ( m_numNodes - m_numMaxHeads + 1 ) ) {
      return false;
    }
    for (unsigned  int i=0 ; i<tempGroup.size(); i++)tempGroup[i].clear(); 
    tempGroup.clear();
    for (int i = 0; i < m_numMaxHeads; i++) {
      vector <int> tempV;
      tempGroup.push_back(tempV);
      tempHeadX[i]    = m_ptrMap->GetNodeXPos(i);
      tempHeadY[i]    = m_ptrMap->GetNodeYPos(i);
      tempHeadList[i] = i + retryTimes;
    }
    while(!convergedFlag)  {
      for (unsigned int i=0 ; i<tempGroup.size(); i++) tempGroup[i].clear(); 
      for (int i = 0; i < m_numNodes; i++) {
        double closetDistance = std::numeric_limits<double>::max();
        int closetHeadIndex = -1;
        for(int j = 0; j < m_numMaxHeads; j++ ) {
          double tempDistance = 
            (tempHeadX[j] - m_ptrMap->GetNodeXPos(i)) * 
            (tempHeadX[j] - m_ptrMap->GetNodeXPos(i)) + 
            (tempHeadY[j] - m_ptrMap->GetNodeYPos(i)) * 
            (tempHeadY[j] - m_ptrMap->GetNodeYPos(i));
          if (closetDistance > tempDistance ) {
            closetDistance = tempDistance;
            closetHeadIndex = j;
          }
        }
        tempGroup[closetHeadIndex].push_back(i);
      }
      convergedFlag = true;
#ifdef DEBUG
      for (int i = 0; i < tempGroup.size(); i++) {
        cout << "cluster: " << i <<"-th ";
        for (int j = 0; j < tempGroup[i].size(); j++) {
          cout << tempGroup[i][j] << ' ';
        }
        cout << endl;
      }
#endif
      for(int i=0; i<m_numMaxHeads; i++) {
        float newHx = 0;
        float newHy = 0;
        arma::vec vecDistance = arma::zeros<arma::vec>(tempGroup[i].size());
        for(unsigned int j=0; j<tempGroup[i].size(); j++) {
          float tempDistance = 0.0;
          for (int k = 0; k < tempGroup[i].size(); k++) {
            if ( j == k ) continue;
            tempDistance += 
              sqrt ( (m_ptrMap->GetNodeXPos(tempGroup[i][j]) - m_ptrMap->GetNodeXPos(tempGroup[i][k]) ) * 
                  (m_ptrMap->GetNodeXPos(tempGroup[i][j]) - m_ptrMap->GetNodeXPos(tempGroup[i][k])) + 
                  (m_ptrMap->GetNodeYPos(tempGroup[i][j]) - m_ptrMap->GetNodeYPos(tempGroup[i][k])) * 
                  (m_ptrMap->GetNodeYPos(tempGroup[i][j]) - m_ptrMap->GetNodeYPos(tempGroup[i][k]))) ;
          }
          vecDistance.at(j) = tempDistance; 
        }
        arma::uword idx;
        vecDistance.min(idx);
        newHx = m_ptrMap->GetNodeXPos(tempGroup[i][idx]);
        newHy = m_ptrMap->GetNodeYPos(tempGroup[i][idx]);
        if ( (abs(newHx-tempHeadX[i]) > 0.01) || (abs(newHy-tempHeadY[i])>0.01) ) convergedFlag = false; 
        tempHeadX[i] = newHx;
        tempHeadY[i] = newHy;
        tempHeadList[i] = tempGroup[i][idx];
      }
    }
    for (int i=0; i<m_numMaxHeads; i++)
      for (int j=i+1; j<m_numMaxHeads; j++)
        if (tempHeadList[i] == tempHeadList[j]) sameHeadFlag = true;
    retryTimes++;
  }

  vecHeadNames.assign(tempHeadList, tempHeadList + m_numMaxHeads);
  listCluMembers.resize(m_numMaxHeads);
  list<list<int> >::iterator iterRows = listCluMembers.begin();
  for (int i = 0; iterRows != listCluMembers.end(); ++iterRows, ++i) {
    iterRows->assign(tempGroup[i].begin(), tempGroup[i].end());
  }

#ifdef DEBUG 
  iterRows = listCluMembers.begin();
  for (int i = 0; iterRows != listCluMembers.end(); ++iterRows, ++i) {
    list<int>::iterator iterCols = iterRows->begin();
    cout << "cluster: " << i <<"-th ";
    for (; iterCols != iterRows->end(); ++iterCols) {
      cout << *iterCols << ' '; 
    }
    cout << endl;
  }
  for (int i = 0; i < m_numMaxHeads; ++i) {
    cout << tempHeadList[i] << ' ';
  }
  cout << endl;
#endif

  delete [] tempHeadX;
  delete [] tempHeadY;
  delete [] tempHeadList; 

  return true;

}

double 
MinResCsFactory::Return1stTotalNcal1stResors_HomoPower() {
    double power1st = m_ptrMap->GetMaxPower();
    double T1=0;
    list <list<int> >::const_iterator it_LiInt = m_ptrCS->GetListCluMemeber().begin();
    for(unsigned int i=0; i < m_ptrCS->GetVecHeadName().size(); ++i, ++it_LiInt) {
        if ( m_ptrCS->GetVecHeadName().at(i) == -1 || it_LiInt->size() == 0 ) continue;
        double information = it_LiInt->size() * m_ptrMap->GetIdtEntropy() + m_ptrMatComputer->computeLog2Det(1.0, m_ptrCS->GetMatClusterStru().at(i));
        double temp = (information/(m_ptrMap->GetBandwidth() *log2(1+power1st*m_ptrMap->GetGi0ByNode(m_ptrCS->GetVecHeadName().at(i)) /m_ptrMap->GetNoise())));
        m_vecClusterHeadBits[i] = information;
        m_vecClusterHeadMS[i]   = temp;
        m_vecClusterHeadWatt[i] = power1st;
        T1+=temp;
    }
    return T1;
}
