#include "minPowerCsFactory.h"
#define SA_INI_TEMP 3.0
#define SA_FIN_TEMP 0.1

MinPowerCsFactory::MinPowerCsFactory( Map const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  CsFactory(myMap, myMatComputer),
  m_fid(NULL),
  m_mapFileName(""),
  m_compressionRatio(-1.0)
{
}

MinPowerCsFactory::~MinPowerCsFactory()
{
  fclose(m_fid);
  if (m_ptrToolSA != NULL) {
    delete m_ptrToolSA;
  }
}

ClusterStructure*
MinPowerCsFactory::CreateClusterStructure()
{
  double SAIter = 5000;
  double alpha = pow (10, -log10(SA_INI_TEMP/SA_FIN_TEMP)/SAIter);
  m_fid = fopen(m_mapFileName.c_str(), "r");
  if(m_fid == NULL) {
    cerr << "Read map file failed!!" << endl;
    return NULL;
  }
  else if ( m_compressionRatio < 0.0){
    cerr << "Compression Not Initialized" << endl;
    return NULL;

  }
  m_ptrToolSA = new MinPowerSACluster(
      m_fid, 
      m_ptrMap->GetNumNodes(), 
      m_ptrMap->GetNumInitHeads(), 
      SAIter, 
      false, 
      false, 
      SA_INI_TEMP, 
      alpha, 
      m_compressionRatio, 
      "NULL",
      m_ptrMap,
      m_ptrMatComputer,
      m_tier1TxTime,
      m_tier2NumSlot
      );
  if(!m_ptrToolSA->setSystem(m_ptrMap->GetMaxPower(), 
        (int)m_ptrMap->GetQBits(), 
        m_ptrMap->GetBandwidth(), 
        m_fidelityRatio)
      ) {
    cerr << "Set parameter failed! " << endl;
    return NULL;
  }
  char iniFlag[]="GraphPartition";
  if (!m_ptrToolSA->setInitialStucture (iniFlag)) {
    cerr<<"The MinRes can't be initialized! " << endl;
    return NULL;
  }

  if (m_ptrCS == 0) {
    std::vector<int> myVecHeadNames;
    std::list<std::list<int> > myListCluMembers;
    m_ptrToolSA->startCool();
    m_ptrCS = new ClusterStructure(
        m_ptrMap->GetNumNodes(), 
        m_ptrMap->GetNumInitHeads() 
        );
    m_ptrCS->SetRecord(
        m_ptrToolSA->GetVecHeadName(), 
        m_ptrToolSA->GetListCluMemeber(),
        m_ptrToolSA->GetAllSupStru()
        );
    return m_ptrCS;
  }
}

bool
MinPowerCsFactory::Kmedoid( std::vector<int>& vecHeadNames, std::list<std::list<int> >& listCluMembers )
{
  int retryTimes = 0;
  double* tempHeadX  = new double [m_numMaxHeads];
  double* tempHeadY  = new double [m_numMaxHeads];
  int* tempHeadList  = new int [m_numMaxHeads];
  bool convergedFlag = false;
  bool sameHeadFlag = true;
  std::vector <std::vector <int> > tempGroup;

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
      std::vector <int> tempV;
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
  std::list<std::list<int> >::iterator iterRows = listCluMembers.begin();
  for (int i = 0; iterRows != listCluMembers.end(); ++iterRows, ++i) {
    iterRows->assign(tempGroup[i].begin(), tempGroup[i].end());
  }

#ifdef DEBUG 
  iterRows = listCluMembers.begin();
  for (int i = 0; iterRows != listCluMembers.end(); ++iterRows, ++i) {
    std::list<int>::iterator iterCols = iterRows->begin();
    std::cout << "cluster: " << i <<"-th ";
    for (; iterCols != iterRows->end(); ++iterCols) {
      std::cout << *iterCols << ' '; 
    }
    std::cout << endl;
  }
  for (int i = 0; i < m_numMaxHeads; ++i) {
    std::cout << tempHeadList[i] << ' ';
  }
  std::cout << endl;
#endif

  delete [] tempHeadX;
  delete [] tempHeadY;
  delete [] tempHeadList; 

  return true;

}

