#include "minResCsFactory.h"

MinResCsFactory::MinResCsFactory( Map const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  CsFactory(myMap, myMatComputer),
  m_vecClusterHeadBits(myMap->GetNumInitHeads()),
  m_vecClusterHeadMS(myMap->GetNumInitHeads()),
  m_vecClusterHeadWatt(myMap->GetNumInitHeads()), 
  m_matBestCluStru(myMap->GetNumInitHeads(), vector<int>(myMap->GetNumNodes(),0)),
  m_vecBestReceivedInterference(myMap->GetNumInitHeads()),
  m_vecBestSINR_forVerification(myMap->GetNumNodes()),
  m_vecBestBpshz_forVerification(myMap->GetNumNodes()),
  m_powerBest(myMap->GetNumNodes()),
  m_vecBestAllSupStru(myMap->GetNumNodes())
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
  vector<double>  vecPower(m_ptrMap->GetNumNodes());
  fill(vecPower.end(), vecPower.end(), 0.0);
  double sysRedundancy = m_ptrMatComputer->computeLog2Det(1.0, m_ptrCS->GetAllSupStru());
  m_wholeSystemEntropy = m_ptrMap->GetNumInitHeads() * m_ptrMap->GetIdtEntropy() +sysRedundancy;
  CSPowerUpdater myCSPowerUpdater(m_ptrMap);
  bool flagAnsFound = false;
  m_curSupNum = m_ptrMap->GetNumNodes(); 
  m_cur2nd_ms = 0.0;
  
  if( m_curSupNum > m_ptrMap->GetNumInitHeads() ) 
    m_cur2nd_ms = myCSPowerUpdater.Solve_withT2Adj_BinerySearch_2(10.0, vecPower, m_ptrCS->GetAllSupStru(), m_ptrCS);
  else
    m_cur2nd_ms = 0;

  m_cur1st_ms       = Return1stTotalNcal1stResors_HomoPower(); 
  m_cur2nd_Joule    = returnTransientJoule(vecPower);
  m_cur1st_Joule    = m_power1st * m_cur1st_ms/1000.0;
#ifdef DEBUG
  cout << m_cur1st_ms <<" " <<m_cur2nd_ms << endl;
#endif

  m_curSupNum         = m_ptrCS->calSupNodes();
  m_curChNum          = m_ptrMap->GetNumInitHeads();
  m_nextChNum         = m_curChNum;
  m_curJEntropy       = m_curSupNum * m_ptrMap->GetIdtEntropy() + 
                        m_ptrMatComputer->computeLog2Det(1.0, m_ptrCS->GetAllSupStru());
  m_curPayoff         = m_cur1st_ms + m_cur2nd_ms;

  cout << "curPayoff: " << m_curPayoff << endl;
  m_bestAllServeFound = false;
  checkBestClusterStructure_DataCentric(0, myCSPowerUpdater, vecPower);


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

/*
    Check if this structure is the best
*/
bool 
MinResCsFactory::checkBestClusterStructure_DataCentric(int inputRound, CSPowerUpdater& myPowerUpdater,  vector<double>& vecPower)
{
    bool curAllServe = ( m_curJEntropy > m_fidelityRatio * m_wholeSystemEntropy?true:false);
    if( curAllServe ) m_bestAllServeFound = true;
    //cout<<"In check Best"<<endl;
    if (( m_curJEntropy > m_bestFeasibleJEntropy) && !curAllServe && !m_bestAllServeFound)
    {
        m_roundBest             = inputRound;
        m_bestFeasibleJEntropy  = m_curJEntropy;
        m_bestFeasibleSupNum    = m_curSupNum ;
        m_best1st_Joule         = m_cur1st_Joule;
        m_best1st_ms            = m_cur1st_ms;
        m_best2nd_Joule         = m_cur2nd_Joule;
        m_best2nd_ms            = m_cur2nd_ms;
        keepBestStructure( vecPower);
        myPowerUpdater.showVerificationResult(vecPower, 
            m_ptrCS->GetAllSupStru(), 
            m_ptrCS, 
            m_vecBestReceivedInterference, 
            m_vecBestSINR_forVerification, 
            m_vecBestBpshz_forVerification
            );
        m_bestChNum             = m_curChNum;

    }
    else if ((m_curPayoff < m_bestFeasiblePayoff) && curAllServe)
    {
        //cout<<"Find new best"<<endl;
        m_roundBest             = inputRound;
        m_bestFeasibleJEntropy  = m_curJEntropy;
        m_bestFeasibleSupNum    = m_curSupNum ;
        m_bestFeasiblePayoff    = m_curPayoff;
        m_best1st_Joule         = m_cur1st_Joule;
        m_best1st_ms            = m_cur1st_ms;
        m_best2nd_Joule         = m_cur2nd_Joule;
        m_best2nd_ms            = m_cur2nd_ms;
        m_bestChNum             = m_curChNum;
        keepBestStructure( vecPower);
        myPowerUpdater.showVerificationResult(vecPower, 
            m_ptrCS->GetAllSupStru(), 
            m_ptrCS,
            m_vecBestReceivedInterference,
            m_vecBestSINR_forVerification,
            m_vecBestBpshz_forVerification
            );
    }
    //cout<<"best FeasiblePayoff="<<bestFeasiblePayoff<<" with headNum="<<bestChNum<<";Info Ratio="<<bestFeasibleJEntropy/wholeSystemEntopy<<endl;
    return false;
}

double 
MinResCsFactory::Return1stTotalNcal1stResors_HomoPower() 
{
    m_power1st = m_ptrMap->GetMaxPower();
    double T1 = 0;
    list <list<int> >::const_iterator it_LiInt = m_ptrCS->GetListCluMemeber().begin();
    for(unsigned int i=0; i < m_ptrCS->GetVecHeadName().size(); ++i, ++it_LiInt) {
        if ( m_ptrCS->GetVecHeadName().at(i) == -1 || it_LiInt->size() == 0 ) continue;
        double information = it_LiInt->size() * m_ptrMap->GetIdtEntropy() + m_ptrMatComputer->computeLog2Det(1.0, m_ptrCS->GetMatClusterStru().at(i));
        double temp = (information/(m_ptrMap->GetBandwidth() *log2(1+m_power1st*m_ptrMap->GetGi0ByNode(m_ptrCS->GetVecHeadName().at(i)) /m_ptrMap->GetNoise())));
        m_vecClusterHeadBits[i] = information;
        m_vecClusterHeadMS[i]   = temp;
        m_vecClusterHeadWatt[i] = m_power1st;
        T1+=temp;
    }
    return T1;
}

double 
MinResCsFactory::returnTransientJoule( const vector<double>& vecPower ) 
{
    list <list<int> >::const_iterator itlist1 = m_ptrCS->GetListCluMemeber().begin();
    double accuJoule = 0;
    for(int i =0; itlist1 != m_ptrCS->GetListCluMemeber().end(); ++itlist1,++i)
    {
        list<int>::const_iterator it1 = itlist1->begin();
        double tempSize = static_cast<double> (itlist1->size());
        if ( tempSize == 1 ) continue;
        for(; it1 != itlist1->end(); it1++) {
            if (*it1 == m_ptrCS->GetVecHeadName().at(i)) continue;
            accuJoule += vecPower.at(*it1) * m_cur2nd_ms/(tempSize-1);
            //  EYESTEVEN
            /*  cout<<"node "<<setw(4)<<*it1<<" : "<<setw(12)<<nextNodePower[(*it1)]<<"(Watt), "<<setw(12)<<nextNodePower[(*it1)]*m_cur2nd_ms/ \
            //  (tempSize-1)/1000<<"(Joule) " <<" to Head "<<setw(4)<<cSystem->vecHeadName[i]<<" with "<<setw(4)<<(int)tempSize << " in same cluster" <<endl;
            */
        }
    }
    accuJoule/=1000;
    return (accuJoule);
}

void 
MinResCsFactory::keepBestStructure( const vector<double>& vecPower)
{
    m_vecHeadNameBest.assign(m_ptrCS->GetVecHeadName().begin(),m_ptrCS->GetVecHeadName().end());
    m_listCluMemBest.assign(m_ptrCS->GetListCluMemeber().begin(),m_ptrCS->GetListCluMemeber().end());
    //1st tier computation has been done in the calculateMatrics_minResors()
    m_vecBestClusterBits.assign(m_vecClusterHeadBits.begin(), m_vecClusterHeadBits.end());
    m_vecBestClusterSize.assign(m_ptrCS->GetVecClusterSize().begin(), m_ptrCS->GetVecClusterSize().end());
    m_vecBestClusterHeadMS.assign(m_vecClusterHeadMS.begin(), m_vecClusterHeadMS.end());
    m_vecBestClusterHeadWatt.assign(m_vecClusterHeadWatt.begin(), m_vecClusterHeadWatt.end());
    for(int i=0; i < m_ptrMap->GetNumInitHeads(); ++i)
    {
        for(int j=0; j < m_ptrMap->GetNumNodes(); ++j)
        {
            m_matBestCluStru[i][j] = m_ptrCS->GetMatClusterStru().at(i).at(j);
        }
    }
    for(int i =0; i < m_ptrMap->GetNumNodes(); ++i)
    {
        m_powerBest[i] = vecPower.at(i);
        // cout<<"InKEEP "<<i<<" Power"<< powerBest[i]<<endl;
    }
    for(int i=0; i < m_ptrMap->GetNumNodes(); ++i)
        m_vecBestAllSupStru[i] = m_ptrCS->GetAllSupStru().at(i);
}
