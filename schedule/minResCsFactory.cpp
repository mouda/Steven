#include "minResCsFactory.h"
#include <cfloat>

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
  m_vecBestAllSupStru(myMap->GetNumNodes()),
  m_aryFlagHRDone(myMap->GetNumInitHeads(), false)
{
  ULAGENT inputNode;
  for (int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
    inputNode.aryConstructor(i, m_ptrMap->GetNodeXPos(i), m_ptrMap->GetNodeYPos(i));
    m_nodes.push_back(inputNode);
  }

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
  m_iniDone = true;
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
            //  (tempSize-1)/1000<<"(Joule) " <<" to Head "<<setw(4)<<m_ptrCS->GetVecHeadName().at[i]<<" with "<<setw(4)<<(int)tempSize << " in same cluster" <<endl;
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

void 
MinResCsFactory::coolOnce_minResors(CSPowerUpdater& myPowerUpdater, vector<double>& vecPower)
{
  int probAdd = ((m_curSupNum<(m_ptrMap->GetNumNodes())) ?20000 :0);
  int probDiscard = ((m_curSupNum < (m_ptrMap->GetNumInitHeads()+1)) ?0:30000);
  bool checkRotateble=false;//check if there are rotatableSet;
  for(int i=0; i<m_ptrMap->GetNumInitHeads(); i++) {
    //cout<<aryFlagHRDone[i]<<" "<<m_ptrCS->GetVecClusterSize().at(i]<<endl;
    if( m_aryFlagHRDone[i] == false && m_ptrCS->GetVecClusterSize().at(i) > 1 )
      checkRotateble=true;
  }
  int probHeadRotate = (((m_curSupNum>2 * m_ptrMap->GetNumInitHeads() ) && 
        checkRotateble) ?10000:0);//Don't do head rotate if there are only a few m_nodes



  int tmpJoinCan = 0;
  bool chkLessCluster = m_ptrCS->returnIfClusterSmall(m_thresholdd,tmpJoinCan);

  //int probJoin = (chkLessCluster&&m_curJEntropy>(fidelityRatio*wholeSystemEntopy))?tmpJoinCan*50:0;
  int probJoin = (chkLessCluster)?tmpJoinCan*30:0;

  //probJoin=((lastJoinPassAccu>thres2-400)?probJoin:0);

  int probIsoltae=((m_curChNum < m_ptrMap->GetNumInitHeads() )?30:0);
  //probIsoltae=((lastJoinPassAccu>thres2)?probIsoltae:0);


  int sumProb = probAdd + probDiscard + probHeadRotate+probJoin+probIsoltae;
  int eventCursor= (int)((double)rand() / ((double)RAND_MAX + 1) * sumProb);
  m_nextEventFlag=-1;// this flag tell add or discard or Headrotate




  //-------------------------------------//
  //Decide event Flag                    //
  //-------------------------------------//

  if (eventCursor<probAdd) m_nextEventFlag = 1;
  else if (eventCursor<(probAdd+probDiscard)) m_nextEventFlag=2;
  else if (eventCursor<(probAdd+probDiscard+probHeadRotate)) m_nextEventFlag=3;
  else if (eventCursor<(probAdd+probDiscard+probHeadRotate+probJoin)) m_nextEventFlag=4;
  else if (eventCursor<sumProb)m_nextEventFlag=5;
  else
  {
    cout<<"Failure in the random step"<<endl;
    cout<<sumProb<<endl;
    cout<<eventCursor<<endl;
  }
  //-------------------------------------//
  // Start the movement                  //
  //-------------------------------------//  if (m_nextEventFlag==1)//Add
  if (m_nextEventFlag==1)
  {
    if (m_ptrCS->GetListUnSupport().size() == 0) 
      cout<<"Error, it should haven't come in here with empty addlist and add."<<endl;
    else
    {
      decideAdd3i_DC_HeadDetMemRan();
      if(m_targetHeadIndex!=-1&&m_targetNode!=-1){
        //cout<<"add "<<m_targetNode<<" to "<<m_ptrCS->GetVecHeadName().at[m_targetHeadIndex]<<endl;
        addMemberSA(m_targetHeadIndex,m_targetNode);
      }
    }
    m_nextChNum=m_curChNum;

  }
  else if (m_nextEventFlag ==2)
  {
    decideDiscard3b();
    discardMemberSA(m_targetHeadIndex,m_targetNode);
    m_nextChNum=m_curChNum;

    //cout<<"discard "<<m_targetNode+1<<" from "<<m_ptrCS->GetVecHeadName().at[m_targetHeadIndex]+1<<endl;
  }
  else if (m_nextEventFlag ==3)
  {
    decideHeadRotate2i_DC_HeadRanMemDet(myPowerUpdater, vecPower);
    //cout<<"HR "<<m_targetNode+1<<" to Replace "<<m_ptrCS->GetVecHeadName().at[m_targetHeadIndex]+1<<endl;
    rotateHeadSA(m_targetHeadIndex,m_targetNode);
    m_nextJEntropy = m_curJEntropy;
    m_nextSupNum = m_curSupNum;
    m_nextChNum=m_curChNum;

  }
  else if (m_nextEventFlag==4){
    m_JoiningHeadIndex=-1;
    decideHeadJoining4b();
    if (m_JoiningHeadIndex==-1 || m_targetHeadIndex==-1) {
      m_nextJEntropy = m_curJEntropy; // entropy unchanged
      m_nextSupNum = m_curSupNum; //support number unchanged
      return;
    }
    join_fromHeadSA(m_JoiningHeadIndex,m_targetHeadIndex);
    m_nextJEntropy = m_curJEntropy;
    m_nextSupNum = m_curSupNum;
  }
  else if (m_nextEventFlag==5){
    m_isolatedHeadIndex=-1;
    m_IsolateNodeName=-1;
    decideIsolate4b();
    isolateHeadSA(m_IsolateNodeName,m_isolatedHeadIndex,m_targetHeadIndex);
    m_nextJEntropy = m_curJEntropy;
    m_nextSupNum = m_curSupNum;
  }
  else
  {
    cout<<"Error. The random Neighbor event "<<m_nextEventFlag<<" choose is wrong"<<endl;
    cout<<"CursupNum "<<m_curSupNum<<" maxChNUm "<< m_ptrMap->GetNumInitHeads() <<endl;
    cout<<"SumProb "<<sumProb<<endl;
  }
}

/*
    m_targetHeadIndex: CLOSET HEAD NOW. (Maybe The Head with least resource usage)
    m_targetNode: randomly proportional to the indepedent information compare to current set

*/
void 
MinResCsFactory::decideAdd3i_DC_HeadDetMemRan() {
  m_targetHeadIndex=-1;
  m_targetNode=-1;

  double curInfo=m_curSupNum*m_ptrMap->GetIdtEntropy()+ m_ptrMatComputer->computeLog2Det(1.0, m_ptrCS->GetAllSupStru());
  double residualInformation = m_wholeSystemEntropy-curInfo;
  double randomCurs=((double)rand() / ((double)RAND_MAX + 1));
  double chooseCurs=0;
  list <int>::const_iterator it_Int = m_ptrCS->GetListUnSupport().begin();

  for(; it_Int != m_ptrCS->GetListUnSupport().end(); it_Int++) {
    m_ptrCS->SetSupport(*it_Int);
    double tmpIndInfo=(m_curSupNum+1)*m_ptrMap->GetIdtEntropy() + m_ptrMatComputer->computeLog2Det(1.0,m_ptrCS->GetAllSupStru())-curInfo;
    chooseCurs += tmpIndInfo;
    m_ptrCS->SetNotSupport(*it_Int);
    if ((chooseCurs/residualInformation)>randomCurs) {
      m_targetNode=*it_Int;
      break;
    }
  }
  if (m_targetNode==-1)return;

  double maxGain=0;
  //find closet head
  for(unsigned int i=0; i < m_ptrCS->GetVecHeadName().size(); i++) {
    if( m_ptrMap->GetGijByPair(m_targetNode,m_ptrCS->GetVecHeadName().at(i))>maxGain) {
      maxGain = m_ptrMap->GetGijByPair(m_targetNode, m_ptrCS->GetVecHeadName().at(i));
      m_targetHeadIndex = i;
    }
  }
}

/*
   Movement Function
   add new member and adjust the "m_nodes" value

*/
void 
MinResCsFactory::addMemberSA(int inputHeadIndex, int inputMemberName)
{
  m_ptrCS->addMemberCs(inputHeadIndex, inputMemberName, m_iniDone);
  m_nodes[inputMemberName].ptrHead = m_ptrCS->returnHeadPtr(inputHeadIndex);
}

/*
    targetHead: (Randomly) Choose a Head uniformly
    m_targetNode: (Deterministically) find The one cause strongest Interference to others
*/
void 
MinResCsFactory::decideDiscard3b()
{
    //Compute the discardable size
    list<list<int> >::const_iterator itli1 = m_ptrCS->GetListCluMemeber().begin();
    int discardableSize=0;
    for(; itli1!=m_ptrCS->GetListCluMemeber().end(); itli1++)if(itli1->size()>1)discardableSize++;
    assert(discardableSize!=0);
    //Randomly choose a cluster
    itli1 = m_ptrCS->GetListCluMemeber().begin();
    m_targetHeadIndex = -1;
    int chooseCursor = (int) ((double)rand() / ((double)RAND_MAX + 1) * discardableSize)+1;//+1 Becasue the intial might
    assert(chooseCursor<=m_ptrMap->GetNumInitHeads()&&chooseCursor>=0);

    itli1 = m_ptrCS->GetListCluMemeber().begin();
    m_targetHeadIndex=0;
    for(int i=0; i<chooseCursor; itli1++,m_targetHeadIndex++)
    {
        if(itli1->size()>1)i++;
        if(i==chooseCursor)break;
    }

    assert(m_ptrCS->GetVecHeadName().at(m_targetHeadIndex)!=-1);
    assert(m_ptrCS->GetVecClusterSize().at(m_targetHeadIndex)==itli1->size());
    assert(m_ptrCS->GetVecClusterSize().at(m_targetHeadIndex)>1);
    double maxInterferenceGen = -1;
    double tempInter = 0;
    int    maxInteredName = -1;
    list<int>::const_iterator it2=itli1->begin();
    for(; it2!=itli1->end(); it2++)
    {
        tempInter = 0;
        if(m_ptrCS->GetVecHeadName().at(m_targetHeadIndex)==(*it2))continue;
        else
        {
            for(int j=0; j<m_ptrMap->GetNumInitHeads(); j++)
            {
                if(j==m_targetHeadIndex)continue;
                else if(m_ptrCS->GetVecHeadName().at(j)==-1)continue;//Cluster Not Exist;
                else
                {
                    tempInter+= (m_nodes[(*it2)].power * m_ptrMap->GetGijByPair(m_ptrCS->GetVecHeadName().at(j), (*it2)) );
                }
            }
            //   cout<<"tempInter: "<<tempInter<<endl;
            if(tempInter>maxInterferenceGen)
            {
                maxInterferenceGen=tempInter;
                maxInteredName = (*it2);
            }
        }
    }
    assert(maxInteredName!=-1);

    m_targetNode = maxInteredName;
}

/*
   Movement Function
   discard the node form the specific head(cluster)
   */
void 
MinResCsFactory::discardMemberSA(int inputHeadIndex, int inputMemberName)
{
  m_ptrCS->discardMemberCs(inputHeadIndex,inputMemberName);
  m_ptrHeadLastDiscard = m_nodes[inputMemberName].ptrHead;
  m_powerLastDiscard = m_nodes[inputMemberName].power;
  m_nodes[inputMemberName].ptrHead=NULL;
  m_nodes[inputMemberName].power =0;
}

void 
MinResCsFactory::decideHeadRotate2i_DC_HeadRanMemDet(CSPowerUpdater& myPowerUpdater, vector<double>& vecPower)
{
	//----Uniformly choosed
	int rotateAbleSize=0;
	for (int i=0; i<m_ptrMap->GetNumInitHeads(); i++)
		if(m_ptrCS->GetVecClusterSize().at(i)>1)rotateAbleSize++;

	//Notice: no handle of "rotateAbleSize == 0", it handle by coolOnce_minResors
	//Then we choose the rotate target member cluster
	int rotateCursor = (int)((double)rand() / ((double)RAND_MAX + 1) * rotateAbleSize)+1;

	list <list<int> >::const_iterator itlist1 = m_ptrCS->GetListCluMemeber().begin();
	m_targetHeadIndex = 0;
	for(int i=0; i<rotateCursor; itlist1++,m_targetHeadIndex++) //Disconsecutive Candidate
	{
		//if(itlist1->size()>1&&aryFlagHRDone[m_targetHeadIndex]==false)i++;
		if(itlist1->size()>1)i++;
		if(i==rotateCursor)break;//break befor move to next iterator
	}
	if (itlist1->size()<=1)cout<<"error, search rotate  error"<<endl;
	assert(m_targetHeadIndex>=0&&m_targetHeadIndex<m_ptrMap->GetNumInitHeads());



	//Choose m_targetNode
	m_targetNode = m_ptrCS->GetVecHeadName().at(m_targetHeadIndex);
	double temp1st_ms=0;
	double temp2nd_ms=0;
	double temp2tiers_ms=0;
	double test2tiers_ms=DBL_MAX;

	int OriginalNode = m_ptrCS->GetVecHeadName().at(m_targetHeadIndex);
	list<int>::const_iterator it1 = itlist1->begin();
	for(; it1!=itlist1->end(); it1++)//find minimum interference received among nodes in the cluster
 	{
          if((*it1)==OriginalNode)continue;
          m_ptrCS->SetVecHeadNameByIdx(m_targetHeadIndex, *it1);//==Head Rotate
          temp2nd_ms    = myPowerUpdater.Solve_withT2Adj_BinerySearch_2(10, vecPower, m_ptrCS->GetAllSupStru(), m_ptrCS);
          temp1st_ms    = Return1stTotalNcal1stResors_HomoPower();
          temp2tiers_ms = temp1st_ms+temp2nd_ms;
          if(temp2tiers_ms<test2tiers_ms)
          {
            test2tiers_ms=temp2tiers_ms;
            m_targetNode=*it1;
          }
	}
	m_ptrCS->SetVecHeadNameByIdx(m_targetHeadIndex, OriginalNode);
}

/*
   Rotate Member SA
   */
void 
MinResCsFactory::rotateHeadSA(int inputHeadIndex, int inputMemberName)
{
  m_rotatedHeadNameLast = m_ptrCS->GetVecHeadName().at(inputHeadIndex);
  m_ptrCS->SetVecHeadNameByIdx(inputHeadIndex, inputMemberName);
}
