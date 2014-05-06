#include "csPowerUpdater.h"
#include <cmath>
#include <cassert>
#include <cfloat>

CSPowerUpdater::CSPowerUpdater(
    Map const * const ptrMap 
    ):
  m_threshold(1e-12),
  m_ptrMap(ptrMap), 
  m_idtEntropy(ptrMap->GetIdtEntropy())
{
  m_maIndexInterference = new int* [m_ptrMap->GetNumInitHeads()];
  for (int i=0; i < m_ptrMap->GetNumInitHeads(); i++) 
    m_maIndexInterference[i] = new int [m_ptrMap->GetNumInitHeads()];

  m_maStrengthInterference = new double* [m_ptrMap->GetNumInitHeads()];
  for (int i=0; i < m_ptrMap->GetNumInitHeads(); i++) 
    m_maStrengthInterference[i] = new double [m_ptrMap->GetNumInitHeads()];
  m_inBandNoise = m_ptrMap->GetNoise() * m_scale;
//  cout << "********idtEntropy " << m_idtEntropy << endl;
//  cout << "********bandwidth " << m_ptrMap->GetBandwidth() << endl;
}

CSPowerUpdater::~CSPowerUpdater()
{
  for(int i = 0; i < m_ptrMap->GetNumInitHeads(); i++) delete[] m_maIndexInterference[i];
  delete [] m_maIndexInterference;
  for(int i = 0; i < m_ptrMap->GetNumInitHeads(); i++) delete[] m_maStrengthInterference[i];
  delete [] m_maStrengthInterference;

}

void
CSPowerUpdater::UpdateInterference( std::vector<double>& vecPower, 
    const std::vector<int>& vecSupport,
    ClusterStructure const * const myPtrCS)
{
  std::list<std::list<int> >::const_iterator itl = myPtrCS->GetListCluMemeber().begin();
  //For each CH we will see there member 1-by-1
  for(int i=0; i<m_ptrMap->GetNumInitHeads(); i++,itl++)  {

    if(myPtrCS->GetVecHeadName().at(i) == -1 ) continue; //Cluster Not exist
    std::list<std::list<int> >::const_iterator itl2 = myPtrCS->GetListCluMemeber().begin(); //search other cluster from start for each CH
    //Search other cluster for each member in cluster (*vecHeadName)[i]
    for (int j = 0; j < m_ptrMap->GetNumInitHeads(); ++j, ++itl2) {
      int interferentSource = -1;
      double maxInterference = -99;
      //Same headIndex which means this is self head.
      if ( i == j ) {
        m_maIndexInterference[i][j]=-1;
        m_maStrengthInterference[i][j]=0;
        continue;
      }
      //empty list ->0
      else if ( itl2->size() <= 1 ) {
        m_maIndexInterference[i][j] = -1;
        m_maStrengthInterference[i][j] = 0;
        continue;
      }

      //Find the strongest interference from other group
      list<int> ::const_iterator itl3 = itl2->begin();

      //search every node(itl3) under each cluster(itl2) != self cluster
      for(; itl3!=itl2->end(); itl3++) {
        double tempInterference = vecPower.at((*itl3)) * m_ptrMap->GetGijByPair(myPtrCS->GetVecHeadName().at(i), (*itl3));
        if( myPtrCS->GetVecHeadName().at(i) == (*itl3) )
          continue;//In Downlink Scenario othe head won't transmit power
        else if ( tempInterference > maxInterference ) {
          interferentSource = (*itl3);
          maxInterference = tempInterference;
        }
      }
      m_maIndexInterference[i][j] = interferentSource;
      //Interference form j(HeadIndex) to i(HeadIndex)
      m_maStrengthInterference[i][j] = maxInterference;
    }
  }

}

void
CSPowerUpdater::ChangeAllMemberPower( std::vector<double>& vecPower, 
    std::vector<double>& vecPowerDiff, 
    std::vector<double>& vecPowerDiffRatio,
    ClusterStructure const * const myPtrCS)
{
  std::list<std::list<int> >::const_iterator it1 = myPtrCS->GetListCluMemeber().begin();
  for(int headIndex = 0; it1!=myPtrCS->GetListCluMemeber().end(); ++headIndex, ++it1) {
    //Compute the Interference headIndex received.
    if (myPtrCS->GetVecHeadName().at(headIndex) == -1) continue;

    double accuInterference = m_inBandNoise;
    for(unsigned int i=0; i<myPtrCS->GetVecHeadName().size(); i++)
    {
      if (myPtrCS->GetVecHeadName().at(i) == -1 ) continue;
      if (m_maIndexInterference[headIndex][i]!=-1)
      {
        accuInterference += m_maStrengthInterference[headIndex][i];
      }
    }
    list <int>::const_iterator it2 =it1->begin();
    // update all the power  in the cluster we are interested.
    int sizeOfCluter = it1->size();
    for(int i=0; it2!=it1->end(); i++,it2++) {
      //if the cluster has only one head or the *it2 is the head(same Name)
      if(sizeOfCluter==1||myPtrCS->GetVecHeadName().at(headIndex) == (*it2) )
      {
        vecPowerDiff[(*it2)]=0;
        vecPowerDiffRatio[(*it2)] =0;
        vecPower.at((*it2)) = 0;
        continue;
      }

      else
      {
        double powerCursor = 0;
        double tempDifference = 0;
//        cout << "m_C2:" << m_C2 << endl;
        powerCursor = (accuInterference *(pow(2.0, m_C2)-1)) / 
          m_ptrMap->GetGijByPair(myPtrCS->GetVecHeadName().at(headIndex), (*it2));
//        cout<<"Power "<<powerCursor<< " accu "<< accuInterference<<endl;
        tempDifference = powerCursor - vecPower.at((*it2));

        if(vecPowerDiff[(*it2)]!=0)vecPowerDiffRatio[(*it2)] = (double)tempDifference / vecPowerDiff[(*it2)];
        vecPowerDiff[(*it2)] = tempDifference;
        vecPower.at((*it2)) = powerCursor;
        if(vecPower.at((*it2)) > m_ptrMap->GetMaxPower()) {
          m_exceedPc = true;
        }
      }
    }
  }

}

double
CSPowerUpdater::PowerUpdate(std::vector<double>& vecPower, const std::vector<int>& vecSupport, ClusterStructure const * const myPtrCS)
{
  m_statusFlag = -1;
  int chNum = myPtrCS->GetNumHeads();
  assert(vecPower.size() == m_ptrMap->GetNumNodes());
  assert(vecSupport.size() == m_ptrMap->GetNumNodes());
  std::vector<double> vecPowerDiff(m_ptrMap->GetNumNodes());
  std::vector<double> vecPowerDiffRatio(m_ptrMap->GetNumNodes());
  std::fill(vecPower.begin(), vecPower.end(), 0.0);
  std::fill(vecPowerDiff.begin(), vecPowerDiff.end(), 0.0);
  std::fill(vecPowerDiffRatio.begin(), vecPowerDiffRatio.end(), 0.0);

  UpdateInterference( vecPower, vecSupport, myPtrCS);
  //------------------------------------//
  //add m_scale to avoid computation error//
  //------------------------------------//
  int loopCounter = 0;//counter to count how many round we updated Power
  m_exceedPc = false;
  while(loopCounter < 2 || (!m_exceedPc) ) {

    if ( loopCounter > 2 && CheckConverged(vecPowerDiffRatio, myPtrCS))
    {
      //cout<<"Ratio Converged"<<endl;
      //break;
    }
    if ( loopCounter > 2 && CheckDifference(vecPowerDiff, myPtrCS))
    {
      //cout<<"Difference is small"<<endl;
      break;
    }
    //Change the uplink node power cluster by cluster, which means that we will change all the members undet cluster i and go i+1
    // Need to update before change all member power to avoid the result from last time
    UpdateInterference( vecPower, vecSupport, myPtrCS);
    ChangeAllMemberPower(vecPower, vecPowerDiff, vecPowerDiffRatio, myPtrCS);
    loopCounter++;
  }

  if (m_avgRatio >=1)//average ratio =1 means converged already
  {
    //cout<<"diverge m_avgRatio= " <<m_avgRatio<<endl;    //PowerDiverged
    for (int i=0; i<m_ptrMap->GetNumNodes(); i++)
    {
      vecPower[i]/=m_scale;
      if (vecPower[i]>(m_ptrMap->GetMaxPower()+1e-6)){
        m_exceedPc=true;
      }
    }
    if(m_exceedPc){
      m_statusFlag=2;
    }
    else {m_statusFlag=3;}
    //cout<<"=========================="<<endl;
    return false;
  }

  for (int i=0; i<m_ptrMap->GetNumNodes(); i++) {
    if(m_avgRatio<1){
      vecPower[i]=vecPower[i] + m_avgRatio * vecPowerDiff.at(i) /(1-m_avgRatio);
      vecPower[i]/=m_scale;
    }
    if (vecPower[i]>(m_ptrMap->GetMaxPower()+1e-6)){
      m_exceedPc=true;
    }
  }
  if (m_exceedPc == true) {
    m_statusFlag=2;
    //cout<<"m_exceedPc"<<endl;    // Infeasible
    //cout<<"================="<<endl;
    return false;
  }
  //cout<<loopCounter<<endl;
  //cout<<"vvvvvvvvvvvvvvvvvvvvvFeasiblevvvvvvvvvvvvvvv"<<endl;
  m_statusFlag=1;
  return true;

}

double 
CSPowerUpdater::Solve_withT2Adj_BinerySearch_2 (double inIniT2, 
    std::vector<double>& vecPower, 
    const std::vector<int>& vecSupport, 
    ClusterStructure const * const myPtrCS){
  double powerConverge=1e-6;
  double tmpPowerMax;
  bool firstFeasibleFound=false;
  double minT2=0;
  double maxT2=inIniT2;
  int tmpNodeIndex=0;
  double tmp=(minT2+maxT2)/2;
  double minFeasibleT=DBL_MAX;
  while(1){
    m_C2=m_idtEntropy/tmp/m_ptrMap->GetBandwidth();
    bool feasFlag = PowerUpdate(vecPower, vecSupport, myPtrCS);
    if (feasFlag){
      firstFeasibleFound=true;
      returnMaxNextNodePower(tmpNodeIndex, tmpPowerMax, vecPower, myPtrCS);
      if(tmp<minFeasibleT)
        minFeasibleT=tmp;

      if (tmpNodeIndex==-1){m_best2nd_ms=0;return m_best2nd_ms;}
//      cout<<"Ori Node: "<<tmpNodeIndex<<" ori power : "<<tmpPowerMax<<endl;

      if (std::abs(m_ptrMap->GetMaxPower() - tmpPowerMax ) < powerConverge ){
//        cout << "power converge"<<endl;
//        cout << std::abs(m_ptrMap->GetMaxPower() -tmpPowerMax) << endl;
//        cout << "Time "<<(minT2+maxT2)/2<<endl;
        break;
      }
      else if((std::abs(minT2-maxT2)<1e-5)){
        //cout<<"Time converge"<<endl;
        int pi=-1;
        double pw=0;
        returnMaxNextNodePower(pi, pw, vecPower, myPtrCS);
        //cout<<"Time "<<(minT2+maxT2)/2<<endl;
        //cout<<"power max="<<pw<<endl;
        break;
      }
      maxT2=tmp;
    }
    else {
      if(firstFeasibleFound){
         if(m_statusFlag==3){
            maxT2=tmp;
         }
         else{
            minT2=tmp;
         }
      }
      else{
          maxT2*=2;
      }

      if((std::abs(minT2-maxT2)<1e-6)){
        //Feasible solution lost because of precision
        tmp  = minFeasibleT;
        m_C2 = m_idtEntropy/tmp/m_ptrMap->GetBandwidth();
        PowerUpdate(vecPower, vecSupport, myPtrCS);
        break;
      }

    }
    returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax, vecPower, myPtrCS);
//    cout<<"T2="<<tmp<<" Power = "<<tmpPowerMax<<" Feasibility="<<feasFlag<<endl;
//    cout<<"             status "<<m_statusFlag<<endl;
//    cout<<"T_u="<<maxT2<<" T_l="<<minT2<<endl;

    tmp=(minT2+maxT2)/2;
  }
  //showVerification
  m_criticalNode      = tmpNodeIndex;
  m_criticalHeadIndex = returnHeadIndex_ByNodeName(m_criticalNode, myPtrCS);
  m_best2nd_ms        = tmp;
  returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax, vecPower, myPtrCS);

  return m_best2nd_ms;
}

void 
CSPowerUpdater::returnMaxNextNodePower(int &nodeIndex, double &nodePower, std::vector<double>& vecPower,ClusterStructure const * const myPtrCS){
  nodePower = 0;
  nodeIndex = -1;
  for(int i=0; i < m_ptrMap->GetNumNodes(); ++i){
    if (vecPower.at(i) > nodePower){
      nodePower = vecPower.at(i);
      nodeIndex = i;
    }
  }
  //assert(nodeIndex!=-1);
}

int 
CSPowerUpdater::returnHeadIndex_ByNodeName(int nodeName, ClusterStructure const * const myPtrCS){
    list<list<int> >::const_iterator it_LiInt = myPtrCS->GetListCluMemeber().begin();
    for(int i=0; it_LiInt != myPtrCS->GetListCluMemeber().end(); ++it_LiInt, ++i){
        list<int>::const_iterator it_Int = it_LiInt->begin();
        for(;it_Int != it_LiInt->end();it_Int++){
            if(*it_Int == nodeName){
                return i;
            }
        }
    }
    cout<<"Head Name Not Found in ULConstraintSolver!"<<endl;
    return -1;
}

bool
CSPowerUpdater::CheckDifference(const std::vector<double>& vecPowerDiff, ClusterStructure const * const myPtrCS ) const
{
  bool allTooSmall = true;
  for(int i = 0; i < m_ptrMap->GetNumNodes(); i++) {
    if( (std::abs(vecPowerDiff.at(i)) > m_threshold) && 
        (vecPowerDiff.at(i) != 0) ) {
      allTooSmall = false;
      break;
    }
  }
  return allTooSmall;
}

bool
CSPowerUpdater::CheckConverged(const std::vector<double>& vecPowerDiffRatio, ClusterStructure const * const myPtrCS) 
{
  bool convergence = true;
  m_avgRatio = 0; // find the m_avgRatio
  double avgMember =0;
  int chNum = myPtrCS->GetNumHeads();

  //compute average power ratio
  for (int i =0; i<m_ptrMap->GetNumNodes(); i++) {
    if (vecPowerDiffRatio.at(i)!=0) {
      m_avgRatio += vecPowerDiffRatio.at(i);
      avgMember++;
    }
  }
//  std::cout << avgMember << std::endl;
//  cout<<"AvGr" << m_avgRatio <<endl;
  if(avgMember!=0) m_avgRatio/=avgMember;
  //check if all the power ratio close to the average

  for (int i = 0; i < chNum; i++)
  {
    assert(m_avgRatio>-1);
    if ( (std::abs(vecPowerDiffRatio.at(i)-m_avgRatio) > m_threshold) && 
        (vecPowerDiffRatio.at(i)!=0) )
    {
      convergence = false;
      break;
    }
  }
  return convergence;
}


void 
CSPowerUpdater::showVerificationResult(
    std::vector<double>& vecPower,
    const std::vector<int>& vecSupport,
    ClusterStructure const * const myPtrCS,
    std::vector<double>& vecBestReceivedInterference,
    std::vector<double>& vecBestSINR_forVerification, 
    std::vector<double>& vecBestBpsHz_forVerification
    )
{
    vecBestReceivedInterference.clear();
    UpdateInterference(vecPower, vecSupport, myPtrCS);
    list <list <int> >::const_iterator it1 = myPtrCS->GetListCluMemeber().begin();
    for(int headIndex=0; it1 != myPtrCS->GetListCluMemeber().end(); headIndex++,it1++)
    {
        //Compute the Interference headIndex received.
        if (myPtrCS->GetVecHeadName().at(headIndex)==-1)continue;
        double accuInterference = m_inBandNoise;
        for(int i=0; i< m_ptrMap->GetNumInitHeads(); i++)
        {
            if(m_maIndexInterference[headIndex][i]!=-1)
            {
                accuInterference+= m_maStrengthInterference[headIndex][i];
            }
        }
        vecBestReceivedInterference.push_back(accuInterference - m_inBandNoise);
        //cout<<"Cluster: "<<headIndex<<" BpsHz "<<C2<<endl;
        list <int>::const_iterator it2 =it1->begin();
        // update all the power  in the cluster we are interested.
        int sizeOfCluter = it1->size();
        for(int i=0; it2!=it1->end(); i++,it2++)
        {
            vecBestSINR_forVerification[*it2] = m_ptrMap->GetGijByPair(myPtrCS->GetVecHeadName().at(headIndex), (*it2))* vecPower.at(*it2) /accuInterference;
            vecBestBpsHz_forVerification[*it2]=(log2(1+ m_ptrMap->GetGijByPair(myPtrCS->GetVecHeadName().at(headIndex), (*it2) ) * vecPower.at(*it2) /accuInterference)/(sizeOfCluter-1));
            /*
            cout.precision(5);
            cout<<scientific<<"Normal: Node "<<*it2<<" SINR "<< Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)] /accuInterference<< \
            "  BpsHz:"<<log2(1+Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)] /accuInterference)/(sizeOfCluter-1)<<endl;
            */
        }
//        cout<<"---------------"<<endl;
    }
}
