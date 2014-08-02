#include "newImagePowerUpdater.h"
#include <cmath>
#include <cassert>
#include <float.h>

NewImagePowerUpdater::NewImagePowerUpdater(
    ImageMap const * const ptrMap, 
    ClusterStructure const * const ptrCS,
    const double txTimeSlot,
    const double txNumSlot
    ):
  m_threshold(1e-12),
  m_ptrMap(ptrMap), 
  m_ptrCS(ptrCS), 
  m_idtEntropy(ptrMap->GetIdtEntropy()),
  m_txTimePerSlot(txTimeSlot)
{
  m_maIndexInterference = new int* [m_ptrMap->GetNumInitHeads()];
  for (int i=0; i < m_ptrMap->GetNumInitHeads(); i++) 
    m_maIndexInterference[i] = new int [m_ptrMap->GetNumInitHeads()];

  m_maStrengthInterference = new double* [m_ptrMap->GetNumInitHeads()];
  for (int i=0; i < m_ptrMap->GetNumInitHeads(); i++) 
    m_maStrengthInterference[i] = new double [m_ptrMap->GetNumInitHeads()];
  m_inBandNoise = m_ptrMap->GetNoise() * m_scale;
  m_C2=m_idtEntropy/m_txTimePerSlot/m_ptrMap->GetBandwidth();
  m_vecC2.resize(m_ptrMap->GetNumNodes());
  double totalTier2Time = m_txTimePerSlot*txNumSlot; 
//  for (int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
//    m_vecC2.at(i) = m_ptrMap->GetIdtEntropy(i)*m_ptrCS->GetChIdxByNamem_ptrCS->GetChIdxByName(i)/totalTier2Time/m_ptrMap->GetBandwidth();
//    cout << m_ptrMap->GetIdtEntropy(i) << endl;
//  }

  std::list<std::list<int> >::const_iterator iterRow = m_ptrCS->GetListCluMemeber().begin();
  for (; iterRow != m_ptrCS->GetListCluMemeber().end(); ++iterRow) {
    std::list<int>::const_iterator iterCol = iterRow->begin();
    double txTime = totalTier2Time/(static_cast<double>(iterRow->size())-1);
    for (; iterCol != iterRow->end(); ++iterCol) {
      m_vecC2.at(*iterCol) = m_ptrMap->GetIdtEntropy(*iterCol)/txTime/m_ptrMap->GetBandwidth();
    }
  }
  cout << "********idtEntropy " << m_idtEntropy << endl;
  cout << "********m_txTimePerSlot " << m_txTimePerSlot << endl;
  cout << "********bandwidth " << m_ptrMap->GetBandwidth() << endl;
}

NewImagePowerUpdater::~NewImagePowerUpdater()
{
  for(int i = 0; i < m_ptrMap->GetNumInitHeads(); i++) delete[] m_maIndexInterference[i];
  delete [] m_maIndexInterference;
  for(int i = 0; i < m_ptrMap->GetNumInitHeads(); i++) delete[] m_maStrengthInterference[i];
  delete [] m_maStrengthInterference;

}

void
NewImagePowerUpdater::UpdateInterference( std::vector<double>& vecPower, const std::vector<int>& vecSupport)
{
  std::list<std::list<int> >::const_iterator itl = m_ptrCS->GetListCluMemeber().begin();
  //For each CH we will see there member 1-by-1
  for(int i=0; i<m_ptrMap->GetNumInitHeads(); i++,itl++)  {

    if(m_ptrCS->GetVecHeadName().at(i) == -1 ) continue; //Cluster Not exist
    std::list<std::list<int> >::const_iterator itl2 = m_ptrCS->GetListCluMemeber().begin(); //search other cluster from start for each CH
    //Search other cluster for each member in cluster (*vecHeadName)[i]
    for (int j = 0; j < m_ptrMap->GetNumInitHeads(); ++j, ++itl2) {
      int interferentSource = -1;
      double maxInterference = -DBL_MAX;
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
        double tempInterference = vecPower.at((*itl3)) * m_ptrMap->GetGijByPair(m_ptrCS->GetVecHeadName().at(i), (*itl3));
        if( m_ptrCS->GetVecHeadName().at(i) == (*itl3) )
          continue;//In Downlink Scenario othe head won't transmit power
        else if (   tempInterference > maxInterference ) {
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
NewImagePowerUpdater::ChangeAllMemberPower( std::vector<double>& vecPower, 
    std::vector<double>& vecPowerDiff, 
    std::vector<double>& vecPowerDiffRatio)
{
  std::list<std::list<int> >::const_iterator it1 = m_ptrCS->GetListCluMemeber().begin();
  for(int headIndex = 0; it1!=m_ptrCS->GetListCluMemeber().end(); ++headIndex, ++it1) {
    //Compute the Interference headIndex received.
    if (m_ptrCS->GetVecHeadName().at(headIndex) == -1) continue;
    double accuInterference = m_inBandNoise;
    for(unsigned int i=0; i<m_ptrCS->GetVecHeadName().size(); i++) {
      if (m_ptrCS->GetVecHeadName().at(i) == -1 ) continue;
      if (m_maIndexInterference[headIndex][i]!=-1) {
        accuInterference += m_maStrengthInterference[headIndex][i];
      }
    }
    list <int>::const_iterator it2 =it1->begin();
    // update all the power  in the cluster we are interested.
    int sizeOfCluter = it1->size();
    for(int i=0; it2!=it1->end(); i++,it2++) {
      //if the cluster has only one head or the *it2 is the head(same Name)
      if(sizeOfCluter==1||m_ptrCS->GetVecHeadName().at(headIndex) == (*it2) ) {
        vecPowerDiff[(*it2)]=0;
        vecPowerDiffRatio[(*it2)] =0;
        vecPower.at((*it2)) = 0;
        continue;
      }
      else {
        double powerCursor = 0;
        double tempDifference = 0;
//        cout << "C2:" << m_C2 << endl;
        powerCursor = (accuInterference *(pow(2.0, m_vecC2.at(*it2))-1)) / 
          m_ptrMap->GetGijByPair(m_ptrCS->GetVecHeadName().at(headIndex), (*it2));
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
NewImagePowerUpdater::Solve(std::vector<double>& vecPower, const std::vector<int>& vecSupport)
{
  int statusFlag = -1;
  int chNum = m_ptrCS->GetNumHeads();
  assert(vecPower.size() == m_ptrMap->GetNumNodes());
  assert(vecSupport.size() == m_ptrMap->GetNumNodes());
  std::vector<double> vecPowerDiff(m_ptrMap->GetNumNodes());
  std::vector<double> vecPowerDiffRatio(m_ptrMap->GetNumNodes());
  std::fill(vecPower.begin(), vecPower.end(), 0.0);
  std::fill(vecPowerDiff.begin(), vecPowerDiff.end(), 0.0);
  std::fill(vecPowerDiffRatio.begin(), vecPowerDiffRatio.end(), 0.0);

  UpdateInterference( vecPower, vecSupport);
  //------------------------------------//
  //add m_scale to avoid computation error//
  //------------------------------------//
  int loopCounter = 0;//counter to count how many round we updated Power
  bool m_exceedPc = false;
  while(loopCounter <2||(!m_exceedPc)) {

    if ( loopCounter > 2 && CheckConverged(vecPowerDiffRatio))
    {
      //cout<<"Ratio Converged"<<endl;
      //break;
    }
    if ( loopCounter > 2 && CheckDifference(vecPowerDiff))
    {
      //cout<<"Difference is small"<<endl;
      break;
    }
    //Change the uplink node power cluster by cluster, which means that we will change all the members undet cluster i and go i+1
    // Need to update before change all member power to avoid the result from last time
    UpdateInterference( vecPower, vecSupport);
    ChangeAllMemberPower(vecPower, vecPowerDiff, vecPowerDiffRatio);
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
      statusFlag=2;
    }
    else {statusFlag=3;}
    //cout<<"=========================="<<endl;
    return false;
  }

  for (int i=0; i<m_ptrMap->GetNumNodes(); i++)
  {
    if(m_avgRatio<1){
      vecPower[i]=vecPower[i] + m_avgRatio * vecPowerDiff.at(i) /(1-m_avgRatio);
      vecPower[i]/=m_scale;
    }
    if (vecPower[i]>(m_ptrMap->GetMaxPower()+1e-6)){
      m_exceedPc=true;
    }
  }
  //showVerificationResult();
  if (m_exceedPc == true)
  {
    statusFlag=2;
    //cout<<"m_exceedPc"<<endl;    // Infeasible
    //cout<<"================="<<endl;
    return false;
  }
  //cout<<loopCounter<<endl;
  //cout<<"vvvvvvvvvvvvvvvvvvvvvFeasiblevvvvvvvvvvvvvvv"<<endl;
  statusFlag=1;
    return true;

  }

bool
NewImagePowerUpdater::CheckDifference(const std::vector<double>& vecPowerDiff ) const
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
NewImagePowerUpdater::CheckConverged(const std::vector<double>& vecPowerDiffRatio) 
{
  bool convergence = true;
  m_avgRatio = 0; // find the m_avgRatio
  double avgMember =0;
  int chNum = m_ptrCS->GetNumHeads();

  //compute average power ratio
  for (int i =0; i<m_ptrMap->GetNumNodes(); i++) {
    if (vecPowerDiffRatio.at(i)!=0) {
      m_avgRatio += vecPowerDiffRatio.at(i);
      avgMember++;
    }
  }
//  std::cout << avgMember << std::endl;
  if(avgMember!=0) m_avgRatio/=avgMember;
  //cout<<"AvGr"<<m_avgRatio<<endl;
  //check if all the power ratio close to the average

  for (int i =0; i<chNum; i++)
  {
    assert(m_avgRatio>-1);
    if ( (abs(vecPowerDiffRatio.at(i)-m_avgRatio) > m_threshold) && 
        (vecPowerDiffRatio.at(i)!=0) )
    {
      convergence = false;
      break;
    }
  }
  return convergence;
}
