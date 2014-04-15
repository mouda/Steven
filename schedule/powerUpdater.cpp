#include "powerUpdater.h"
#include <cmath>
#include <cassert>

PowerUpdater::PowerUpdater(
    Map const * const ptrMap, 
    ClusterStructure const * const ptrCS
    ):
  m_threshold(1e-12),
  m_ptrMap(ptrMap), 
  m_ptrCS(ptrCS) 
{
  m_maIndexInterference = new int* [m_ptrMap->GetNumInitHeads()];
  for (int i=0; i < m_ptrMap->GetNumInitHeads(); i++) 
    m_maIndexInterference[i] = new int [m_ptrMap->GetNumInitHeads()];

  m_maStrengthInterference = new double* [m_ptrMap->GetNumInitHeads()];
  for (int i=0; i < m_ptrMap->GetNumInitHeads(); i++) 
    m_maStrengthInterference[i] = new double [m_ptrMap->GetNumInitHeads()];
}

PowerUpdater::~PowerUpdater()
{
  for(int i = 0; i < m_ptrMap->GetNumInitHeads(); i++) delete[] m_maIndexInterference[i];
  delete [] m_maIndexInterference;
  for(int i = 0; i < m_ptrMap->GetNumInitHeads(); i++) delete[] m_maStrengthInterference[i];
  delete [] m_maStrengthInterference;

}

void
PowerUpdater::UpdateInterference( std::vector<double>& vecPower)
{
  std::list<std::list<int> >::const_iterator itl = m_ptrCS->GetListCluMemeber().begin();
  //For each CH we will see there member 1-by-1
  for(int i=0; i<m_ptrMap->GetNumInitHeads(); i++,itl++)  {

    if(m_ptrCS->GetVecHeadName().at(i) == -1 ) continue; //Cluster Not exist
    std::list<std::list<int> >::const_iterator itl2 = m_ptrCS->GetListCluMemeber().begin(); //search other cluster from start for each CH
    //Search other cluster for each member in cluster (*vecHeadName)[i]
    for (int j = 0; j < m_ptrMap->GetNumInitHeads(); j++, itl2++) {
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
        double tempInterference = vecPower.at((*itl3)) * m_ptrMap->GetGijByPair(m_ptrCS->GetVecHeadName().at(i), (*itl3));
        if( m_ptrCS->GetVecHeadName().at(i) == (*itl3) )
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
PowerUpdater::ChangeAllMemberPower( std::vector<double>& vecPower) const
{
  std::list<std::list<int> >::const_iterator it1 = m_ptrCS->GetListCluMemeber().begin();
  for(int headIndex = 0; it1!=m_ptrCS->GetListCluMemeber().end(); ++headIndex, ++it1) {
    //Compute the Interference headIndex received.
    if ((*vecHeadName)[headIndex]==-1)continue;

    double accuInterference = inBandNoise;
    for(unsigned int i=0; i<(*vecHeadName).size(); i++)
    {
      if ((*vecHeadName)[i]==-1)continue;
      if(maIndexInterference[headIndex][i]!=-1)
      {
        accuInterference+= maStrengthInterference[headIndex][i];
      }
    }
    //cout<<"Cluster: "<<headIndex<<" ->"<<accuInterference<<" C2 "<<C2<<endl;
    list <int>::const_iterator it2 =it1->begin();
    // update all the power  in the cluster we are interested.
    int sizeOfCluter = it1->size();
    for(int i=0; it2!=it1->end(); i++,it2++) {
      if(sizeOfCluter==1||(*vecHeadName)[headIndex]==(*it2))//if the cluster has only one head or the *it2 is the head(same Name)
      {
        powerDifference[(*it2)]=0;
        powerDifferenceRatio[(*it2)] =0;
        vecPower.at((*it2)) = 0;
        continue;
      }

      else
      {
        float powerCursor = 0;
        double tempDifference = 0;
        //cout<<"c2 "<<C2<<endl;
        //cout<<"Channel Gain from"<<(*vecHeadName)[headIndex]<<" "<<(*it2)<<" "<<Gij[(*vecHeadName)[headIndex]][(*it2)];
        powerCursor = (accuInterference *(pow((float)2.0, (sizeOfCluter-1)*C2)-1))/Gij[(*vecHeadName)[headIndex]][(*it2)];
        //cout<<"Power "<<powerCursor/scale<< " accu "<< accuInterference/scale<<endl;
        tempDifference = powerCursor - vecPower.at((*it2));

        if(powerDifference[(*it2)]!=0)powerDifferenceRatio[(*it2)] = (double)tempDifference / powerDifference[(*it2)];
        powerDifference[(*it2)] = tempDifference;
        vecPower.at((*it2)) = powerCursor;
        if(vecPower.at((*it2)) > powerMax) {
          exceedPc = true;
        }
      }
    }
  }

}

double
PowerUpdater::Solve(const std::vector<int>& vecSupport)
{
  int statusFlag = -1;
  int chNum = m_ptrCS->GetNumHeads();
  std::vector<double> vecPower(m_ptrMap->GetNumNodes());
  std::vector<double> vecPowerDiff(m_ptrMap->GetNumNodes());
  std::vector<double> vecPowerDiffRatio(m_ptrMap->GetNumNodes());
  std::fill(vecPower.begin(), vecPower.end(), 0.0);
  std::fill(vecPowerDiff.begin(), vecPowerDiff.end(), 0.0);
  std::fill(vecPowerDiffRatio.begin(), vecPowerDiffRatio.end(), 0.0);

  UpdateInterference( vecPower);
  //------------------------------------//
  //add m_scale to avoid computation error//
  //------------------------------------//
  int loopCounter = 0;//counter to count how many round we updated Power
  bool exceedPc = false;
  while(loopCounter <2||(!exceedPc)) {
    //cout<<"loop "<<loopCounter<<endl;
    //cout<<"avgR"<< m_avgRatio<<endl;
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
    UpdateInterference( vecPower);// Need to update before change all member power to avoid the result from last time
    ChangeAllMemberPower(vecPower);
    loopCounter++;
    //cout<<loopCounter<<"-th round"<<endl;

  }

  if (m_avgRatio >=1)//average ratio =1 means converged already
  {
    //cout<<"diverge m_avgRatio= " <<m_avgRatio<<endl;    //PowerDiverged
    for (int i=0; i<m_ptrMap->GetNumNodes(); i++)
    {
      vecPower[i]/=m_scale;
      if (vecPower[i]>(m_ptrMap->GetMaxPower()+1e-6)){// cout<<vecPower[i]<<endl;
        exceedPc=true;
      }
      //cout<<"node "<<i <<": "<<vecPower[i]<<endl;
      //cout<<"  =>"<<vecPower[i]<<endl;
    }
    if(exceedPc){
      statusFlag=2;
    }
    else {statusFlag=3;}
    //cout<<"=========================="<<endl;
    return false;
  }

  for (int i=0; i<m_ptrMap->GetNumNodes(); i++)
  {
    //cout<<"node i :"<<i <<" Ori= "<<vecPower[i];
    if(m_avgRatio<1){
      vecPower[i]=vecPower[i] + m_avgRatio * vecPowerDiff.at(i) /(1-m_avgRatio);
      vecPower[i]/=m_scale;
    }
    //cout<<"node "<<i <<": "<<vecPower[i]<<endl;
    if (vecPower[i]>(m_ptrMap->GetMaxPower()+1e-6)){// cout<<vecPower[i]<<endl;
      exceedPc=true;}

      //cout<<"  =>"<<vecPower[i]<<endl;
  }
  //showVerificationResult();
  if (exceedPc == true)
  {
    statusFlag=2;
    //cout<<"exceedPc"<<endl;    // Infeasible
    //cout<<"================="<<endl;
    return false;
  }
  //cout<<loopCounter<<endl;
  //cout<<"vvvvvvvvvvvvvvvvvvvvvFeasiblevvvvvvvvvvvvvvv"<<endl;
  statusFlag=1;
  return true;

}

bool
PowerUpdater::CheckDifference(const std::vector<double>& vecPowerDiff ) const
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
PowerUpdater::CheckConverged(const std::vector<double>& vecPowerDiffRatio) 
{
  bool convergence = true;
  double thre = 1e-5;
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
  if(avgMember!=0) m_avgRatio/=avgMember;
  //cout<<"AvGr"<<m_avgRatio<<endl;
  //check if all the power ratio close to the average

  for (int i =0; i<chNum; i++)
  {
    assert(m_avgRatio>-1);
    if (abs(vecPowerDiffRatio.at(i)-m_avgRatio)>thre&&vecPowerDiffRatio.at(i)!=0)
    {
      convergence = false;
      break;
    }
  }
  return convergence;
}
