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

}

PowerUpdater::~PowerUpdater()
{

}

void
PowerUpdater::UpdateInterference( std::vector<double>& vecPower)
{

}

void
PowerUpdater::ChangeAllMemberPower( std::vector<double>& vecPower) const
{

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
