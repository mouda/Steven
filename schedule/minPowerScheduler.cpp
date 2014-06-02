#include "minPowerScheduler.h"

using Ipopt::Index;

MinPowerScheduler::MinPowerScheduler( const double txTime, 
    const int tier2NumSlot,
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), 
    m_tier2NumSlot(tier2NumSlot),
    m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("MinPower")
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numHeads =  m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();
  m_vecNodePower.resize(m_numNodes);
  fill(m_vecNodePower.begin(), m_vecNodePower.end(), m_maxPower);

//  m_matA = Eigen::MatrixXd::Zero(m_numNodes, m_numNodes);
//  m_matB = Eigen::MatrixXd::Zero(m_numNodes, m_tier2NumSlot);
//  m_matO = Eigen::MatrixXd::Zero(m_numHeads, 1);
//  m_matE = Eigen::MatrixXd::Zero(m_numNodes, m_numHeads);

  int numConstraints = m_numNodes * m_tier2NumSlot + m_numNodes + m_numHeads + m_numHeads * m_tier2NumSlot;
  int numVariables = 2 * m_numNodes * m_tier2NumSlot;

  InitSolution();
}

MinPowerScheduler::~MinPowerScheduler()
{

}

double
MinPowerScheduler::ScheduleOneSlot( std::vector<int>& vecSupport )
{
  /* NULL */
  return 0.0;
}

double
MinPowerScheduler::ScheduleOneSlot( std::vector<int>& vecSupport, std::vector<double>& vecPower, const std::vector<double>& vecVariance)
{

  return 0.0;
}


void 
MinPowerScheduler::InitSolution()
{

  Eigen::IOFormat CleanFmt(2, 0, " ", "\n", "", ";");
  Index numVariables = 2* m_numNodes * m_tier2NumSlot;
  Index numConstraints = m_numNodes * m_tier2NumSlot + m_numNodes + m_numHeads + m_numHeads * m_tier2NumSlot;
  Index numNz_jac_g = m_numNodes * m_tier2NumSlot + m_numNodes + m_numHeads + m_numHeads * m_tier2NumSlot;
  Index numNz_h_lag = 0;
  
  SmartPtr<TMINLP> tminlp = 
    new MinPowerMILP( numVariables, numConstraints, numNz_jac_g, numNz_h_lag,
        m_txTimePerSlot, m_tier2NumSlot, m_bandwidthKhz,  
        m_ptrCS, m_ptrMap);
//  MinPowerMILP* rawPtr = dynamic_cast<MinPowerMILP*>(GetRawPtr(tminlp));
//  std::vector<int> tmp(m_ptrMap->GetNumNodes(), 0);
//  rawPtr->SetExtraConstraints(tmp);
//  //rawPtr->SetExtraConstraints(m_vecSched);

//  FILE * fp = fopen("log.out","w");
//  CoinMessageHandler handler(fp);
//  BonminSetup bonmin(&handler);
//  bonmin.initializeOptionsAndJournalist();
//  // Here we can change the default value of some Bonmin or Ipopt option
//  bonmin.options()->SetNumericValue("bonmin.time_limit", 1000); //changes bonmin's time limit
//  bonmin.options()->SetStringValue("mu_oracle","loqo");
//  //Here we read several option files
//  bonmin.readOptionsFile("Mybonmin.opt");
//  bonmin.readOptionsFile();// This reads the default file "bonmin.opt"
//  // Options can also be set by using a string with a format similar to the bonmin.opt file
//  bonmin.readOptionsString("bonmin.algorithm B-BB\n");
//
//  //Now initialize from tminlp
//  bonmin.initialize(GetRawPtr(tminlp));
//  //Set up done, now let's branch and bound
//  try {
//    Bab bb;
//    bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc
//  }
//  catch(TNLPSolver::UnsolvedError *E) {
//    //There has been a failure to solve a problem with Ipopt.
//    std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
//  }
//  catch(OsiTMINLPInterface::SimpleError &E) {
//    std::cerr << E.className() << "::" << E.methodName()
//	      << std::endl
//	      << E.message() << std::endl;
//  }
//  catch(CoinError &E) {
//    std::cerr << E.className() << "::" << E.methodName()
//	      << std::endl
//	      << E.message() << std::endl;
//  }
//
//  m_matSolution = rawPtr->GetMatSolution();
}
