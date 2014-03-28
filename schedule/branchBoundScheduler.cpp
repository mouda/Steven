#include "branchBoundScheduler.h"

using Ipopt::Index;
BranchBoundScheduler::BranchBoundScheduler( const double txTime, 
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("BranchBound")
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads = m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();
  Eigen::IOFormat CleanFmt(2, 0, ", ", "\n", "[", "]");
  m_A = Eigen::MatrixXd::Zero(m_numMaxHeads, m_numNodes); 
  m_B = Eigen::MatrixXd::Zero(m_numMaxHeads, m_numNodes); 
  m_C = Eigen::MatrixXd::Zero(m_numMaxHeads, m_numNodes);
  m_X = Eigen::MatrixXd::Zero(m_numMaxHeads, m_numNodes);
  for (int i = 0; i < m_numMaxHeads; ++i) {
    for (int j = 0; j < m_numNodes; ++j) {
     if (i == m_ptrCS->GetChIdxByName(j)) {
       m_A(i,j) = OmegaValue(j);
     }  
    }
  }
  double exponent = pow(2.0,m_ptrMap->GetIdtEntropy()/m_bandwidthKhz/m_txTimePerSlot);
  for (int i = 0; i < m_numMaxHeads; ++i) {
    for (int j = 0; j < m_numNodes; ++j) {
     if (i != m_ptrCS->GetChIdxByName(j)) {
       m_B(i,j) = (exponent-1) * m_maxPower * 
         m_ptrMap->GetGijByPair(m_ptrCS->GetVecHeadName()[i],j);
     }  
    }
  }
  for (int i = 0; i < m_numMaxHeads; ++i) {
    for (int j = 0; j < m_numNodes; ++j) {
      if (i == m_ptrCS->GetChIdxByName(j)) {
        int headName = m_ptrCS->GetVecHeadName()[i];
        m_C(i,j) = m_maxPower * m_ptrMap->GetGijByPair(headName,j) - 
          (exponent -1)*m_bandwidthKhz * m_ptrMap->GetNoise() - OmegaValue(j) ;
      }
    }
  }

  Eigen::MatrixXd tmp = m_A+m_B-m_C;
 
#ifdef DEBUG 
  cout << "================ m_A ================" << endl;
  cout << m_A.format(CleanFmt) << endl;
  cout << "================ m_B ================" << endl;
  cout << m_B.format(CleanFmt) << endl;
  cout << "================ m_C ================" << endl;
  cout << m_C.format(CleanFmt) << endl;
  cout << "================ m_A + m_B-m_C ================" << endl;
  cout << tmp.format(CleanFmt) << endl;
  cout << "Varince: " << m_ptrMatComputer->GetCorrationFactor() << endl;
#endif

  
  double corrFactor = m_ptrMatComputer->GetCorrationFactor();
  double variance = 1.0;
  m_Signma = Eigen::MatrixXd::Zero(m_numNodes, m_numNodes);
  for (int i = 0; i < m_numNodes; ++i) {
    for (int j = 0; j < m_numNodes; ++j) {
      m_Signma(i,j) = variance * exp(-1*(m_ptrMatComputer->GetDijSQByPair(i,j))/corrFactor) ;
    }
  }

  m_vecSched.resize(m_numNodes);
  fill(m_vecSched.begin(), m_vecSched.end(), 0);
}

BranchBoundScheduler::~BranchBoundScheduler()
{

}

double
BranchBoundScheduler::ScheduleOneSlot( std::vector<int>& vecSupport )
{
  Eigen::MatrixXd tmp = m_Signma;
  int activeNodes = this->SolverHook(vecSupport, tmp);
  double result = activeNodes * m_ptrMap->GetIdtEntropy() + m_ptrMatComputer->computeLog2Det(1.0, vecSupport);
  return result;
}

bool
BranchBoundScheduler::ScheduleOneSlot( std::vector<int>& vecSupport, const std::vector<double>& vecVariance)
{
  /* construct the covariance matrix */
  Eigen::MatrixXd mySigma = Eigen::MatrixXd::Zero(m_numNodes, m_numNodes);
  for (int i = 0; i < m_numNodes; ++i) {
    for (int j = 0; j < m_numNodes; ++j) {
      mySigma(i,j) = vecVariance.at(i) * vecVariance.at(j) * exp(-1*(m_ptrMatComputer->GetDijSQByPair(i,j))/m_ptrMatComputer->GetSpatialCorrelationFactor()) ;
    }
  }
  int activeNodes = this->SolverHook(vecSupport, mySigma);
  return true;
}

int
BranchBoundScheduler::SolverHook(std::vector<int>& vecSupport, Eigen::MatrixXd& matSigma)
{
  Eigen::IOFormat CleanFmt(2, 0, " ", "\n", "", ";");
  Index numVariables = m_numNodes;
  Index numConstraints = 2 * m_numMaxHeads;
  Index numNz_jac_g = 2 * m_numMaxHeads * m_numNodes;
  Index numNz_h_lag = 0;
  
  SmartPtr<TMINLP> tminlp = 
    new MyTMINLP( numVariables, numConstraints, numNz_jac_g, numNz_h_lag, matSigma,
        (m_A+m_B-m_C), m_ptrCS, m_ptrMap);
  MyTMINLP* rawPtr = dynamic_cast<MyTMINLP*>(GetRawPtr(tminlp));
  std::vector<int> tmp(m_ptrMap->GetNumNodes(), 0);
  rawPtr->SetExtraConstraints(tmp);
  //rawPtr->SetExtraConstraints(m_vecSched);

  FILE * fp = fopen("log.out","w");
  CoinMessageHandler handler(fp);
  BonminSetup bonmin(&handler);
  bonmin.initializeOptionsAndJournalist();
  // Here we can change the default value of some Bonmin or Ipopt option
  bonmin.options()->SetNumericValue("bonmin.time_limit", 600); //changes bonmin's time limit
  bonmin.options()->SetStringValue("mu_oracle","loqo");
  //Here we read several option files
  bonmin.readOptionsFile("Mybonmin.opt");
  bonmin.readOptionsFile();// This reads the default file "bonmin.opt"
  // Options can also be set by using a string with a format similar to the bonmin.opt file
  bonmin.readOptionsString("bonmin.algorithm B-BB\n");

  //Now initialize from tminlp
  bonmin.initialize(GetRawPtr(tminlp));
  //Set up done, now let's branch and bound
  try {
    Bab bb;
    bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc
  }
  catch(TNLPSolver::UnsolvedError *E) {
    //There has been a failure to solve a problem with Ipopt.
    std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
  }
  catch(OsiTMINLPInterface::SimpleError &E) {
    std::cerr << E.className() << "::" << E.methodName()
	      << std::endl
	      << E.message() << std::endl;
  }
  catch(CoinError &E) {
    std::cerr << E.className() << "::" << E.methodName()
	      << std::endl
	      << E.message() << std::endl;
  }

  vecSupport.assign(rawPtr->GetVecSolution().begin(), rawPtr->GetVecSolution().end());
  for (int i = 0; i < m_vecSched.size(); ++i) {
    m_vecSched.at(i) = m_vecSched.at(i) + vecSupport.at(i);
  }

  int activeNodes = 0;
  for (int i = 0; i < m_numNodes; i++) {
    if (vecSupport[i] == 1) ++activeNodes;
  }
  return activeNodes;
}

double
BranchBoundScheduler::OmegaValue(const int nodeName)
{
  Eigen::MatrixXd matOwnership = Eigen::MatrixXd::Zero(m_numMaxHeads, m_numNodes);
  list<list<int> >::const_iterator iterRow = m_ptrCS->GetListCluMemeber().begin();
  for (int i = 0; i < m_numMaxHeads; ++i, ++iterRow) {
    list<int>::const_iterator iterCol= iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      matOwnership(i,*iterCol) = 1;
    }
  }
  Eigen::MatrixXd matGij(m_numNodes, m_numNodes);
  for (int i = 0; i < m_numNodes; ++i) {
    for (int j = 0; j < m_numNodes; ++j) {
      matGij(i,j) = m_ptrMap->GetGijByPair(i,j);
    }
  }
  double interference = 0.0;
  int index = 0;
  double rate = pow(2.0,m_ptrMap->GetIdtEntropy());
  for (int i = 0; i < m_numMaxHeads; i++) {
    if ( m_ptrCS->GetChIdxByName(nodeName) == i ) continue; 
    Eigen::Matrix<double,Eigen::Dynamic,1> vecTmp(m_numNodes,1);
    int headName =m_ptrCS->GetVecHeadName()[i];
    vecTmp = matOwnership.row(i).cwiseProduct(matGij.row(nodeName));
    double maxValue = vecTmp.maxCoeff(&index); 
    interference += maxValue * m_maxPower;
  }
  double lhs = (pow(2.0,m_ptrMap->GetIdtEntropy()/m_bandwidthKhz/m_txTimePerSlot) - 1 ) * 
    (m_ptrMap->GetNoise() + interference);
  double rhs = m_maxPower * 
    matGij(nodeName, m_ptrCS->GetChNameByName(nodeName));

  if (lhs > rhs ) {
    return -1*lhs;
  }
  else {
    return 0.0;
  }
  return 0.0;
}
