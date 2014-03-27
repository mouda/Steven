#include "bruteForceScheduler.h"

BruteForceScheduler::BruteForceScheduler( const double txTime, 
    const double bandwidthKhz, 
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrMatComputer, 
    ClusterStructure const * const ptrCS): 
    m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),
    m_type("BruteForceScheduler")
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads =  m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();

  /* initialize the fixed power */
  m_vecNodePower.resize(m_numNodes);
  fill(m_vecNodePower.begin(), m_vecNodePower.end(), m_maxPower);

  Init();
  /* initialize the scheduled record */
  m_vecSched.resize(m_numNodes);
  fill(m_vecSched.begin(), m_vecSched.end(), 0);

}

BruteForceScheduler::~BruteForceScheduler()
{

}

double 
BruteForceScheduler::ScheduleOneSlot(std::vector<int>& vecSupport )
{
  fill(vecSupport.begin(), vecSupport.end(), 1);
  return 0.0;
}

bool 
BruteForceScheduler::ScheduleOneSlot(std::vector<int>& vecSupport, const std::vector<double>& vecVariance)
{
  fill(vecSupport.begin(), vecSupport.end(), 1);
  return true;
}

void
BruteForceScheduler::Init()
{
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

  Eigen::MatrixXd tmp = m_A+m_B-m_C;
  
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

double
BruteForceScheduler::OmegaValue(const int nodeName)
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
