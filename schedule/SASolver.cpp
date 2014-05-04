#include "SASolver.h"
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <iostream>
#include <numeric>

SASolver::SASolver(
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrGField, 
    ClusterStructure const * const ptrCS,
    const double txTimeSlot
    ):
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), 
    m_ptrGField(ptrGField),
    m_maxIter(10000),
    m_txTimePerSlot(txTimeSlot),
    m_powerUpdater(ptrMap,ptrCS,txTimeSlot)
{
  Init();
  std::srand(std::time(NULL));
}

SASolver::~SASolver()
{

}

void
SASolver::Init()
{
  m_minPayoff = DBL_MAX;
  m_minVecSolution.resize(m_ptrMap->GetNumNodes());
  std::fill(m_minVecSolution.begin(), m_minVecSolution.end(), 0);
}


void
SASolver::InitSolution(double& objective, std::vector<int>& vecSolution)
{
  objective = DBL_MAX;
  vecSolution.resize(m_ptrMap->GetNumNodes());
  std::fill(vecSolution.begin(), vecSolution.end(), 0);
  std::list<std::list<int> >::const_iterator iterRow = m_ptrCS->GetListCluMemeber().begin();
  for (int i = 0 ;iterRow != m_ptrCS->GetListCluMemeber().end(); ++iterRow, ++i) {
    if (iterRow->size() > 1) {
      vecSolution.at(RandomSelectMember(i, *iterRow)) = 1;
    }
  }
}

int 
SASolver::RandomSelectMember( const int clusterIdx, const std::list<int>& listCluster)
{
  int supportName = -1;
  while(true) {
    int i = 0;
    std::list<int>::const_iterator iterCol = listCluster.begin();
    int idx = std::rand() % listCluster.size();
    while( i < idx){ 
      ++i;
      ++iterCol;
    }
    if (*iterCol == m_ptrCS->GetChNameByName(*iterCol)) {
      idx = std::rand() % listCluster.size();
    } else {
      supportName = *iterCol;
      break;
    }
  }
  return supportName;
}


int 
SASolver::RandomSelectCluster(std::list<std::list<int> >::const_iterator& iterCluster)
{
  iterCluster = m_ptrCS->GetListCluMemeber().begin();
  int chIdx = std::rand() % m_ptrCS->GetNumHeads();
  for (int i = 0; i < chIdx; ++i) ++iterCluster; 
  return chIdx;
}

double
SASolver::Solve(std::vector<int>& vecSolution)
{
  std::vector<int> vecTmpSolu;
  double objective = 0.0;
  m_minPayoff = DBL_MAX;
  InitSolution(objective, vecTmpSolu);
  for (int i = 0; i < m_maxIter; ++i) {
    Move(vecTmpSolu);
    if (IsFeasible(vecTmpSolu)) {
      objective = Optimize(vecTmpSolu);
    }
    else if(CoolProcess(objective)){
      objective = Optimize(vecTmpSolu);
    }
    UpdateBest(objective, vecTmpSolu);
  }
  vecSolution.assign(m_minVecSolution.begin(), m_minVecSolution.end());
}

void
SASolver::Move(std::vector<int>& vecSupport)
{
  std::list<std::list<int> >::const_iterator iterCluster;
  int chIdx = RandomSelectCluster(iterCluster);
  while(iterCluster->size() == 1) chIdx = RandomSelectCluster(iterCluster);
  int memberIdx = RandomSelectMember(chIdx, *iterCluster);
  /* de-activate the origin select node */
  for (int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
    if ((m_ptrCS->GetChNameByName(i) == m_ptrCS->GetChNameByName(memberIdx)) &&
        (m_ptrCS->GetChNameByName(i) != i) && 
        (vecSupport.at(i) == 1)) {
      vecSupport.at(i) = 0;
    }
  }
  /* activate the random select node */
  vecSupport.at(memberIdx) = 1;
}

template< class T>
string
SASolver::VecToString( const vector<T>& vec)
{
  stringstream ss;
  for (int i = 0; i < vec.size(); ++i) {
    ss << vec[i] << ' ';
  }
  return ss.str();
}

double
SASolver::Optimize(const std::vector<int>& vecSolution)
{
  std::vector<double> myVecPower(m_ptrMap->GetNumNodes());
  m_powerUpdater.Solve(myVecPower, vecSolution);
  return std::accumulate(myVecPower.begin(), myVecPower.end(), 0.0);
}

bool
SASolver::IsFeasible(const std::vector<int>& vecSolution)
{
  return true;
}

void
SASolver::UpdateBest(const double objective, const std::vector<int>& vecSolution)
{
  if (objective < m_minPayoff) {
    m_minPayoff = objective;
    m_minVecSolution.assign(vecSolution.begin(), vecSolution.end());
  }
}
