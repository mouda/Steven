#include "SASolver.h"
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <iostream>

SASolver::SASolver(
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrGField, 
    ClusterStructure const * const ptrCS
    ):
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), 
    m_ptrGField(ptrGField),
    m_maxIter(1000)
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
  m_payoff = -DBL_MAX;
  m_minPayoff = -DBL_MAX;
  m_vecSolution.resize(m_ptrMap->GetNumNodes());
  std::fill(m_vecSolution.begin(), m_vecSolution.end(), 0);

  std::list<std::list<int> >::const_iterator iterRow = m_ptrCS->GetListCluMemeber().begin();
  for (int i = 0 ;iterRow != m_ptrCS->GetListCluMemeber().end(); ++iterRow, ++i) {
    if (iterRow->size() > 1) {
      m_vecSolution.at(i) = RandomSelectMember(i, *iterRow);
    }
  }
  m_minVecSolution.resize(m_ptrMap->GetNumNodes());
  std::fill(m_minVecSolution.begin(), m_minVecSolution.end(), 0);
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
  for (int i = 0; i < m_maxIter; ++i) {
    Move();
    Optimize();
    CheckIfFeasible();
    CheckIfBest();
  }
}

void
SASolver::Move()
{
  std::list<std::list<int> >::const_iterator iterCluster;
  int chIdx = RandomSelectCluster(iterCluster);
  while( iterCluster->size() == 1 ) chIdx = RandomSelectCluster(iterCluster);
  int memberIdx = RandomSelectMember(chIdx, *iterCluster);

}

void
SASolver::Optimize()
{

}

void
SASolver::CheckIfFeasible()
{

}

void
SASolver::CheckIfBest()
{

}
