#include "SASolver.h"
#include <cfloat>
#include <limits>

SASolver::SASolver(
    Map const * const ptrMap, 
    CORRE_MA_OPE* ptrGField, 
    ClusterStructure const * const ptrCS
    ):
    m_ptrMap(ptrMap), 
    m_ptrCS(ptrCS), 
    m_ptrGField(ptrGField)
{
  Init();
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
  m_minVecSolution.resize(m_ptrMap->GetNumNodes());
  std::fill(m_minVecSolution.begin(), m_minVecSolution.end(), 0);
}

void
SASolver::Move()
{

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
