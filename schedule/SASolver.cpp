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
