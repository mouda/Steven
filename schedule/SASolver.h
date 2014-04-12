#ifndef _SASOLVER_
#define _SASOLVER_

#include <vector>
#include "powerUpdater.h"

#include "map.h"
#include "CORRE_MA_OPE.h"
#include "clusterStructure.h"

class SASolver
{
  public:
    SASolver(
        Map const * const, 
        CORRE_MA_OPE* , 
        ClusterStructure const * const
        );
    ~SASolver();
    double Solve(std::vector<int>& );

  private:
    void Init();
    void Move();
    void Optimize();
    void CheckIfFeasible();
    void CheckIfBest();

    PowerUpdater m_powerUpdater;
    std::vector<int>                m_vecSolution;
    std::vector<int>                m_minVecSolution;
    double                          m_payoff;
    double                          m_minPayoff;

    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE*                   m_ptrGField;
};
#endif
