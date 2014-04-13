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
    int  RandomSelectMember(const int clusterIdx, const std::list<int>& listCluster);
    int RandomSelectCluster(std::list<std::list<int> >::const_iterator& iterCluster );
    void Move();
    void Optimize();
    void CheckIfFeasible();
    void CheckIfBest();

    PowerUpdater m_powerUpdater;
    std::vector<int>                m_vecSolution;
    std::vector<int>                m_minVecSolution;
    double                          m_payoff;
    double                          m_minPayoff;
    int                             m_maxIter;

    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE*                   m_ptrGField;
};
#endif
