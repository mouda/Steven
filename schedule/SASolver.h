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
    void InitSolution(double& objective, std::vector<int>& vecSolution);
    int  RandomSelectMember(const int clusterIdx, const std::list<int>& listCluster);
    int  RandomSelectCluster(std::list<std::list<int> >::const_iterator& iterCluster );
    void Move();
    double Optimize(const std::vector<int>& vecSolution);
    bool CoolProcess();
    bool IsFeasible(const std::vector<int>& vecSolution);
    bool CheckIfBest(const double obj);

    PowerUpdater m_powerUpdater;
    std::vector<int>                m_minVecSolution;
    double                          m_minPayoff;
    int                             m_maxIter;

    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE*                   m_ptrGField;
};
#endif
