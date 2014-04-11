#ifndef _SASOLVER_
#define _SASOLVER_

#include <vector>
#include "powerUpdater.h"

class SASolver
{
  public:
    SASolver();
    ~SASolver();
    double Solve(std::vector<int>& );
  private:
    void Init();
    void Move();
    void Optimize();
    void CheckIfFeasible();
    void CheckIfBest();

    PowerUpdater m_powerUpdater;
};
#endif
