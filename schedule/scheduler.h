#ifndef _SCHEDULER_
#define _SCHEDULER_ 
#include <vector>
#include "map.h"
#include "CORRE_MA_OPE.h"
#include "clusterStructure.h"

class Scheduler
{
  public:
    Scheduler(){};
    Scheduler(const double txTime, Map const * const, CORRE_MA_OPE const * const, 
        ClusterStructure const * const);
    virtual ~Scheduler();
    virtual double ScheduleOneSlot(std::vector<int>& solution) = 0;
    virtual std::string PrintSelf() = 0;
  private:
};
#endif