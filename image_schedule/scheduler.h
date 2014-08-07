#ifndef _SCHEDULER_
#define _SCHEDULER_ 
#include <vector>
#include "imageMap.h"
#include "imageSource.h"
#include "clusterStructure.h"

class Scheduler
{
  public:
    Scheduler(){};
    Scheduler(const double txTime, ImageMap const * const, ImageSource const * const, 
        ClusterStructure const * const);
    virtual void SetGaussianField(ImageSource* ) = 0;
    virtual ~Scheduler();
    virtual double ScheduleOneSlot(std::vector<int>& solution) = 0;
    virtual double ScheduleOneSlot(std::vector<int>& solution, std::vector<double>& vecPower, const std::vector<double>& vecVariance) = 0;
    virtual std::string PrintSelf() = 0;
    virtual double GetTxTimePerSlot() const = 0;
};
#endif
