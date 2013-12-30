#ifndef _SIMULATOR_
#define _SIMULATOR_ 
#include "map.h"
#include "clusterStructure.h"
#include "scheduler.h"
class Simulator
{
  public:
    Simulator();
    Simulator(const Map& myMap, const ClusterStructure& myCS, const Scheduler& myScheduler );
    ~Simulator();

    void Run();

  private:

};
#endif
