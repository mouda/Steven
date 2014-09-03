#ifndef _SCHEDULERFACTORY_
#define _SCHEDULERFACTORY_ 
#include "scheduler.h"
#include "minPowerImageScheduler.h"

#include <string>
#include <iostream>
#include <list>

using std::string;
using std::vector;
using std::list;
using std::pair;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::stringstream;
using std::make_pair;

class SchedulerFactory
{
  public:
    SchedulerFactory(const double txTime, 
        const int tier2NumSlot,
        const double bandwidthKhz, 
        ImageMap const * const, 
        ImageSource* , 
        ClusterStructure const * const);

    ~SchedulerFactory();
    Scheduler* CreateScheduler(const string& );
  private:
    Scheduler*                      m_ptrSched;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    const double                    m_maxPower;
    const int                       m_tier2NumSlot;
    ImageMap const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    ImageSource*                    m_ptrImageSource;
};
#endif
