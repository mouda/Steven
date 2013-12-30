#ifndef _SCHEDULERFACTORY_
#define _SCHEDULERFACTORY_ 
#include "scheduler.h"

class SchedulerFactory
{
  public:
    Scheduler* CreateScheduler();
  private:
    Scheduler* m_ptrSched;
};
#endif
