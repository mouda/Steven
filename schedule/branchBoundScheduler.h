#ifndef _BBSECHEDULER_
#define _BBSECHEDULER_ 
#include "scheduler.h"

class BranchBoundScheduler: public Scheduler
{
  public:
    BranchBoundScheduler(const double txTime, 
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE const * const, 
        ClusterStructure const * const);
    ~BranchBoundScheduler();
    double ScheduleOneSlot( vector<int>& );
    string PrintSelf(){ return m_type; }
  private:
    int                             m_numNodes;
    int                             m_numMaxHeads;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    double                          m_maxPower;
    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE const * const      m_ptrMatComputer;
    vector<double>                  m_vecNodePower;
    const string                    m_type;

};
#endif
