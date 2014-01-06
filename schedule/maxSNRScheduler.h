#ifndef _MAXSNRSCHEDULER_
#define _MAXSNRSCHEDULER_ 
#include "scheduler.h"
#include <vector>

class MaxSNRScheduler: public Scheduler
{
  public:
    MaxSNRScheduler(const double txTime, 
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE const * const, 
        ClusterStructure const * const);
    ~MaxSNRScheduler();
    double ScheduleOneSlot( std::vector<int>& );
    string PrintSelf(){ return m_type; }
  private:
    bool CheckFeasible( bool const * const supStru, double txTime2nd);
    bool CheckFeasible( const std::vector<int>& supStru, double txTime2nd);

    int                             m_numNodes;
    int                             m_numMaxHeads;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    double                          m_maxPower;
    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE const * const      m_ptrMatComputer;
    std::vector<double>                  m_vecNodePower;
    std::vector<int>                     m_vecSched;
    const string                    m_type;

};

#endif
