#ifndef _MAXSNRSCHEDULER_
#define _MAXSNRSCHEDULER_ 
#include "scheduler.h"
#include <vector>

class MaxSNRScheduler: public Scheduler
{
  public:
    MaxSNRScheduler(
        const double txTime, 
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE* , 
        ClusterStructure const * const);
    ~MaxSNRScheduler();

    void SetGaussianField(CORRE_MA_OPE* myGField) { m_ptrMatComputer = myGField;}
    double ScheduleOneSlot( std::vector<int>& );
    double ScheduleOneSlot( std::vector<int>& , std::vector<double>&, const std::vector<double>& );
    string PrintSelf(){ return m_type; }
    double GetTxTimePerSlot() const { return m_txTimePerSlot; }
  private:
    bool CheckFeasible( bool const * const supStru, double txTime2nd);
    bool CheckFeasible( const std::vector<int>& supStru, double txTime2nd);
    /* to be refactored  */
    bool CheckGroupScheduled( const int );
    bool CheckAllScheduled();
    bool ResetGroupScheduled( const int);

    int                             m_numNodes;
    int                             m_numMaxHeads;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    double                          m_maxPower;
    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE*                   m_ptrMatComputer;
    std::vector<double>             m_vecNodePower;
    std::vector<int>                m_vecSched;
    const string                    m_type;

};

#endif
