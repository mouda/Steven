#ifndef _NONSIMPLIFIEDSCHEDULER_
#define _NONSIMPLIFIEDSCHEDULER_

#include "scheduler.h"
class NonSimplifiedScheduler: public Scheduler
{
  public:
    NonSimplifiedScheduler(const double txTime, 
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE* , 
        ClusterStructure const * const);
    ~NonSimplifiedScheduler();

    void SetGaussianField(CORRE_MA_OPE* myGField) { m_ptrMatComputer = myGField;}
    double ScheduleOneSlot( std::vector<int>& );
    double ScheduleOneSlot( std::vector<int>& , const std::vector<double>& );
    string PrintSelf(){ return m_type; }
    double GetTxTimePerSlot() const { return m_txTimePerSlot; }

  private:
    int                             m_numNodes;
    int                             m_numMaxHeads;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    double                          m_maxPower;
    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE*                   m_ptrMatComputer;
    std::vector<double>             m_vecNodePower;
    std::vector<int>                m_vecSched; /* record the scheduled node */
    const string                    m_type;

}; 
#endif
