#include "schedFactory.h"

SchedulerFactory::SchedulerFactory(const double txTime, 
  const double bandwidthKhz, 
  Map const * const ptrMap, 
  CORRE_MA_OPE* ptrMatComputer, 
  ClusterStructure const * const ptrCS): 
  m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
  m_ptrMap(ptrMap), m_maxPower(0),
  m_ptrCS(ptrCS), m_ptrMatComputer(ptrMatComputer),m_ptrSched(0)
{

}

SchedulerFactory::~SchedulerFactory()
{
  if (m_ptrSched != 0) {
    delete m_ptrSched;
  }
}


Scheduler*
SchedulerFactory::CreateScheduler( const string& scheduleType)
{
  if (m_ptrSched != 0) {
    return m_ptrSched;
  }
  if (scheduleType == "Baseline") {
    m_ptrSched = new MaxSNRScheduler(m_txTimePerSlot, m_bandwidthKhz, 
        m_ptrMap, m_ptrMatComputer, m_ptrCS); 
    return m_ptrSched; 
  }
  if (scheduleType == "Branchbound") {
    m_ptrSched = new BranchBoundScheduler(m_txTimePerSlot, m_bandwidthKhz, 
        m_ptrMap, m_ptrMatComputer, m_ptrCS); 
    return m_ptrSched; 
  }
  else{
    return NULL;
  }

}
