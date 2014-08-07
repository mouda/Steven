#include "schedFactory.h"

SchedulerFactory::SchedulerFactory(const double txTime, 
  const int tier2NumSlot,
  const double bandwidthKhz, 
  ImageMap const * const ptrMap, 
  ImageSource* ptrImageSource, 
  ClusterStructure const * const ptrCS): 
  m_txTimePerSlot(txTime), m_bandwidthKhz(bandwidthKhz),
  m_tier2NumSlot(tier2NumSlot),
  m_ptrMap(ptrMap), m_maxPower(0),
  m_ptrCS(ptrCS), m_ptrImageSource(ptrImageSource),m_ptrSched(0)
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
  else if (scheduleType == "ImageSource") {
    m_ptrSched = new MinPowerImageScheduler(m_txTimePerSlot, m_tier2NumSlot, m_bandwidthKhz,
        m_ptrMap, m_ptrImageSource, m_ptrCS);
    return m_ptrSched;
  }
  else{
    return NULL;
  }

}
