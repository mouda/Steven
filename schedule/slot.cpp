#include "slot.h"

Slot::Slot(Map* myMap, Scheduler* myScheduler): 
  m_ptrMap(myMap), m_ptrSched(myScheduler)
{
  m_vecSupport.resize(m_ptrMap->GetNumNodes());
  std::fill(m_vecSupport.begin(), m_vecSupport.end(), 0);
  m_ptrSched->ScheduleOneSlot(m_vecSupport);
}

Slot::~Slot()
{

}

Slot*
Slot::GetNextSlot() const  
{

}
