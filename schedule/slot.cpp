#include "slot.h"

Slot::Slot(
    const std::vector<int>& vecSupport,
    const std::vector<double>& vecVariance
    ): 
  m_vecSupport(vecSupport),
  m_vecVariance(vecVariance)
{
}

Slot::~Slot()
{

}

Slot*
Slot::GetNextSlot() const  
{

}
