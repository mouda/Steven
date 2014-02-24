#include "slot.h"

Slot::Slot(
    const std::vector<int>& vecSupport,
    const std::vector<double>& vecVariance,
    const double entropy
    ): 
  m_vecSupport(vecSupport),
  m_vecVariance(vecVariance),
  m_entropy(entropy)
{
}

Slot::~Slot()
{

}

Slot*
Slot::GetNextSlot() const  
{

}


