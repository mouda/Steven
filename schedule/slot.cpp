#include "slot.h"

Slot::Slot(
    const std::vector<int>& vecSupport,
    const std::vector<double>& vecVariance,
    const double entropy,
    const double MSE
    ): 
  m_vecSupport(vecSupport),
  m_vecVariance(vecVariance),
  m_entropy(entropy),
  m_MSE(MSE)
{
}

Slot::~Slot()
{

}

Slot*
Slot::GetNextSlot() const  
{

}


