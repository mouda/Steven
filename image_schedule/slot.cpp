#include "slot.h"

Slot::Slot(
    const std::vector<int>& vecSupport,
    const std::vector<double>& vecVariance,
    const std::vector<double>& vecPower,
    const double entropy,
    const double totalEntropy,
    const double MSE,
    const double totalPower
    ): 
  m_vecSupport(vecSupport),
  m_vecVariance(vecVariance),
  m_entropy(entropy),
  m_totalEntropy(totalEntropy),
  m_MSE(MSE),
  m_vecPower(vecPower),
  m_totalPower(totalPower)
{
}

Slot::~Slot()
{

}

Slot*
Slot::GetNextSlot() const  
{

}


