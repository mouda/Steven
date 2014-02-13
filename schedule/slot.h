#ifndef _SLOT_
#define _SLOT_

#include <vector>
#include "gaussianField.h"
#include "map.h"
#include "scheduler.h"

class Slot
{
  public:
    Slot(Map* myMap, Scheduler* myScheduler);
    ~Slot();

    Slot* GetNextSlot() const;

  private:
    std::vector<int>    m_vecSupport;
    std::vector<double> m_vecEntropy;
    GaussianField*      m_ptrGaussianField;
    const Map*          m_ptrMap;
    Scheduler*          m_ptrSched; 

};

#endif
