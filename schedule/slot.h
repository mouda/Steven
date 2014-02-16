#ifndef _SLOT_
#define _SLOT_

#include <vector>
#include "map.h"
#include "scheduler.h"

class Slot
{
  public:
    Slot(
        const std::vector<int>&, 
        const std::vector<double>& 
        );
    ~Slot();

    Slot* GetNextSlot() const;

  private:
    std::vector<int>    m_vecSupport;
    std::vector<double> m_vecVariance;
};

#endif
