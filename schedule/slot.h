#ifndef _SLOT_
#define _SLOT_

#include <vector>
#include "CORRE_MA_OPE.h"

class Slot
{
  public:
    Slot();
    ~Slot();

    Slot* GetNextSlot() const;

  private:
    std::vector<int>  m_vecSched;
    CORRE_MA_OPE*     m_fieldComputer;

};

#endif
