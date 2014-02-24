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
        const std::vector<double>&,
        const double
        );
    ~Slot();


    Slot* GetNextSlot() const;
    const std::vector<int>& GetSupport() const { return m_vecSupport; }
    const std::vector<double>& GetVariance() const { return m_vecVariance;}
    double GetEntropy() const { return m_entropy; }

  private:

    std::vector<int>    m_vecSupport;
    std::vector<double> m_vecVariance;
    double              m_entropy;
};

#endif
