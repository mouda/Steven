#ifndef _SLOT_
#define _SLOT_

#include <vector>
#include "scheduler.h"

class Slot
{
  public:
    Slot(
        const std::vector<int>&, 
        const std::vector<double>&,
        const std::vector<double>&,
        const double,
        const double,
        const double,
        const double
        );
    ~Slot();


    Slot* GetNextSlot() const;
    const std::vector<int>& GetSupport() const { return m_vecSupport; }
    const std::vector<double>& GetVariance() const { return m_vecVariance;}
    const std::vector<double>& GetVecPower() const { return m_vecPower; }
    double GetEntropy() const { return m_entropy; }
    double GetTotalEntropy() const { return m_totalEntropy;}
    double GetMSE() const { return m_MSE; } 

  private:

    std::vector<int>    m_vecSupport;
    std::vector<double> m_vecVariance;
    std::vector<double> m_vecPower;
    double              m_totalPower; /* total power per slot */
    double              m_entropy;
    double              m_totalEntropy;
    double              m_MSE;
};

#endif
