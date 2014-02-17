#ifndef _SIMULATOR_
#define _SIMULATOR_ 
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "map.h"
#include "slot.h"
#include "event.h"
#include "clusterStructure.h"
#include "scheduler.h"
#include "maxSNRScheduler.h"
#include "CORRE_MA_OPE.h"

using std::string;
using std::pair;
using std::cout;
using std::endl;
using std::cerr;
using std::fstream;
using std::ios;
using std::setfill;
using std::setw;
using std::stringstream;
using std::make_pair;

class Simulator
{
  public:
    Simulator();
    Simulator(Map* myMap,ClusterStructure* myCS, Scheduler* myScheduler );
    Simulator(Map* myMap,ClusterStructure* myCS, Scheduler* myScheduler, CORRE_MA_OPE* myField);
    ~Simulator();

    void SetEvents(double t_ms);
    void SequentalRun(double t_ms);
    void SequentalFree();
    void Run();
    void Run(const int numSlots);


    bool SelfCheck();
    std::vector<int>  CheckConnection(const std::vector<int>& );

  private:
    bool CheckFeasible(const std::vector<int>& supStru, double txTime2nd);
    double GetTotalEntropy(const std::vector<int>& vecSupport) const;
    double Get1stSlotEntropy(const std::vector<int>& vecSupport) const;
    void Print(const std::vector<int>& );
    string toString( const std::vector<int>& );
    Slot* GetNextSlot(Slot*);

    Map*                m_ptrMap;
    ClusterStructure*   m_ptrCS;
    Scheduler*          m_ptrSched;
    CORRE_MA_OPE*       m_ptrGaussianField;
    std::vector<int>*   m_vecSupport;

    /* event driven simulation data structure */
    std::list<Slot*>    m_listSlot;
    Slot*               m_ptrSlotHead;
    std::list<Event*>   m_listEvent;
};
#endif
