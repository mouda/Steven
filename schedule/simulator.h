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
#include "fileHandler.h"

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
    Simulator(Map* myMap,ClusterStructure* myCS, Scheduler* myScheduler, CORRE_MA_OPE* myField, 
        const string&,
        const string&,
        const string&,
        const string&, 
        const string&);
    ~Simulator();

    void SetEvents(double t_ms);
    void SequentialRun(int tier2NumSlot);
    void SequentialFree();
    void Run();

    void WriteCS( const string& );
    void WriteEntropy();
    void WriteMSE();
    void WriteSolution();
    void WriteSupport();
    void WriteTotalEntropy();


    bool SelfCheck();
    std::vector<int>  CheckConnection(const std::vector<int>& );

  private:
    bool CheckFeasible(const std::vector<int>& supStru, double txTime2nd);
    double GetTotalEntropy(const std::vector<int>& vecSupport) const;
    double Get1stSlotEntropy(const std::vector<int>& vecSupport) const;
    void Print(const std::vector<int>& );
    template < class T>
    string VecToString( const std::vector<T>& );
    Slot* GetNextSlot(Slot*, std::vector<int>& );

    Map*                m_ptrMap;
    ClusterStructure*   m_ptrCS;
    Scheduler*          m_ptrSched;
    CORRE_MA_OPE*       m_ptrGField;
    std::vector<int>*   m_vecSupport;

    std::vector<int>    m_vecTotal;

    /* event driven simulation data structure */
    std::list<Slot*>    m_listSlot;
    Slot*               m_ptrSlotHead;
    std::list<Event*>   m_listEvent;

    /* output  */
    FileHandler         m_entropyFHandler;
    FileHandler         m_MSEFHandler;
    FileHandler         m_solutionFHandler;
    FileHandler         m_supportNumFHandler;
    FileHandler         m_totalEntropyFHandler;
};
#endif
