#ifndef _SIMULATOR_
#define _SIMULATOR_ 
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "map.h"
#include "clusterStructure.h"
#include "scheduler.h"
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

    void Run();
    void Run(const int numSlots);
    bool SelfCheck();

  private:
    void Print(const vector<int>& );
    string toString( const vector<int>& );
    Map*                m_ptrMap;
    ClusterStructure*   m_ptrCS;
    Scheduler*          m_ptrSched;
    CORRE_MA_OPE*       m_ptrGaussianField;
    

};
#endif
