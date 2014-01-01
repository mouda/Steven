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
    ~Simulator();

    void Run();
    bool SelfCheck();

  private:
    Map*                m_ptrMap;
    ClusterStructure*   m_ptrCS;
    Scheduler*          m_ptrSched;
    

};
#endif
