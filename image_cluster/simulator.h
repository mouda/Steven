#ifndef _SIMULATOR_
#define _SIMULATOR_ 
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "imageMap.h"
#include "clusterStructure.h"
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
    Simulator(ImageMap* myMap,ClusterStructure* myCS, CORRE_MA_OPE* myField
        );
    ~Simulator();


    void WriteCS( const string& );
    void WriteWorseCaseTier2Power( const string&, const double, const double, const double);

    bool SelfCheck();

  private:
    bool CheckFeasible(const std::vector<int>& supStru, double txTime2nd);
    void Print(const std::vector<int>& );
    template < class T>
    string VecToString( const std::vector<T>& );

    ImageMap*                m_ptrMap;
    ClusterStructure*   m_ptrCS;
    CORRE_MA_OPE*       m_ptrGField;
    std::vector<int>*   m_vecSupport;
    std::vector<int>    m_vecTotal;

};
#endif
