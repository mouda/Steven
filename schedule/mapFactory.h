#ifndef _MAPFACTORY_
#define _MAPFACTORY_ 
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "map.h"
#include "CORRE_MA_OPE.h"

using std::string;
using std::pair;
using std::cout;
using std::endl;
using std::cerr;
using std::fstream;
using std::ios;
using std::stringstream;
using std::make_pair;

class MapFactory 
{
  public:
    MapFactory(const string& mapFileName, const double maxPower, 
        const double corrFactor, const double quantizationBits,
        const double bandwidthKhz,  
        const int maxNumHead, const int numNodes);
    ~MapFactory();
    Map* CreateMap();
    CORRE_MA_OPE* CreateMatrixComputer();

  private:
    string                        m_mapFileName;
    const double                  m_maxPower;
    const int                     m_maxNumHead;
    const int                     m_numNodes;
    const double                  m_corrFactor;
    const double                  m_quantizationBits;
    const double                  m_bandwidthKhz;
    int                           m_mapId;
    Map*                          m_ptrMap;
    vector<pair<double, double> > m_vecPairPos;
    CORRE_MA_OPE*                 m_ptrMatComputer;

};
#endif