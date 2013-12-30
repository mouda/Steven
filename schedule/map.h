#ifndef _MAP_
#define _MAP_
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "simSystem.h"
#include "CORRE_MA_OPE.h"

using std::vector;
using std::pair;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::fstream;
using std::ios;
using std::stringstream;
using std::make_pair;

class Map 
{
  public:
    Map(const int numNodes, const int numHeads, const double maxPower, const double corrFactor, const int mapId);
    ~Map();
    void SetChannelByXYPair(const vector<pair<double,double> >& );
    const vector<vector<double> >& GetGij();
    double** const GetMatDistance(){ return m_matDistance;}

  private:
    const int               m_mapId;
    const int               m_numNodes;
    const int               m_numInitHeads;
    const double            m_maxPower;
    const double            m_corrFactor;
    vector<vector<double> > m_matGij;
    double**                m_matDistance;
    vector<double>          m_vecPower;
    SimSystem*              m_systemComputer;

    

};
#endif
