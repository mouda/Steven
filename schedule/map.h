#ifndef _MAP_
#define _MAP_
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "simSystem.h"
#include "CORRE_MA_OPE.h"

#define PI 3.1415
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
    Map(const int numNodes, const int numHeads, const double maxPower, 
        const double corrFactor, const double quantizationBits,
        const double bandwidthKhz, const int mapId);
    ~Map();
    void SetChannelByXYPair(const vector<pair<double,double> >& );
    const vector<vector<double> >& GetGij();
    double** const GetMatDistance(){ return m_matDistance;}
    int GetMapId() const { return m_mapId; }
    int GetNumNodes() const { return m_numNodes; }
    int GetNumInitHeads() const { return m_numInitHeads; }
    double GetQBits() const { return m_quantizationBits; }
    double GetNodeXPos(const int idx) const { return m_vecPairNodePos[idx].first;}
    double GetNodeYPos(const int idx) const { return m_vecPairNodePos[idx].second;}
    double GetGijByPair( const int lhs, const int rhs) const {return m_matGij[lhs][rhs]; } 
    double GetGi0ByNode( const int idx ) const {return m_vecGi0.at(idx); } 
    double GetMaxPower() const {return m_maxPower; }
    double GetIdtEntropy() const { return m_idtEntropy; }
    double GetNoise() const { return m_realNoise; } 
    double GetBandwidth() const { return m_bandwidthKhz; }

  private:
    const int                     m_mapId;
    const int                     m_numNodes;
    const int                     m_numInitHeads;
    const double                  m_maxPower;
    const double                  m_corrFactor;
    const double                  m_bandwidthKhz;
    double                        m_quantizationBits;
    double                        m_idtEntropy;
    double                        m_realNoise;
    vector<vector<double> >       m_matGij;
    vector<double>                m_vecGi0;
    vector<pair<double, double> > m_vecPairNodePos;
    double**                      m_matDistance;
    vector<double>                m_vecPower;
    SimSystem*                    m_systemComputer;

};
#endif
