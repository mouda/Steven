#ifndef _MAPFACTORY_
#define _MAPFACTORY_ 
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "imageMap.h"
#include "imageSource.h"

using std::string;
using std::pair;
using std::cout;
using std::endl;
using std::cerr;
using std::fstream;
using std::ios;
using std::stringstream;
using std::make_pair;

class ImageMapFactory 
{
  public:
    ImageMapFactory(
        const string& mapFileName,
        const double maxPower, 
        const double spatialCorrFactor, 
        const double bandwidthKhz,  
        const int maxNumHead, 
        const int numNodes
        );
    ~ImageMapFactory();
    ImageMap* CreateMap( bool myImageFlag );
    ImageSource* CreateImageSource();

  private:
    string                        m_mapFileName;
    const double                  m_maxPower;
    const int                     m_maxNumHead;
    const int                     m_numNodes;
    const double                  m_spatialCorrFactor;
    const double                  m_bandwidthKhz;
    int                           m_mapId;
    ImageMap*                          m_ptrMap;
    vector<pair<double, double> > m_vecPairPos;
    ImageSource*                  m_ptrImageSource;

};
#endif
