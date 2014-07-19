#include "mapFactory.h"
#include "utility.h"
#include <cstdio>

MapFactory::MapFactory(
    const string& mapFileName, 
    const double maxPower, 
    const double spatialCorrFactor, 
    const double temporalCorrFactor,
    const double quantizationBits,
    const double bandwidthKhz,  
    const int maxNumHead, const int numNodes):
  m_ptrMap(0), m_ptrMatComputer(0), m_maxPower(maxPower), 
  m_spatialCorrFactor(spatialCorrFactor), 
  m_temporalCorrFactor(temporalCorrFactor),
  m_quantizationBits(quantizationBits),
  m_bandwidthKhz(bandwidthKhz), 
  m_maxNumHead(maxNumHead), m_numNodes(numNodes)
{
  
#ifdef DEBUG
  cout << mapFileName << endl;
#endif
  m_mapFileName = mapFileName;
//  vector<string> macroTokens = split(mapFileName,'_'); 
//  vector<string> microTokens = split(macroTokens[6], '.');
  m_mapId = 0;

}

MapFactory::~MapFactory()
{
  if (m_ptrMap != 0) {
    delete m_ptrMap;
  }
  if (m_ptrMatComputer != 0) {
    delete m_ptrMatComputer;
  }
}

Map*
MapFactory::CreateMap(bool myImageFlag)
{
  fstream mapFile;
  mapFile.open(m_mapFileName.c_str(), fstream::in);
  string strLine;
  getline(mapFile, strLine);
  vector<string> tokens = split(strLine,' ');
  int numNodes = atoi(tokens[0].c_str());
  int varRarius = atof(tokens[1].c_str());
#ifdef DEBUG 
  cout << "Nodes: " << m_numNodes << " Power: " << m_maxPower << endl;
  cout << "CHs: " << m_maxNumHead  << " mapId: " << m_mapId << endl;
#endif
  m_vecPairPos.clear();
  while(getline(mapFile, strLine)){
    vector<string> posTokens = split(strLine, ' ');
    stringstream ss(strLine);
    double x, y;
    ss >> x >> y;
    m_vecPairPos.push_back(make_pair(x,y));
  }
  if ( m_vecPairPos.size() != numNodes || m_numNodes != numNodes) {
    cerr << "Error: not match declarion and # of position pairs " << endl;
    return NULL;
  }
  mapFile.close();
  m_ptrMap = new Map(numNodes, m_maxNumHead, m_maxPower, m_spatialCorrFactor, m_quantizationBits, m_bandwidthKhz, m_mapId);
  m_ptrMap->SetChannelByXYPair(m_vecPairPos);
  if (myImageFlag) {
    fstream myImageIdtFile("paper720_30cam_indepByte.txt",std::ios::in);
    vector<double> vecIdtCodingBits;
    for (int i = 0; i < m_numNodes; ++i) {
      double tmp = 0;
      myImageIdtFile >> tmp;
      vecIdtCodingBits.push_back(tmp*8);
    }
    m_ptrMap->SetIdtImageCodingBits(vecIdtCodingBits);

    myImageIdtFile.close();
  }
  return m_ptrMap;
}

CORRE_MA_OPE*
MapFactory::CreateMatrixComputer()
{
  m_ptrMatComputer = new CORRE_MA_OPE(m_numNodes, m_spatialCorrFactor, m_temporalCorrFactor, m_ptrMap->GetMatDistance(), m_quantizationBits);
  return m_ptrMatComputer;
}

ImageSource*
MapFactory::CreateImageSource()
{
  m_ptrImageSource = new ImageSource(m_numNodes, m_spatialCorrFactor, m_temporalCorrFactor, m_ptrMap->GetMatDistance(), m_quantizationBits);
  return m_ptrImageSource;

}
