#include "mapFactory.h"
#include "utility.h"
#include <cstdio>

MapFactory::MapFactory(const string& mapFileName, const double maxPower, 
    const double corrFactor, const double quantizationBits,
    const double bandwidthKhz,  
    const int maxNumHead, const int numNodes):
  m_ptrMap(0), m_ptrMatComputer(0), m_maxPower(maxPower), 
  m_corrFactor(corrFactor), m_quantizationBits(quantizationBits),
  m_bandwidthKhz(bandwidthKhz), 
  m_maxNumHead(maxNumHead), m_numNodes(numNodes)
{
  
#ifdef DEBUG
  cout << mapFileName << endl;
#endif
  m_mapFileName = mapFileName;
  vector<string> macroTokens = split(mapFileName,'_'); 
  vector<string> microTokens = split(macroTokens[6], '.');
  m_mapId = atoi(microTokens[0].c_str());

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
MapFactory::CreateMap()
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
  m_ptrMap = new Map(numNodes, m_maxNumHead, m_maxPower, m_corrFactor, m_quantizationBits, m_bandwidthKhz, m_mapId);
  m_ptrMap->SetChannelByXYPair(m_vecPairPos);
  return m_ptrMap;
}

CORRE_MA_OPE*
MapFactory::CreateMatrixComputer()
{
  m_ptrMatComputer = new CORRE_MA_OPE(m_numNodes, m_corrFactor, m_ptrMap->GetMatDistance());
  cout << m_numNodes << endl;
  m_ptrMatComputer->returnNSetCorrelationFactorByCompressionRatio(m_corrFactor , m_ptrMap->GetIdtEntropy() ,static_cast<double>(m_numNodes));
  return m_ptrMatComputer;
}
