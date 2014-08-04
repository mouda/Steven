#include "imageMapFactory.h"
#include "utility.h"
#include <cstdio>

ImageMapFactory::ImageMapFactory(
    const string& mapFileName, 
    const string& corrFileName,
    const string& idtFileName,
    const double maxPower, 
    const double bandwidthKhz,  
    const int maxNumHead, const int numNodes):
  m_ptrMap(0),  
  m_ptrImageSource(0), 
  m_maxPower(maxPower), 
  m_bandwidthKhz(bandwidthKhz), 
  m_maxNumHead(maxNumHead), 
  m_numNodes(numNodes),
  m_strIdtFName(idtFileName),
  m_strCorrFName(corrFileName)
{
  
#ifdef DEBUG
  cout << mapFileName << endl;
#endif
  m_mapFileName = mapFileName;
//  vector<string> macroTokens = split(mapFileName,'_'); 
//  vector<string> microTokens = split(macroTokens[6], '.');
  
  m_mapId = 0;
}

ImageMapFactory::~ImageMapFactory()
{
  if (m_ptrMap != 0) {
    delete m_ptrMap;
  }
  if (m_ptrImageSource != 0) {
    delete m_ptrImageSource;
  }
}

ImageMap*
ImageMapFactory::CreateMap(bool myImageFlag)
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
  m_ptrMap = new ImageMap(numNodes, m_maxNumHead, m_maxPower, m_bandwidthKhz, m_mapId);
  m_ptrMap->SetChannelByXYPair(m_vecPairPos);
  if (myImageFlag) {
    fstream myImageIdtFile(m_strIdtFName.c_str(),std::ios::in);
    vector<double> vecIdtCodingBits;
    for (int i = 0; i < m_numNodes; ++i) {
      double tmp = 0;
      myImageIdtFile >> tmp;
      vecIdtCodingBits.push_back(tmp*8);
    }
    m_ptrMap->SetIdtImageCodingBits(vecIdtCodingBits);
    m_ptrMap->SetIdtImageFlag();

    myImageIdtFile.close();
  }
  return m_ptrMap;
}


ImageSource*
ImageMapFactory::CreateImageSource()
{
  m_ptrImageSource = new ImageSource(m_numNodes, m_strIdtFName, m_strCorrFName, m_ptrMap->GetMatDistance());
  return m_ptrImageSource;
}
