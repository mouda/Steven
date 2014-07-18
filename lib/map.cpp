#include "map.h"

Map::Map(
    const int numNodes, 
    const int numHeads, 
    const double maxPower, 
    const double corrFactor, 
    const double quantizationBits, 
    const double bandwidthKhz, 
    const int mapId):
  m_numNodes(numNodes), m_numInitHeads(numHeads), m_maxPower(maxPower), m_mapId(mapId),
  m_corrFactor(corrFactor), m_quantizationBits(quantizationBits), 
  m_bandwidthKhz(bandwidthKhz), m_systemComputer(0)
{
  m_systemComputer = new SimSystem;
  m_idtEntropy = 0.5*log2(2*PI*exp(1)) + quantizationBits;
  m_realNoise = m_systemComputer->returnInBandThermalNoise(bandwidthKhz);
}

Map::~Map()
{
  if (m_systemComputer != 0) {
    delete m_systemComputer;
  }
  if (m_matDistance != 0) {
    for(int i = 0; i < m_numNodes; ++i){
      if (m_matDistance[i] != 0) {
        delete [] m_matDistance[i];
      }
    }
    delete [] m_matDistance;
  }

}

void
Map::SetChannelByXYPair(const vector<pair<double, double> >& vecPosPair)
{
#ifdef DEBUG
  for (int i = 0; i < vecPosPair.size(); i++) {
    cout << vecPosPair[i].first << ' ' << vecPosPair[i].second << endl;
  }
#endif
  m_matDistance = new double* [m_numNodes]; 
  m_matGij.resize(m_numNodes, vector<double>(m_numNodes,0.0));
  m_vecGi0.resize(m_numNodes);
  m_vecPairNodePos = vecPosPair;
  for (int i = 0; i < m_numNodes; i++) {
    double lhsX = vecPosPair[i].first;
    double lhsY = vecPosPair[i].second;
    m_matDistance[i] = new double [m_numNodes];
    m_vecGi0.at(i) = m_systemComputer->returnChannelGain_BS_ByPos(lhsX, lhsY ); 
    for (int j = 0; j < m_numNodes; j++) {
      double rhsX = vecPosPair[j].first;
      double rhsY = vecPosPair[j].second;
      if ( i == j ) {
        m_matGij[i][j] = 1.0;
        m_matDistance[i][j] = 0.0;
      }
      else if ( i > j ) {
        m_matGij[i][j] = m_matGij[j][i];
        m_matDistance[i][j] = m_matDistance[j][i];
      }
      else {
        m_matGij[i][j] = 
          m_systemComputer->returnChannelGainByPos(lhsX, lhsY, rhsX, rhsY);
        m_matDistance[i][j] = 
          pow(lhsX - rhsX, 2) + pow(lhsY - rhsY,2);
      }
    }
  }
  for (int i = 0; i < m_numNodes; i++) {
    for (int j = 0; j < m_numNodes; j++) {
      //cout << "i: " << i << "j: " << j << " " << m_matGij[i][j] << ' ';
//      cout << "i: " << i << vecPosPair[i].first <<' ' <<vecPosPair[i].second<< " j: " << j << " "<< vecPosPair[j].first <<' ' <<vecPosPair[j].second<<' ' <<m_matDistance[i][j] << endl;
    }
  }
}

const vector<vector<double> >&
Map::GetGij()
{
  return m_matGij;
}


void 
Map::SetIdtImageCodingBits(const vector<double>& myIdtEntropy)
{
  m_vecIdtEntropy.assign(myIdtEntropy.begin(),myIdtEntropy.end());
}
