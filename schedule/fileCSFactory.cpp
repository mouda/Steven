#include "fileCSFactory.h"
#include <fstream>

FileCSFactory::FileCSFactory( Map const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  CsFactory(myMap, myMatComputer)
{

}

FileCSFactory::~FileCSFactory()
{

}

ClusterStructure*
FileCSFactory::CreateClusterStructure()
{
  if (m_ptrCS == 0) {
    std::vector<int> myVecHeadNames(m_ptrMap->GetNumInitHeads());
    std::list<std::list<int> > myListCluMembers;
    if (ReadCSFromFile(myVecHeadNames, myListCluMembers)) {
      m_ptrCS = new ClusterStructure(m_ptrMap->GetNumNodes(), 
          m_ptrMap->GetNumInitHeads() );
      m_ptrCS->SetRecord(myVecHeadNames, myListCluMembers, vector<int>(m_ptrMap->GetNumNodes(),1));
      return m_ptrCS;
    }
    else{
      return NULL;
    }
  }
}

bool
FileCSFactory::ReadCSFromFile(std::vector<int>& vecHeadNames, std::list<list<int> >& lisCluMembers)
{
  
  fstream CSFile;
  CSFile.open(m_CSFName.c_str(), std::fstream::in);
  int numNodes = 0; 
  int numHeads = 0;
  CSFile >> numNodes; 
  CSFile >> numHeads;
  if ( numNodes != m_ptrMap->GetNumNodes()) {
    cerr << "Error: mismatch number of nodes: " << numNodes << ' ' << m_ptrMap->GetNumNodes() << endl;
    return false;
  }
  vecHeadNames.resize(m_ptrMap->GetNumInitHeads());
  if ( numHeads != m_ptrMap->GetNumInitHeads() ) {
    cerr << "Error: mismatch number of heads: " << numHeads << ' ' << m_ptrMap->GetNumInitHeads() << endl;
    return false;
  }
  string dummy;
  for (int i = 0; i < 8; ++i) {
    CSFile >> dummy; 
  }
  int indicator = 0;
  for (int i = 0; i < numHeads; ++i) {
    for (int j = 0; j < numNodes; ++j) {
      CSFile >> indicator;
    }
  }
  int headName = -1;
  for (int i = 0; i < numHeads; ++i) {
    CSFile >> headName; 
    vecHeadNames.at(i)= headName - 1;
    cout << headName  - 1<< endl;
  }
  for (int i = 0; i < vecHeadNames.size(); ++i) {
    cout << vecHeadNames.at(i) << ' ';
  }
  cout << endl;
  double nodePower = 0.0;
  for (int i = 0; i < numNodes; ++i) {
    CSFile >> nodePower;
  }
  int nodeName;
  lisCluMembers.resize(m_ptrMap->GetNumInitHeads());
  list<list<int> >::iterator iterRow = lisCluMembers.begin(); 
  string line;
  std::getline(CSFile,line,'\n');
  for (int i = 0; i < numHeads; ++i, ++iterRow ) {
    std::getline(CSFile, line, '\n');
    stringstream ss(line);
    int ch;
    while(true){
      ss >> nodeName;
      if (ss.eof()) break; 
      iterRow->push_back(nodeName -1 );
    }
    line.clear();
  }
  return true;
}
