#include "csFactory.h"
#include <limits>
CsFactory::CsFactory( Map const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  m_ptrMap(myMap), m_ptrMatComputer(myMatComputer),
  m_ptrCS(0)
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads = m_ptrMap->GetNumInitHeads();
}

CsFactory::~CsFactory()
{
  if (m_ptrCS != 0) {
    delete m_ptrCS;
  }
}


ClusterStructure*
CsFactory::CreateClusterStructure()
{
  if (m_ptrCS == 0) {
    m_ptrCS = new ClusterStructure(m_ptrMap->GetNumNodes(), 
        m_ptrMap->GetNumInitHeads() );
  }
}

void
CsFactory::Kmedoid()
{
  int retryTimes = 0;
  double tempHeadX [m_numMaxHeads];
  double tempHeadY [m_numMaxHeads];
  int tempHeadList [m_numMaxHeads];
  vector <vector <int> > tempGroup;
  bool convergedFlag = false;
  bool sameHeadFlag = true;
  while(sameHeadFlag)
  {
    sameHeadFlag = false;
    convergedFlag = false;
    //Clear before added members
    if (retryTimes>(m_numNodes-m_numMaxHeads+1)) {
      return false;
    }
    for (unsigned  int i=0 ; i<tempGroup.size(); i++)tempGroup[i].clear(); 
    //clear all the eixsted group members
    tempGroup.clear();
    for (int i=0; i<m_numMaxHeads; i++)
    {
      vector <int> tempV;
      tempGroup.push_back(tempV);
      tempHeadX[i] = nodes[i+retryTimes].locX;
      tempHeadY[i] = nodes[i+retryTimes].locY;
      tempHeadList[i]= nodes[i+retryTimes].nodeIndex;
    }
    while(!convergedFlag) // This loop want to find a new K-means coordinate
    {
      for (unsigned int i=0 ; i<tempGroup.size(); i++) tempGroup[i].clear(); 
      for (int i=0; i<m_numNodes; i++)
      {
        double closetDistance = std::numeric_limits<double>::max( );
        int closetHeadIndex = -1;
        for(int j=0; j<m_numMaxHeads; j++)
        {
          float tempDistance = (tempHeadX[j] - nodes[i].locX)*(tempHeadX[j] - nodes[i].locX)
                               +(tempHeadY[j] - nodes[i].locY)*(tempHeadY[j] - nodes[i].locY);
          if (closetDistance > tempDistance )
          {
            closetDistance = tempDistance;
            closetHeadIndex = j;
          }
        }
        tempGroup[closetHeadIndex].push_back(i);
      }
      convergedFlag = true;
      //find the k-means coordinate of each cluster
      for (int i = 0; i < tempGroup.size(); i++) {
        cout << "cluster: " << i <<"-th ";
        for (int j = 0; j < tempGroup[i].size(); j++) {
          cout << tempGroup[i][j] << ' ';
        }
        cout << endl;
      }
      for(int i=0; i<m_numMaxHeads; i++)
      {
        float newHx = 0;
        float newHy = 0;
        arma::vec vecDistance = arma::zeros<arma::vec>(tempGroup[i].size());
        for(unsigned int j=0; j<tempGroup[i].size(); j++)
        {
          float tempDistance = 0.0;
          for (int k = 0; k < tempGroup[i].size(); k++) {
            if ( j == k ) continue;
            tempDistance += 
              sqrt ( (nodes[tempGroup[i][j]].locX - nodes[tempGroup[i][k]].locX ) * 
                  (nodes[tempGroup[i][j]].locX - nodes[tempGroup[i][k]].locX) + 
                  (nodes[tempGroup[i][j]].locY - nodes[tempGroup[i][k]].locY) * 
                  (nodes[tempGroup[i][j]].locY - nodes[tempGroup[i][k]].locY) ) ;
          }
          vecDistance.at(j) = tempDistance; 
        }
        arma::uword idx;
        vecDistance.min(idx);
        newHx = nodes[tempGroup[i][idx]].locX;
        newHy = nodes[tempGroup[i][idx]].locY;
        if ( (abs(newHx-tempHeadX[i]) > 0.01) || (abs(newHy-tempHeadY[i])>0.01) ) convergedFlag = false; 
        // checkcheck if the original head close enough
        //find the new approriate location of the head
        tempHeadX[i] = newHx;
        tempHeadY[i] = newHy;
        tempHeadList[i] = tempGroup[i][idx];
      }
    }
    for (int i=0; i<m_numMaxHeads; i++)
      for (int j=i+1; j<m_numMaxHeads; j++)
        if (tempHeadList[i] == tempHeadList[j]) sameHeadFlag = true;
    retryTimes++;
  }

  return true;

}
