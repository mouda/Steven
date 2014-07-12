#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <limits>
#include <algorithm>//for sort
#include <cfloat>
#include <list>
#include <cassert>
#include <armadillo>
#include <iterator>

#define PENALTYSIZE 10e-6
using namespace std;
#include "minPowerSACluster.h"
bool pairCompare(const std::pair<int,double>& lPair, const std::pair<int,double>& rPair)
{
  return lPair.second > rPair.second;
}

MinPowerSACluster::MinPowerSACluster(
    FILE *fileReadCursor, 
    int inputTotalNodes, 
    int inputMaxChNum,
    int inputSAFac,  
    int inOutputControl,
    int isStrucOuput,
    double inputTemprature, 
    double InputSaAlpha, 
    double inCorrelationFactor, 
    std::string ipAddr, 
    Map const * const myPtrMap,
    CORRE_MA_OPE const * const myPtrGField,
    double tier1TxTime, 
    int tier2NumSlot,
    bool logFlag):
  m_ptrMap(myPtrMap),
  m_tier1TxTime(tier1TxTime),
  m_tier2NumSlot(tier2NumSlot),
  matrixComputer(myPtrGField),
  m_logFlag(logFlag)
{
  sysComputing = new SimSystem;
  terminated=false;
  //----------------------//
  //Simulation Control    //`
  //----------------------//

  SAIter = inputSAFac;
  temparature  = inputTemprature;
  constantIniTemprature=inputTemprature;
  alpha = InputSaAlpha;
  outCtrl = inOutputControl;

  isDetailOutputOn = isStrucOuput;
  strIpAddr = ipAddr;
  //----------------------//
  //System Parameter      //
  //----------------------//
  totalNodes = inputTotalNodes;
  maxChNum = inputMaxChNum;
  //correlationFactor = inCorrelationFactor;
  compRatio=inCorrelationFactor;//temporary name
  bestFeasibleJEntropy=0;
  bestFeasibleSupNum=0;
  float dx, dy;
  int inputIndex = 0;
  ULAGENT inputNode;
  //*********************//
  //Read Topology        //
  //*********************//
  while(inputIndex<totalNodes) {
    if( fscanf(fileReadCursor, "%f", &dx)<0 || fscanf(fileReadCursor, "%f", &dy)<0 ) {
      cout << "File map's content is wrong!!\n";
      assert(0);
    }
    else {
      inputNode.aryConstructor(inputIndex,myPtrMap->GetNodeXPos(inputIndex), myPtrMap->GetNodeYPos(inputIndex));
      nodes.push_back(inputNode);
    }
    inputIndex++;
  }

  //********************************//
  //Compute Gij distanceOf2Nodes Gib//
  //********************************//
  Gib = new float   [totalNodes];
  rateibMax = new double [totalNodes];
  distanceOf2Nodes = new double* [totalNodes];
  Gij = new float* [totalNodes];
  for(int i=0; i<totalNodes; i++) {
    distanceOf2Nodes[i] = new double [totalNodes];
    Gij[i] = new float [totalNodes];
    Gib[i] = sysComputing->returnChannelGain_BS(nodes[i]);
    for(int j = 0; j<totalNodes; j++) {
      if (i==j) {
        distanceOf2Nodes[i][j]=0.0;
        Gij[i][j] = 1.0;
      }
      else if (i>j) {
        distanceOf2Nodes[i][j]=distanceOf2Nodes[j][i];
        Gij[i][j] = Gij[j][i];
      }
      else {
        Gij[i][j] = sysComputing->returnChannelGain_2Nodes(nodes[i],nodes[j]);
        double tempDib = pow(nodes[i].locX-nodes[j].locX,2) +pow(nodes[i].locY - nodes[j].locY,2);
        distanceOf2Nodes[i][j] = tempDib;
      }
    }
  }

  cSystem = new ULCS1b(inputTotalNodes, inputMaxChNum);
  powerBest = new double [totalNodes];
  nextNodePower = new double [totalNodes];
  for (int i=0; i<maxChNum; i++ )
    for(int j=0; j<totalNodes; j++)
      cSystem->clusterStru[i][j] = false;
  sortIndex = -1;//Just a check point
  //Keep the sturcture for best Structure
  aryFlagHRDone = new bool [maxChNum];
  for(int i=0; i<maxChNum; i++)aryFlagHRDone[i]=false;

  //-----Best Performance Index-----//

  //Interference Matrix Keep the interference index/Name set from othe cluster for each cluster.
  bestMaClusterStru = new bool *[maxChNum];
  for (int i=0; i<maxChNum; i++) {
    bestMaClusterStru[i]= new bool [totalNodes];
  }

  bestAllSupStru = new bool[totalNodes];
  prevAllSupStru = new bool[totalNodes];
  vecBestBpshz_forVerification.resize(totalNodes);
  vecBestSINR_forVerification.resize(totalNodes);
  vecChooseIndex.reserve(totalNodes);


  maBestInterference= new double* [maxChNum];
  for (int i=0; i<maxChNum; i++)
    maBestInterference[i] = new double[maxChNum];

  maBestInterfernceIndex= new int* [maxChNum];
  for (int i=0; i<maxChNum; i++)
    maBestInterfernceIndex[i] = new int[maxChNum];

  bestFeasiblePayoff = DBL_MAX;

}
MinPowerSACluster::~MinPowerSACluster()
{
    if(!terminated)releaseMemory();
    m_logFile.close();
}

std::vector<int>
MinPowerSACluster::GetAllSupStru() const {
  std::vector<int> myAllSupStru(totalNodes,0);
  std::list<std::list<int> >::const_iterator iterRow = listCluMemBest.begin();
  for (;iterRow != listCluMemBest.end(); ++iterRow) {
    std::list<int>::const_iterator iterCol = iterRow->begin();
    if (*iterCol > 0 ) {
      for (; iterCol != iterRow->end(); ++iterCol) {
        myAllSupStru.at(*iterCol) = 1;
      }
    }
  }
  return myAllSupStru;
}

void MinPowerSACluster::releaseMemory()
{

    //release distanceOf2Nodes
    for(int i=0; i<totalNodes; i++) delete [] distanceOf2Nodes[i];
    delete [] distanceOf2Nodes;
    //release Gij
    for(int i=0; i<totalNodes; i++) delete [] Gij[i];
    delete [] Gij;
    //release Gb
    delete [] Gib;
    delete [] rateibMax;
    delete [] bestAllSupStru;
    delete [] prevAllSupStru;
    //release bestMaClusterStru
    for(int i=0; i<maxChNum; i++) delete [] bestMaClusterStru[i];
    delete [] bestMaClusterStru;

    //release matrix of interference of best peformance
    for(int i=0; i<maxChNum; i++) delete [] maBestInterference[i];
    delete [] maBestInterference;
    for(int i=0; i<maxChNum; i++) delete [] maBestInterfernceIndex[i];
    delete [] maBestInterfernceIndex;


    delete cSystem;
    //delete  cSystem and  share same pointer we ony need to delete once.
    delete [] powerBest;
    delete [] nextNodePower;
    delete [] aryFlagHRDone;
    delete sysComputing;
    terminated=true;
}


//-------------------------------------------------------------------//
// @Purpose: Read Topology File and Calculate Gij,distanceOf2Nodes
// @Called: by main
//-------------------------------------------------------------------//
bool MinPowerSACluster::setSystem(float inPowerMaxWatt, int inQuantizationBits,double inBandwidthKhz, double inFidelity)
{
    powerMax = inPowerMaxWatt;
    quantizationBits = inQuantizationBits;
    bandwidthKhz = inBandwidthKhz;
    realNoise = sysComputing->returnInBandThermalNoise(bandwidthKhz);
    fidelityRatio=inFidelity;
    power1st=powerMax;
     for (int i=0; i<totalNodes; i++) {
        rateibMax[i] = sysComputing->returnRate_BS(nodes[i],bandwidthKhz,powerMax);//bps
    }
    return true;
}

/*
    Purpose:coonstruct a initialized Clustering Stucture
    Flag: "kemans",..,
    and set the initial interference std::list for each node
*/
bool MinPowerSACluster::setInitialStucture(char* iniFlag)
{
    iniDone=false;
    bool normalFlag = true;
    cSystem->resetSystem();
    resetSA3iSystem();
    for(int i=0; i<totalNodes; i++) {
        nodes[i].power=0;
        nodes[i].ptrHead=NULL;
    }


    if (!strcmp(iniFlag, "kmeans")) 
      normalFlag = setIniStruKmeans();
    else if (!strcmp(iniFlag, "kmedoids_distance")) 
      normalFlag = setIniStruDistanceKmedoids();
    else if (!strcmp(iniFlag, "GraphPartition"))
      normalFlag = setIniGraphPartition();
    else if (!strcmp(iniFlag, "BalancedModelCluster"))
      normalFlag = setIniBanancedModelCluster();


    //--------------------------------------------------------------//
    //Check the initial constraint and drop node form the furtherest//
    //--------------------------------------------------------------//
    curSupNum = totalNodes;
    iniDone = true;
    return normalFlag;
}

/*
    @Purpose: Set the initial clustering stucture by kmeans algorithm

*/
bool MinPowerSACluster::setIniStruKmeans()
{
  int retryTimes = 0;
  float tempHeadX [maxChNum];
  float tempHeadY [maxChNum];
  int tempHeadList [maxChNum];
  std::vector <std::vector <int> > tempGroup;
  bool convergedFlag = false;
  bool sameHeadFlag = true;
  while(sameHeadFlag)
  {
    sameHeadFlag = false;
    convergedFlag = false;
    //Clear before added members
    if (retryTimes>(totalNodes-maxChNum+1))
    {
      return false;
    }
    for (unsigned  int i=0 ; i<tempGroup.size(); i++)tempGroup[i].clear(); //clear all the eixsted group members
    tempGroup.clear();

    for (int i=0; i<maxChNum; i++)
    {
      std::vector <int> tempV;
      tempGroup.push_back(tempV);
      tempHeadX[i] = nodes[i+retryTimes].locX;
      tempHeadY[i] = nodes[i+retryTimes].locY;
      tempHeadList[i]= nodes[i+retryTimes].nodeIndex;
    }
    while(!convergedFlag) // This loop want to find a new K-means coordinate
    {
      //choose "maxChNum" nember of node to use as the intial Cluster head
      //Find the closet head to form a cluster

      //Same number cluster head but clear in the converge process
      for (unsigned int i=0 ; i<tempGroup.size(); i++)tempGroup[i].clear(); //clear all the eixsted group members
      for (int i=0; i<totalNodes; i++)
      {
        float closetDistance = numeric_limits<float>::max( );
        int closetHeadIndex = -1;
        for(int j=0; j<maxChNum; j++)
        {
          float tempDistance = (tempHeadX[j] - nodes[i].locX)*(tempHeadX[j] - nodes[i].locX)\
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
      for(int i=0; i<maxChNum; i++)
      {
        float newHx = 0;
        float newHy = 0;
        for(unsigned int j=0; j<tempGroup[i].size(); j++)
        {
          newHx += nodes[tempGroup[i][j]].locX;
          newHy += nodes[tempGroup[i][j]].locY;
        }
        newHx /= tempGroup[i].size();
        newHy /= tempGroup[i].size();
        if ( (abs(newHx-tempHeadX[i]) > 0.01) || (abs(newHy-tempHeadY[i])>0.01) ) convergedFlag = false; // checkcheck if the original head close enough
        //find the new approriate location of the head
        tempHeadX[i] = newHx;
        tempHeadY[i] = newHy;
      }
    }
    // leave this loop if 'converged = 1;'
    for (int i=0; i<maxChNum; i++) tempHeadList[i] = \
      tempGroup[i][returnClosetNodeIndexInGroup(tempHeadX[i], tempHeadY[i], tempGroup[i])];
    //check there is same head exist
    for (int i=0; i<maxChNum; i++)
      for (int j=i+1; j<maxChNum; j++)if (tempHeadList[i] == tempHeadList[j]) sameHeadFlag = true;

    retryTimes++;
  }
  //Construct the initial Structure in the Group from
  //tempHeadList[] tempGroup[][]

  //---------------------------------------------------------
  //Intial Structure: Full Set and connection set according to K-means

  for (int i=0; i<maxChNum; i++)
  {
    cSystem->addNewHeadCs(tempHeadList[i]);
    for(unsigned int j=0 ; j<tempGroup[i].size(); j++)
    {
      addMemberSAIni(i, tempGroup[i][j]);//we correct the ptrHead later, Because the address will change
    }
  }
  //Re assign the ptrHead NOW
  for (int i=0; i<maxChNum; i++)
  {
    for(int j=0; j<totalNodes; j++)
    {
      if(cSystem->clusterStru[i][j]==true)
        nodes[j].ptrHead = &(cSystem->vecHeadName[i]);
    }
  }
  return true;
}



bool MinPowerSACluster::setIniStruDistanceKmedoids() 
{
  int retryTimes = 0;
  float tempHeadX [maxChNum];
  float tempHeadY [maxChNum];
  int tempHeadList [maxChNum];
  std::vector <std::vector <int> > tempGroup;
  bool convergedFlag = false;
  bool sameHeadFlag = true;
  while(sameHeadFlag)
  {
    sameHeadFlag = false;
    convergedFlag = false;
    //Clear before added members
    if (retryTimes>(totalNodes-maxChNum+1))
    {
      return false;
    }
    for (unsigned  int i=0 ; i<tempGroup.size(); i++)tempGroup[i].clear(); //clear all the eixsted group members
    tempGroup.clear();

    for (int i=0; i<maxChNum; i++)
    {
      std::vector <int> tempV;
      tempGroup.push_back(tempV);
      tempHeadX[i] = nodes[i+retryTimes].locX;
      tempHeadY[i] = nodes[i+retryTimes].locY;
      tempHeadList[i]= nodes[i+retryTimes].nodeIndex;
    }
    while(!convergedFlag) // This loop want to find a new K-means coordinate
    {
      //choose "maxChNum" nember of node to use as the intial Cluster head
      //Find the closet head to form a cluster

      //Same number cluster head but clear in the converge process
      for (unsigned int i=0 ; i<tempGroup.size(); i++) tempGroup[i].clear(); //clear all the eixsted group members
      for (int i=0; i<totalNodes; i++)
      {
        float closetDistance = numeric_limits<float>::max( );
        int closetHeadIndex = -1;
        for(int j=0; j<maxChNum; j++)
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
//#ifdef DEBUG
      //find the k-means coordinate of each cluster
      for (int i = 0; i < tempGroup.size(); i++) {
        cout << "cluster: " << i <<"-th ";
        for (int j = 0; j < tempGroup[i].size(); j++) {
          cout << tempGroup[i][j] << ' ';
        }
        cout << endl;
      }
//#endif
      for(int i=0; i<maxChNum; i++)
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
        if ( (abs(newHx-tempHeadX[i]) > 0.01) || (abs(newHy-tempHeadY[i])>0.01) ) convergedFlag = false; // checkcheck if the original head close enough
        //find the new approriate location of the head
        tempHeadX[i] = newHx;
        tempHeadY[i] = newHy;
        tempHeadList[i] = tempGroup[i][idx];
      }
    }
    // leave this loop if 'converged = 1;'
//    for (int i=0; i<maxChNum; i++) tempHeadList[i] = \
//      tempGroup[i][returnClosetNodeIndexInGroup(tempHeadX[i], tempHeadY[i], tempGroup[i])];
    //check there is same head exist
    for (int i=0; i<maxChNum; i++)
      for (int j=i+1; j<maxChNum; j++)if (tempHeadList[i] == tempHeadList[j]) sameHeadFlag = true;

    retryTimes++;
  }
// -------------------------------------------------------------------------- //
// @Description: confirm initialization of structure
// @Provides: 
// -------------------------------------------------------------------------- //
  for (int i=0; i<maxChNum; i++)
  {
    cSystem->addNewHeadCs(tempHeadList[i]);
    for(unsigned int j=0 ; j<tempGroup[i].size(); j++)
    {
      addMemberSAIni(i, tempGroup[i][j]);//we correct the ptrHead later, Because the address will change
    }
  }
  //Re assign the ptrHead NOW
  for (int i=0; i<maxChNum; i++)
  {
    for(int j=0; j<totalNodes; j++)
    {
      if(cSystem->clusterStru[i][j]==true)
        nodes[j].ptrHead = &(cSystem->vecHeadName[i]);
    }
  }

  return true;
}

bool
MinPowerSACluster::setIniGraphPartition()
{
  /* read file */
  fstream fpart;
  fpart.open("WeightMatrix_N100.out_part_25",std::ios::in);
  std::vector<std::vector<int> > tmpGroup(maxChNum);
  int nodeName = 0;
  while(nodeName < m_ptrMap->GetNumNodes()) {
    int k;
    fpart >> k;
    tmpGroup.at(k).push_back(nodeName);
    ++nodeName;
  }
  /* determine the initial head */
  std::vector<int> tmpVecHead(maxChNum);
  std::vector<std::vector<int> >::const_iterator constIterRow = tmpGroup.begin();
  int clusterIdx = 0;
  for (; constIterRow != tmpGroup.end(); ++constIterRow, ++clusterIdx) {
    std::vector<int>::const_iterator constIterCol = constIterRow->begin();
    double  maxGain = -DBL_MAX;
    int     maxGainIdx = -1;
    for (; constIterCol != constIterRow->end(); ++constIterCol) {
      if (maxGain < m_ptrMap->GetGi0ByNode(*constIterCol)) {
        maxGain = m_ptrMap->GetGi0ByNode(*constIterCol);
        maxGainIdx = *constIterCol;
      }
    }
    assert(maxGainIdx != -1 );
    tmpVecHead.at(clusterIdx) = maxGainIdx;
  }

  for (int i = 0; i < maxChNum; ++i) {
    cSystem->addNewHeadCs(tmpVecHead.at(i));
    for(unsigned int j=0 ; j < tmpGroup.at(i).size(); ++j) {
      addMemberSAIni(i, tmpGroup[i][j]);//we correct the ptrHead later, Because the address will change
    }
  }
  //Re assign the ptrHead NOW
  for (int i=0; i<maxChNum; i++) {
    for(int j=0; j<totalNodes; j++) {
      if(cSystem->clusterStru[i][j]==true)
        nodes[j].ptrHead = &(cSystem->vecHeadName[i]);
    }
  }
  return true;

}

bool
MinPowerSACluster::setIniBanancedModelCluster()
{
  int retryTimes = 0;
  float tempHeadX [maxChNum];
  float tempHeadY [maxChNum];
  int tempHeadList [maxChNum];
  std::vector <std::vector <int> > tempGroup;
  bool convergedFlag = false;
  bool sameHeadFlag = true;
  while(sameHeadFlag) {
    sameHeadFlag = false;
    convergedFlag = false;
    //Clear before added members
    if (retryTimes>(totalNodes-maxChNum+1)) {
      return false;
    }
    for (unsigned  int i=0 ; i<tempGroup.size(); i++)tempGroup[i].clear(); //clear all the eixsted group members
    tempGroup.clear();

    for (int i=0; i<maxChNum; i++) {
      std::vector <int> tempV;
      tempGroup.push_back(tempV);
      tempHeadX[i] = nodes[i+retryTimes].locX;
      tempHeadY[i] = nodes[i+retryTimes].locY;
      tempHeadList[i]= nodes[i+retryTimes].nodeIndex;
    }
    while(!convergedFlag) { // This loop want to find a new K-means coordinate 
      //choose "maxChNum" nember of node to use as the intial Cluster head
      //Find the closet head to form a cluster

      //Same number cluster head but clear in the converge process
      for (unsigned int i=0 ; i<tempGroup.size(); i++) tempGroup.at(i).clear(); //clear all the eixsted group members
      std::vector<bool> vecAssigned(totalNodes, false);
      for (int j = 0; j < maxChNum - 1; ++j) {
        std::vector<std::pair<int,double> > vecDVi;
        for (int i = 0; i < totalNodes; ++i) {
          if (vecAssigned.at(i) ) continue; 
          double minDistance = DBL_MAX;
          int    minDisChIdx = -1;
          for (int k = j+1; k < maxChNum; ++k) {
            if (GetNodeDistance(k,i) < minDistance) {
              minDistance = GetNodeDistance(k,i);
              minDisChIdx = k;
            }
          }
          double dVi = minDistance - GetNodeDistance(j,i);
          vecDVi.push_back(std::make_pair(i, dVi));
        }
        std::sort(vecDVi.begin(), vecDVi.end(), pairCompare);
        for (int l = 0; l <= m_tier2NumSlot; ++l) {
          tempGroup.at(j).push_back(vecDVi.at(l).first);
          vecAssigned.at(vecDVi.at(l).first) = true;
        }
      }
      for (int i = 0; i < totalNodes; ++i) {
        if (vecAssigned.at(i) != true) {
          tempGroup.at(maxChNum-1).push_back(i);
        }
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
      for(int i=0; i<maxChNum; i++) {
        float newHx = 0;
        float newHy = 0;
        arma::vec vecDistance = arma::zeros<arma::vec>(tempGroup[i].size());
        for(unsigned int j=0; j<tempGroup[i].size(); j++) {
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
        if ( (abs(newHx-tempHeadX[i]) > 0.01) || (abs(newHy-tempHeadY[i])>0.01) ) convergedFlag = false; // checkcheck if the original head close enough
        //find the new approriate location of the head
        tempHeadX[i] = newHx;
        tempHeadY[i] = newHy;
        tempHeadList[i] = tempGroup[i][idx];
      }
    }
    for (int i=0; i<maxChNum; i++)
      for (int j=i+1; j<maxChNum; j++)if (tempHeadList[i] == tempHeadList[j]) sameHeadFlag = true;

    retryTimes++;
  }
  // -------------------------------------------------------------------------- //
  // @Description: confirm initialization of structure
  // @Provides: 
  // -------------------------------------------------------------------------- //
  for (int i=0; i<maxChNum; i++) {
    cSystem->addNewHeadCs(tempHeadList[i]);
    for(unsigned int j=0 ; j<tempGroup[i].size(); j++) {
      addMemberSAIni(i, tempGroup[i][j]);//we correct the ptrHead later, Because the address will change
    }
  }
  //Re assign the ptrHead NOW
  for (int i=0; i<maxChNum; i++) {
    for(int j=0; j<totalNodes; j++)
    {
      if(cSystem->clusterStru[i][j]==true)
        nodes[j].ptrHead = &(cSystem->vecHeadName[i]);
    }
  }

  return true;
}

double
MinPowerSACluster::GetNodeDistance( const int lName, const int rName)
{
  return sqrt( 
      (nodes[lName].locX - nodes[rName].locX) * (nodes[lName].locX - nodes[rName].locX) +
      (nodes[lName].locY - nodes[rName].locY) * (nodes[lName].locY - nodes[rName].locY)
      );
}

double 
MinPowerSACluster::returnComprRatio()
{

    bool tempAry2[totalNodes];
    for(int i=0; i<totalNodes; i++)tempAry2[i]=true;
    double sysRedundancyLocal=matrixComputer->computeLog2Det(1.0,tempAry2);
    double totalEntropyLocal =totalNodes*indEntropy;
    double compressRatio = 1-(sysRedundancyLocal+totalEntropyLocal)/totalEntropyLocal;
    //cout<<"done"<<endl;
    return compressRatio;
}

bool MinPowerSACluster::startCool()
{
  begin = clock();
  indEntropy = m_ptrMap->GetIdtEntropy();
  tempAddT=0;
  tempDisT=0;
  tempHRT=0;
  bool inClu[totalNodes];
  for(int i=0; i<totalNodes; i++)inClu[i]=true;
  double sysRedundancy =matrixComputer->computeLog2Det(1.0, inClu);
  wholeSystemEntopy = totalNodes*indEntropy+sysRedundancy;

  flagAnsFound  =  false;
  //-----------------------------//
  //Initialize performance matrix//
  //-----------------------------//

  curSupNum     =   cSystem->calSupNodes();
  curChNum      =   maxChNum;
  nextChNum     =   curChNum;
  curJEntropy   =   curSupNum*indEntropy + matrixComputer->computeLog2Det(1.0,cSystem->allSupStru);

  m_cur1st_watt   = OptimalRateControl();

  vector<double>  sizePenalty(m_ptrMap->GetNumInitHeads(), 10.0);
  vector<double>  tier2Penalty(m_ptrMap->GetNumNodes(), 1.0);
  double          entropyPenalty = 10.0;
  vector<double>  tmpSizePenalty(m_ptrMap->GetNumInitHeads(), 10.0);
  vector<double>  tmpTier2Penalty(m_ptrMap->GetNumNodes(), 10.0);
  double          tmpEntropyPenalty = 10.0;

  m_curPayoff     = GetPayOff(sizePenalty, tier2Penalty, entropyPenalty);
  m_prevVecClusterSize.assign(cSystem->vecClusterSize.begin(), cSystem->vecClusterSize.end());
  m_prevVecHeadName.assign(cSystem->vecHeadName.begin(), cSystem->vecHeadName.end());
  cur2nd_Joule    = returnTransientJoule();
  cur1st_Joule    = power1st*cur1st_ms/1000.0;
  m_logFile.open("minPoweriLog.out",std::ios::out);

  cout << "m_curPayoff: " << m_curPayoff << endl;
  bestAllServeFound = false;

  if ( checkBestClusterStructure_DataCentric(0) ) return true;
  cout << "Compression Ratio " << returnComprRatio() << " indEntropy " << indEntropy << endl;

  for(int i = 1; i < SAIter; ++i)
  {
    if (i%2 == 0) {
      GetNeighbor1(i);
      bool curAllServe = (curJEntropy>(fidelityRatio*wholeSystemEntopy)?true:false) && CheckTier2Feasible() ; 
      for (int m = 0; m < m_prevVecClusterSize.size(); ++m) {
        if (m_prevVecHeadName.at(m) >= 0 && m_prevVecClusterSize.at(m) > m_tier2NumSlot + 1) {
          curAllServe = false;
        }
      }
      if (m_logFlag) {
        m_logFile << curAllServe << ' ' << m_curPayoff <<' ' << bestFeasiblePayoff << endl;
      }

      if( targetHeadIndex == -1 || targetHeadIndex == -1 ) {
        if ( nextEventFlag == 4 ) {
          passNext2Cur();
        }
        continue;
      }

      calculateMatrics_minResors(sizePenalty, tier2Penalty, entropyPenalty);
      ConfirmNeighbor1();
    }
    else {
      GetNeighbor2(i, sizePenalty, tier2Penalty, entropyPenalty, tmpSizePenalty, tmpTier2Penalty, tmpEntropyPenalty );
      ConfirmNeighbor2(sizePenalty, tier2Penalty, entropyPenalty, tmpSizePenalty, tmpTier2Penalty, tmpEntropyPenalty );
    }

    if ( curJEntropy > ( fidelityRatio * wholeSystemEntopy ) ) {
      flagAnsFound=true;
    }
    assert(curSupNum>=0);
    if(checkBestClusterStructure_DataCentric(i)) {
      cout<<"Congratulation All nodes are served"<<endl;
      break;
    }
    int tempche = (signed)cSystem->listUnSupport->size();

    assert(( curSupNum + tempche) <= totalNodes);
    temparature*=alpha;

  }

  end = clock();
  computingTimes = ((float)(end-begin))/CLOCKS_PER_SEC;
  cout << "best "<< bestFeasibleSupNum << "   Information Ratio:" << bestFeasibleJEntropy / wholeSystemEntopy << endl;
  if(!flagAnsFound) {
    cout<<"Not Found the answer Yet"<<endl;
    //return false ; // We don't care if the result is feasible.
  }
  else cout<<"SA end up correctly"<<endl;
  cout<<"Add Times="<<tempAddT<<endl;
  cout<<"Discard Times="<<tempDisT<<endl;
  cout<<"HR Times="<<tempHRT<<endl;

  char timeBuf[32];
  TimeStamp obj_time;
  obj_time.returnRealWordTime(timeBuf,32);

  char str[500];
  char str2[500];
  char str3[500];
  if( bestFeasibleJEntropy >= ( wholeSystemEntopy * fidelityRatio) ) {
    cout<<timeBuf<<endl;
    return false;
  }
}

double
MinPowerSACluster::
OptimalRateControl() const
{
  Index numVariables = nextChNum;
  Index numConstraints = 1;
  Index numNz_jac_g = numVariables;
  Index numNz_h_lag = numVariables;

  // Create an instance of your nlp...
  SmartPtr<MyTier1NLP> mynlp = 
    new MyTier1NLP( numVariables, numConstraints, numNz_jac_g, numNz_h_lag,
        m_ptrMap,
        cSystem,
        matrixComputer,
        m_tier1TxTime
        );
  MyTier1NLP* rawPtr = dynamic_cast<MyTier1NLP*>(GetRawPtr(mynlp));

  // Create an instance of the IpoptApplication
  //
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int) status;
  }

  status = app->OptimizeTNLP(mynlp);

  if (status == Solve_Succeeded) {
    // Retrieve some statistics about the solve
    Index iter_count = app->Statistics()->IterationCount();
    //std::cout << std::endl << std::endl << "*** The problem solved in " << iter_count << " iterations!" << std::endl;

    Number final_obj = app->Statistics()->FinalObjective();
    //std::cout << std::endl << std::endl << "*** The final value of the objective function is " << final_obj << '.' << std::endl;
  }
  return rawPtr->GetMinimalPower();
}


/*
    Movement Function
    use in setInitialStructure
*/
void MinPowerSACluster::addMemberSAIni(int inputHeadIndex, int inputMemberName)
{
    cSystem->addMemberCs(inputHeadIndex,inputMemberName,iniDone);
    nodes[inputMemberName].ptrHead = cSystem->returnHeadPtr(inputHeadIndex);
    nextSupNum = curSupNum + 1;//Serve one more node
}

/*
    do one step movement add/discard/
*/
void MinPowerSACluster::GetNeighbor1( const int iterSA)
{
  int probAdd = 0;
  int probDiscard = 0;
  if (iterSA < SAIter/5) {
    probAdd = ((curSupNum<(totalNodes)) ?20000 :0);
    probDiscard = ((curSupNum<(maxChNum+1)) ?0:25000);
  }
  else {
    probAdd = ((curSupNum<(totalNodes)) ?20000 :0);
    probDiscard = ((curSupNum<(maxChNum+1)) ?0:20000);
  }
//  int probAdd = 0;
//  int probDiscard = 0;
  int probExchange = 0;
  bool checkRotateble=false;//check if there are rotatableSet;
  for(int i=0; i<maxChNum; i++) {
    //cout<<aryFlagHRDone[i]<<" "<<cSystem->vecClusterSize[i]<<endl;
    if( aryFlagHRDone[i] == false && cSystem->vecClusterSize[i] > 1 )
      checkRotateble=true;
  }
  int probHeadRotate = ((checkRotateble) ?1000:0);//Don't do head rotate if there are only a few nodes
  //int probHeadRotate = 0;



  int tmpJoinCan=0;
  bool chkLessCluster=cSystem->returnIfClusterSmall(thresholdd,tmpJoinCan);

  //int probJoin = (chkLessCluster&&curJEntropy>(fidelityRatio*wholeSystemEntopy))?tmpJoinCan*50:0;
  int probJoin = 0;
  if (iterSA < SAIter*3/5) {
     probJoin = (chkLessCluster)?tmpJoinCan*125:0;
  }
  else {
     probJoin = (chkLessCluster)?tmpJoinCan*1250:0;
  } 
  //int probJoin = 0;

  //probJoin=((lastJoinPassAccu>thres2-400)?probJoin:0);

  int probIsoltae=((curChNum<maxChNum)?1250:0);
  //probIsoltae=((lastJoinPassAccu>thres2)?probIsoltae:0);
  //int probIsoltae = 0;

  if (iterSA < 10) {
    bestFeasiblePayoff = DBL_MAX;
  }

//  int sumProb = probAdd + probDiscard + probExchange + probHeadRotate+probJoin+probIsoltae;
  int sumProb = probAdd + probDiscard + probHeadRotate+probJoin+probIsoltae;
  int eventCursor= (int)((double)rand() / ((double)RAND_MAX + 1) * sumProb);
  nextEventFlag=-1;// this flag tell add or discard or Headrotate

  //-------------------------------------//
  //Decide event Flag                    //
  //-------------------------------------//

  if (eventCursor<probAdd) nextEventFlag = 1;
  else if (eventCursor<(probAdd+probDiscard)) nextEventFlag=2;
  else if (eventCursor<(probAdd+probDiscard+probHeadRotate)) nextEventFlag=3;
  else if (eventCursor<(probAdd+probDiscard+probHeadRotate+probJoin)) nextEventFlag=4;
  else if (eventCursor<sumProb)nextEventFlag=5;
  else
  {
    cout<<"Failure in the random step"<<endl;
    cout<<sumProb<<endl;
    cout<<eventCursor<<endl;
  }
  //-------------------------------------//
  // Start the movement                  //
  //-------------------------------------//  if (nextEventFlag==1)//Add
  if (nextEventFlag==1)
  {
    if (cSystem->listUnSupport->size()==0) cout<<"Error, it should haven't come in here with empty addlist and add."<<endl;
    else
    {
      //decideAdd3i_DC_HeadDetMemRan();
      decideAddRandSelectCluster();
      if(targetHeadIndex!=-1&&targetNode!=-1){
        addMemberSA(targetHeadIndex,targetNode);
      }
    }
    nextChNum=curChNum;

  }
  else if (nextEventFlag ==2)
  {
    //decideDiscard3b();
    decideDiscard3o();
    discardMemberSA(targetHeadIndex,targetNode);
    nextChNum=curChNum;
  }
  else if (nextEventFlag ==3)
  {
    decideHeadRotate2i_DC_HeadRanMemDet();
    rotateHeadSA(targetHeadIndex,targetNode);
    nextJEntropy = curJEntropy;
    nextSupNum = curSupNum;
    nextChNum=curChNum;

  }
  else if (nextEventFlag==4){
    JoiningHeadIndex=-1;
    decideHeadJoining4b();
    if (JoiningHeadIndex==-1 || targetHeadIndex==-1) {
      nextJEntropy = curJEntropy; // entropy unchanged
      nextSupNum = curSupNum; //support number unchanged
      return;
    }
    join_fromHeadSA(JoiningHeadIndex,targetHeadIndex);
    nextJEntropy = curJEntropy;
    nextSupNum = curSupNum;
  }
  else if (nextEventFlag==5){
    isolatedHeadIndex=-1;
    IsolateNodeName=-1;
    decideIsolate4b();
    isolateHeadSA(IsolateNodeName,isolatedHeadIndex,targetHeadIndex);
    nextJEntropy = curJEntropy;
    nextSupNum = curSupNum;
  }
  else if ( nextEventFlag == 6) {
    decideExchangeNode();
    cout << "Exchane: " << targetNode << endl;

  }
  else
  {
    cout<<"Error. The random Neighbor event "<<nextEventFlag<<" choose is wrong"<<endl;
    cout<<"CursupNum "<<curSupNum<<" maxChNUm "<<maxChNum<<endl;
    cout<<"SumProb "<<sumProb<<endl;
  }
}
void
MinPowerSACluster::GetNeighbor2( const int iterSA, 
    const vector<double>& sizePenalty, 
    const vector<double>& tier2Penalty, 
    const double& entropyPenalty,
    vector<double>& tmpSizePenalty,
    vector<double>& tmpTier2Penalty,
    double& tmpEntropyPenalty)
{
  /* Neighbor of Cluster Size */
  tmpSizePenalty.assign(sizePenalty.begin(), sizePenalty.end());
  tmpTier2Penalty.assign(tier2Penalty.begin(), tier2Penalty.end());
  tmpEntropyPenalty = entropyPenalty;
  for (int k = 0; k < cSystem->vecClusterSize.size(); ++k) {
    if (cSystem->vecClusterSize.at(k) != 0
        && (static_cast<double>(cSystem->vecClusterSize.at(k)) - static_cast<double>(m_tier2NumSlot) - 1.0) > 0.0 ) {
      if (sizePenalty.at(k) > 0 ) {
        if (rand()%2  == 1) {
          tmpSizePenalty.at(k) = sizePenalty.at(k) + 1 * PENALTYSIZE; 
        }
        else {
          tmpSizePenalty.at(k) = sizePenalty.at(k) - 1 * PENALTYSIZE; 
        }
      }
      else {
        tmpSizePenalty.at(k) = sizePenalty.at(k) + 1 * PENALTYSIZE; 
      } 
    }
    else {
      tmpSizePenalty.at(k) = sizePenalty.at(k);
    }
  }
  /* Neighbor of Link feasibity */
  list<list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();  
  for (int chIdx = 0; iterRow != cSystem->listCluMember->end(); ++iterRow, ++chIdx) {
    if (cSystem->vecHeadName.at(chIdx) == -1 ) continue; 
    list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      if ( GetTier2ExpectPower(*iterCol, cSystem->vecHeadName.at(chIdx)) > m_ptrMap->GetMaxPower()) {
        if (tier2Penalty.at(*iterCol) > 0 ) {
          if (rand()%2  == 1) {
            tmpTier2Penalty.at(*iterCol) = tier2Penalty.at(*iterCol) + 1 * PENALTYSIZE; 
          }
          else {
            tmpTier2Penalty.at(*iterCol) = tier2Penalty.at(*iterCol) - 1 * PENALTYSIZE; 
          }
        }
        else {
          tmpTier2Penalty.at(*iterCol) = tier2Penalty.at(*iterCol) + 1 * PENALTYSIZE;
        } 
      }
      else {
        tmpTier2Penalty.at(*iterCol) = tier2Penalty.at(*iterCol);
      }
    }
  }

  /* Neighbor of Entropy */
  double tmpJEntropy = nextSupNum*indEntropy + matrixComputer->computeLog2Det(1.0, cSystem->allSupStru);
  if (fidelityRatio*wholeSystemEntopy - tmpJEntropy > 0) {
    if (entropyPenalty > 0 ) {
      if (rand()%2  == 1) {
        tmpEntropyPenalty = entropyPenalty + 1 * PENALTYSIZE; 
      }
      else {
        tmpEntropyPenalty = entropyPenalty - 1 * PENALTYSIZE; 
      }
    }
    else {
      tmpEntropyPenalty = entropyPenalty + 1 * PENALTYSIZE;
    } 
  }
  else {
    tmpEntropyPenalty = entropyPenalty;
  }
}

double
MinPowerSACluster::GetPayOff( 
    const vector<double>& sizePenalty, 
    const vector<double>& tier2Penalty, 
    const double& entropyPenalty)
{
  cout << "S: " << GetSizePenalty(sizePenalty) << ", T: " << GetTier2Penalty(tier2Penalty) << ", E: " << GetEntropyPenalty(entropyPenalty) << endl;
  for (int i = 0; i < sizePenalty.size(); ++i) {
    cout << sizePenalty.at(i) << ' ';
  }
  cout << endl;
  return OptimalRateControl()+GetSizePenalty(sizePenalty) + GetTier2Penalty(tier2Penalty) + GetEntropyPenalty(entropyPenalty);
}


/*
 * to exchange the node randomly
 */
void 
MinPowerSACluster::decideExchangeNode()
{
  decideDiscard3b();

}

/*
    targetHeadIndex: CLOSET HEAD NOW. (Maybe The Head with least resource usage)
    targetNode: randomly proportional to the indepedent information compare to current set

*/
void 
MinPowerSACluster::decideAddSmallestSize() {
  targetHeadIndex=-1;
  targetNode=-1;
  double curInfo=curSupNum*indEntropy+matrixComputer->computeLog2Det(1.0,cSystem->allSupStru);
  double residualInformation=wholeSystemEntopy-curInfo;
  double randomCurs=((double)rand() / ((double)RAND_MAX + 1));
  double chooseCurs=0;
  std::vector<int>::const_iterator maxIter = std::min_element(cSystem->vecClusterSize.begin(), cSystem->vecClusterSize.end());
  int chIdx = 0;
  for (std::vector<int>::const_iterator iter = cSystem->vecClusterSize.begin(); iter != maxIter; ++iter) ++chIdx; 
  cout << "Min Ch idx: " << chIdx << "with size: " << cSystem->vecClusterSize.at(chIdx) <<endl;

  double maxGain = -DBL_MAX;
  int    maxName = -1; 
  list <int>::iterator itList=cSystem->listUnSupport->begin();

  for(; itList!=cSystem->listUnSupport->end(); ++itList) {
    double tmpGain = m_ptrMap->GetGijByPair(*itList, cSystem->vecHeadName.at(chIdx) );
    if (tmpGain > maxGain) {
      maxGain = tmpGain;
      maxName = *itList;
    }
  }
  /* Decide Result */
  targetNode = maxName;
  targetHeadIndex = chIdx;

}
void 
MinPowerSACluster::decideAdd3i_DC_HeadDetMemRan() {
  targetHeadIndex=-1;
  targetNode=-1;

  double curInfo=curSupNum*indEntropy+matrixComputer->computeLog2Det(1.0,cSystem->allSupStru);
  double residualInformation=wholeSystemEntopy-curInfo;
  double randomCurs=((double)rand() / ((double)RAND_MAX + 1));
  double chooseCurs=0;
  list <int>::iterator it_Int=cSystem->listUnSupport->begin();

  for(; it_Int!=cSystem->listUnSupport->end(); it_Int++) {
    cSystem->allSupStru[*it_Int]=true;
    double tmpIndInfo=(curSupNum+1)*indEntropy+matrixComputer->computeLog2Det(1.0,cSystem->allSupStru)-curInfo;
    chooseCurs+=tmpIndInfo;
    cSystem->allSupStru[*it_Int]=false;
    if ((chooseCurs/residualInformation)>randomCurs) {
      targetNode=*it_Int;
      break;
    }
  }
  if (targetNode==-1)return;

  double maxGain=0;
  //find closet head
  for(unsigned int i=0; i<cSystem->vecHeadName.size(); i++) {
    if (cSystem->vecHeadName.at(i) > 0 ) {
      if( cSystem->vecClusterSize.at(i) < m_tier2NumSlot + 1 && Gij[targetNode][cSystem->vecHeadName[i]] > maxGain) {
        maxGain=Gij[targetNode][cSystem->vecHeadName[i]];
        targetHeadIndex=i;
      }
    }
  }
}

void 
MinPowerSACluster::decideAddClosetAddableNode() {
  targetHeadIndex=-1;
  targetNode=-1;

  double maxGain = -DBL_MAX;
  int    maxName = -1; 
  int    chIdx = -1;
  for (int i = 0; i < cSystem->vecClusterSize.size(); ++i) {
    if (cSystem->vecHeadName.at(i) >= 0 && cSystem->vecClusterSize.at(i) < m_tier2NumSlot + 1) {
      std::list<int>::const_iterator itList = cSystem->listUnSupport->begin();
      for (; itList != cSystem->listUnSupport->end(); ++itList) {
        double tmpGain = m_ptrMap->GetGijByPair(*itList, cSystem->vecHeadName.at(i) );
        if (tmpGain > maxGain) {
          maxGain = tmpGain;
          maxName = *itList;
          chIdx = i;
        }
      }
    }
  }
  /* Decide Result */
  targetNode = maxName;
  targetHeadIndex = chIdx;

}

void
MinPowerSACluster::decideAddRandSelectCluster()
{
  targetHeadIndex = -1;
  targetNode = -1;
  list<list<int> >::iterator itli1 = cSystem->listCluMember->begin();
  int addableSize=0;
  for(; itli1 != cSystem->listCluMember->end(); itli1++) if(itli1->size()> 1 && itli1->size() < m_tier2NumSlot + 1)addableSize++;

  if (addableSize == 0) {
    return;
  }
  //Randomly choose a cluster
  itli1 = cSystem->listCluMember->begin();
  int chooseCursor = (int) ((double)rand() / ((double)RAND_MAX + 1) * addableSize)+1;//+1 Becasue the intial might
  assert(chooseCursor <= maxChNum && chooseCursor >=0 );
  itli1 = cSystem->listCluMember->begin();
  targetHeadIndex=0;
  for(int i=0; i<chooseCursor; itli1++,targetHeadIndex++) {
    if(itli1->size()> 1  && itli1->size() < m_tier2NumSlot + 1) ++i;
    if(i==chooseCursor)break;
  }
  double maxGain = -DBL_MAX;
  double tempInter = 0;
  int    maxGainName = -1;
  list <int>::iterator itList=cSystem->listUnSupport->begin();
  for(; itList!=cSystem->listUnSupport->end(); ++itList) {
    double tmpGain = m_ptrMap->GetGijByPair(*itList, cSystem->vecHeadName.at(targetHeadIndex) );
    if (tmpGain > maxGain) {
      maxGain = tmpGain;
      maxGainName = *itList;
    }
  }
  //cout << "Feasible: " << CheckLinkFeasible(cSystem->vecHeadName.at(targetHeadIndex), maxGainName) << endl;
//  if (!CheckLinkFeasible(cSystem->vecHeadName.at(targetHeadIndex), maxGainName)) {
//    maxGainName = -1;
//    targetHeadIndex = -1;
//  }

  targetNode = maxGainName;

}

/*
    targetHeadIndex: (Randomly) Choose a Head uniformly
    targetNode: (Deterministically) find The one cause strongest Interference to others
*/
void 
MinPowerSACluster::decideDiscard3b()
{
  //Compute the discardable size

  targetNode = -1;
  std::vector<int>::const_iterator minIter = std::max_element(cSystem->vecClusterSize.begin(), cSystem->vecClusterSize.end());
  int chIdx = 0;
  for (std::vector<int>::const_iterator iter = cSystem->vecClusterSize.begin(); iter != minIter; ++iter) ++chIdx; 
  std::list<std::list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();
  for (int i = 0; i < chIdx; ++i) { ++iterRow; }
  std::list<int>::const_iterator iterCol = iterRow->begin();
  double minGain = DBL_MAX;
  int    minName = -1; 
  for (; iterCol != iterRow->end(); ++iterCol) {
    double tmpGain = m_ptrMap->GetGijByPair(*iterCol, cSystem->vecHeadName.at(chIdx) );
    if (tmpGain < minGain) {
      minGain = tmpGain;
      minName = *iterCol; 
    }
  }

  /* Decide Result */
  targetNode = minName;
  targetHeadIndex = chIdx;
  assert(targetNode != -1 && targetHeadIndex != -1);
}
void 
MinPowerSACluster::decideDiscard3o()
{
    //Compute the discardable size
    list<list<int> >::iterator itli1 = cSystem->listCluMember->begin();
    int discardableSize=0;
    for(; itli1!=cSystem->listCluMember->end(); itli1++)if(itli1->size()>1)discardableSize++;
    assert(discardableSize!=0);
    //Randomly choose a cluster
    itli1 = cSystem->listCluMember->begin();
    targetHeadIndex = -1;
    int chooseCursor = (int) ((double)rand() / ((double)RAND_MAX + 1) * discardableSize)+1;//+1 Becasue the intial might
    assert(chooseCursor<=maxChNum&&chooseCursor>=0);

    itli1 = cSystem->listCluMember->begin();
    targetHeadIndex=0;
    for(int i=0; i<chooseCursor; itli1++,targetHeadIndex++)
    {
        if(itli1->size()>1)i++;
        if(i==chooseCursor)break;
    }

    assert(cSystem->vecHeadName[targetHeadIndex]!=-1);
    assert(cSystem->vecClusterSize[targetHeadIndex]==itli1->size());
    assert(cSystem->vecClusterSize[targetHeadIndex]>1);
    double minGain = DBL_MAX;
    double tempInter = 0;
    int    minGainName = -1;
    list<int>::iterator it2=itli1->begin();
    for(; it2!=itli1->end(); it2++) {
        tempInter = 0;
        if(cSystem->vecHeadName[targetHeadIndex]==(*it2))continue;
        else {
            if( m_ptrMap->GetGijByPair(cSystem->vecHeadName[targetHeadIndex], (*it2) ) < minGain) {
                minGain = m_ptrMap->GetGijByPair(cSystem->vecHeadName[targetHeadIndex], (*it2) ) ;
                minGainName = (*it2);
            }
        }
    }
    assert(minGainName!=-1);

    targetNode = minGainName;
}

void
MinPowerSACluster::decideDiscard3f()
{
    //Compute the discardable size
    list<list<int> >::iterator itli1 = cSystem->listCluMember->begin();
    int discardableSize=0;
    for(; itli1!=cSystem->listCluMember->end(); itli1++)if(itli1->size()>1)discardableSize++;
    assert(discardableSize!=0);
    //Randomly choose a cluster
    itli1 = cSystem->listCluMember->begin();
    targetHeadIndex = -1;
    int chooseCursor = (int) ((double)rand() / ((double)RAND_MAX + 1) * discardableSize)+1;//+1 Becasue the intial might
    assert(chooseCursor<=maxChNum&&chooseCursor>=0);

    itli1 = cSystem->listCluMember->begin();
    targetHeadIndex=0;
    for(int i=0; i<chooseCursor; itli1++,targetHeadIndex++)
    {
        if(itli1->size()>1)i++;
        if(i==chooseCursor)break;
    }

    assert(cSystem->vecHeadName[targetHeadIndex]!=-1);
    assert(cSystem->vecClusterSize[targetHeadIndex]==itli1->size());
    assert(cSystem->vecClusterSize[targetHeadIndex]>1);
    double minGain = DBL_MAX;
    double tempInter = 0;
    int    minGainName = -1;
    list<int>::iterator it2=itli1->begin();
    for(; it2!=itli1->end(); it2++) {
        tempInter = 0;
        if(cSystem->vecHeadName[targetHeadIndex]==(*it2))continue;
        else {
          if (!CheckLinkFeasible(cSystem->vecHeadName.at(targetHeadIndex), *it2)) {
            minGainName = *it2;
            break;
          }
        }
    }
    assert(minGainName!=-1);

    targetNode = minGainName;

}

void
MinPowerSACluster::decideDiscardMinGain()
{
  //Compute the discardable size

  targetNode = -1;
  targetHeadIndex = -1;

  double  minGain   = DBL_MAX;
  int     minChIdx  = -1;
  int     minNode   = -1;
  std::list<std::list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();
  for (int i = 0; iterRow != cSystem->listCluMember->end(); ++iterRow, ++i) {
    if (iterRow->size() > 1) {
      std::list<int>::const_iterator iterCol = iterRow->begin();
      for (; iterCol != iterRow->end(); ++iterCol) {
        double tmpGain = m_ptrMap->GetGijByPair(*iterCol,i);
        if (tmpGain < minGain) {
          minGain   = tmpGain;
          minChIdx  = i;
          minNode   = *iterCol; 
        }
      }
    }
  }

  /* Decide Result */
  targetNode = minNode;
  targetHeadIndex = minChIdx;
  assert(targetNode != -1 && targetHeadIndex != -1);
}

void 
MinPowerSACluster::decideHeadRotate2i_DC_HeadRanMemDet()
{
  //----Uniformly choosed
  int rotateAbleSize=0;
  for (int i=0; i<maxChNum; i++)
    if(cSystem->vecClusterSize[i]>1)rotateAbleSize++;
  //if(cSystem->vecClusterSize[i]>1&&aryFlagHRDone[i]==false)rotateAbleSize++;

  //Notice: no handle of "rotateAbleSize == 0", it handle by GetNeighbor1
  //Then we choose the rotate target member cluster
  int rotateCursor = (int)((double)rand() / ((double)RAND_MAX + 1) * rotateAbleSize)+1;

  list <list<int> > ::iterator itlist1 = cSystem->listCluMember->begin();
  targetHeadIndex = 0;
  for(int i=0; i<rotateCursor; itlist1++,targetHeadIndex++) //Disconsecutive Candidate
  {
    //if(itlist1->size()>1&&aryFlagHRDone[targetHeadIndex]==false)i++;
    if(itlist1->size()>1)i++;
    if(i==rotateCursor)break;//break befor move to next iterator
  }
  if (itlist1->size()<=1)cout<<"error, search rotate  error"<<endl;
  assert(targetHeadIndex>=0&&targetHeadIndex<maxChNum);



  //Choose targetNode
  targetNode=cSystem->vecHeadName[targetHeadIndex];
  double maxTier1Gain = -DBL_MAX;

  int OriginalNode= cSystem->vecHeadName[targetHeadIndex];
  list<int>::iterator it1=itlist1->begin();
  for(; it1!=itlist1->end(); it1++)//find minimum interference received among nodes in the cluster
  {
    if((*it1)==OriginalNode)continue;
    double tmpTier1Gain = m_ptrMap->GetGi0ByNode(*it1);//==Head Rotate
    if( tmpTier1Gain > maxTier1Gain )
    {
      maxTier1Gain = tmpTier1Gain;
      targetNode=*it1;
    }
  }
  assert(targetHeadIndex != -1 && targetNode != -1);
}

// -------------------------------------------------------------------------- //
// @Description: decideHeadJoining4b
// @Provides: 
// -------------------------------------------------------------------------- //


void 
MinPowerSACluster::decideHeadJoining4b(){
    int threshold=thresholdd; // To determine the cluster size is large enough 
                              // to perform join operation 
    JoiningHeadIndex=-1;
    targetHeadIndex=-1;
    /* maximal rate reduction */
    double maxEntropyDiff = -DBL_MAX;
    int lChIdx = -1;
    int rChIdx = -1;
    for (int i = 0; i < maxChNum; ++i) {
      for (int j = 0; j < maxChNum; ++j) {
        if ( (i != j) &&
              ( cSystem->vecClusterSize.at(i) > 0) &&
              ( cSystem->vecClusterSize.at(j) > 0) &&
              ( cSystem->vecClusterSize.at(i) + cSystem->vecClusterSize.at(j) <= m_tier2NumSlot + 1) 
              ) {
          double lEntropy = GetClusterEntropy(i);
          double rEntropy = GetClusterEntropy(j);
          double jEntropy = GetJoinHeadEntropy(i,j);
          if (lEntropy + rEntropy - jEntropy > maxEntropyDiff) {
            maxEntropyDiff = lEntropy + rEntropy - jEntropy;
            lChIdx = i;
            rChIdx = j;
          }
        }
      }
    }
    JoiningHeadIndex = lChIdx;
    targetHeadIndex = rChIdx; 

}

double
MinPowerSACluster::GetClusterEntropy(const int chIdx)
{
  assert ( chIdx >= 0 && chIdx < maxChNum );
  std::list<std::list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();
  std::vector<int> tmpIndicator(m_ptrMap->GetNumNodes(), 0);
  std::vector<double> tmpVariance(m_ptrMap->GetNumNodes(), 1.0);

  for (int i = 0; i < chIdx; ++i) ++iterRow; 
  std::list<int>::const_iterator iterCol = iterRow->begin();
  for (; iterCol != iterRow->end(); ++iterCol) {
    if (*iterCol > 0 ) {
      tmpIndicator.at(*iterCol) = 1;
    }
  }
  return matrixComputer->GetJointEntropy(tmpIndicator, tmpVariance, 0, m_ptrMap->GetQBits());
}

double
MinPowerSACluster::GetJoinHeadEntropy(const int lChIdx, const int rChIdx)
{
  assert ( lChIdx >= 0 && lChIdx < maxChNum );
  assert ( rChIdx >= 0 && rChIdx < maxChNum );
  std::list<std::list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();
  std::vector<int>    tmpIndicator(m_ptrMap->GetNumNodes(), 0);
  std::vector<double> tmpVariance(m_ptrMap->GetNumNodes(), 1.0);
  /* l ChIdx */
  for (int i = 0; i < lChIdx; ++i) ++iterRow; 
  std::list<int>::const_iterator iterCol = iterRow->begin();
  for (; iterCol != iterRow->end(); ++iterCol) {
    if (*iterCol > 0 ) {
      tmpIndicator.at(*iterCol) = 1;
    }
  }
  /* r ChIdx*/
  iterRow = cSystem->listCluMember->begin();
  for (int i = 0; i < rChIdx; ++i) ++iterRow; 
  iterCol = iterRow->begin();
  for (; iterCol != iterRow->end(); ++iterCol) {
    if (*iterCol > 0 ) {
      tmpIndicator.at(*iterCol) = 1;
    }
  }
  return matrixComputer->GetJointEntropy(tmpIndicator, tmpVariance, 0, m_ptrMap->GetQBits());
}

bool
MinPowerSACluster::CheckTwoLinkFeasible(const int lChName, const int lMemberName, const int rChName, const int rMemberName)
{
  double Gamma = 1.0;
  double snr_require = Gamma * (pow(2, m_ptrMap->GetIdtEntropy()/m_tier1TxTime/m_ptrMap->GetBandwidth()) - 1.0);  
  double lDecisionValue = (snr_require * m_ptrMap->GetNoise() * m_ptrMap->GetGijByPair(lMemberName, lChName) + 
      snr_require * snr_require * m_ptrMap->GetNoise() * m_ptrMap->GetGijByPair(lMemberName,rChName)) /
    (m_ptrMap->GetGijByPair(lMemberName, lChName)* m_ptrMap->GetGijByPair(rMemberName, rChName) - 
     snr_require * snr_require  * m_ptrMap->GetGijByPair(lMemberName, rChName) * m_ptrMap->GetGijByPair(rMemberName, lChName));
  double rDecisionValue = (snr_require * m_ptrMap->GetNoise() * m_ptrMap->GetGijByPair(rMemberName, rChName) + 
      snr_require * snr_require * m_ptrMap->GetNoise() * m_ptrMap->GetGijByPair(rMemberName,lChName)) /
    (m_ptrMap->GetGijByPair(lMemberName, lChName)* m_ptrMap->GetGijByPair(rMemberName, rChName) - 
     snr_require * snr_require * m_ptrMap->GetGijByPair(lMemberName, rChName) * m_ptrMap->GetGijByPair(rMemberName, lChName));
//  cout <<m_ptrMap->GetMaxPower()<<' ' << lDecisionValue << ' ' << rDecisionValue << endl;
  if (lDecisionValue > 0 && lDecisionValue < m_ptrMap->GetMaxPower() && rDecisionValue > 0 && rDecisionValue < m_ptrMap->GetMaxPower()) {
    return true;
  }
  else {
    return false;
  }
}

bool
MinPowerSACluster::CheckLinkFeasible(const int chName, const int name)
{
  std::list<std::list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();
  for (int idx = 0; iterRow != cSystem->listCluMember->end(); ++iterRow, ++idx) {
    std::list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      if ((*iterCol) != name && *iterCol != cSystem->vecHeadName.at(idx) && cSystem->vecHeadName.at(idx) != chName) {
        if(!CheckTwoLinkFeasible(chName, name, cSystem->vecHeadName.at(idx), *iterCol))
            return false;
      }
    }
  }
  return true;

}

bool
MinPowerSACluster::CheckAllFeasible()
{
  std::list<std::list<int> >::const_iterator iterRow =  cSystem->listCluMember->begin();
  for (int idx = 0; iterRow != cSystem->listCluMember->end(); ++iterRow, ++idx) {
    std::list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      if (cSystem->vecHeadName.at(idx) != *iterCol) {
        if (!CheckLinkFeasible(cSystem->vecHeadName.at(idx), *iterCol)) {
          return false;
        }
      }
    }
  }
  return true;

}

double
MinPowerSACluster::GetTier2ExpectPower(const int Name, const int chName)
{
  double denomiator = m_ptrMap->GetNoise();
  double Gamma = 1.0;
  double snr_require = Gamma * (pow(2, m_ptrMap->GetIdtEntropy()/m_tier1TxTime/m_ptrMap->GetBandwidth()) - 1.0);  
  list<list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();  
  for (int chIdx = 0; iterRow !=cSystem->listCluMember->end(); ++iterRow, ++chIdx) {
    list<int>::const_iterator iterCol = iterRow->begin();
    if ( cSystem->vecHeadName.at(chIdx) == -1  || 
        cSystem->vecHeadName.at(chIdx) == chName || 
        cSystem->vecClusterSize.at(chIdx) == 1) continue;
    double maxInterfere = -DBL_MAX;
    int maxIterfereNode = -1;
    for (; iterCol != iterRow->end(); ++iterCol) {
      if ( cSystem->vecHeadName.at(chIdx) != *iterCol &&  m_ptrMap->GetGijByPair(chName, *iterCol) > maxInterfere) {
        maxInterfere = m_ptrMap->GetGijByPair(chName, *iterCol);
        maxIterfereNode = *iterCol;
      }
    }
    denomiator += maxInterfere* m_ptrMap->GetMaxPower();
  }
  return snr_require*denomiator/m_ptrMap->GetGijByPair(Name, chName);
}

bool
MinPowerSACluster::CheckTier2Feasible()
{
  list<list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();  
  for (int chIdx = 0; iterRow != cSystem->listCluMember->end(); ++iterRow, ++chIdx) {
    if (cSystem->vecHeadName.at(chIdx) == -1 ) continue; 
    list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
//        cout << GetTier2ExpectPower(*iterCol, cSystem->vecHeadName.at(chIdx)) <<' ' <<m_ptrMap->GetMaxPower() << ' '
//          <<*iterCol << ' ' << cSystem->vecHeadName.at(chIdx) <<endl; 
      if ( GetTier2ExpectPower(*iterCol, cSystem->vecHeadName.at(chIdx)) > m_ptrMap->GetMaxPower()) {
        return false;
      }
    }
  }
  return true;
}

double
MinPowerSACluster::GetSizePenalty( const vector<double>& sizePenalty)
{
  double tmpSizePenalty = 0.0;
  for (int k = 0; k < cSystem->vecClusterSize.size(); ++k) {
    if (cSystem->vecClusterSize.at(k) != 0
        && (static_cast<double>(cSystem->vecClusterSize.at(k)) - static_cast<double>(m_tier2NumSlot) - 1.0) > 0.0 ) {
      tmpSizePenalty += sizePenalty.at(k) *
        pow(static_cast<double>(cSystem->vecClusterSize.at(k)) - static_cast<double>(m_tier2NumSlot) - 1.0,5); 
    }
  }
  return tmpSizePenalty;

}
double
MinPowerSACluster::GetTier2Penalty( const vector<double>& tier2Penalty )
{
  double penalty = 0.0;
  list<list<int> >::const_iterator iterRow = cSystem->listCluMember->begin();  
  for (int chIdx = 0; iterRow != cSystem->listCluMember->end(); ++iterRow, ++chIdx) {
    if (cSystem->vecHeadName.at(chIdx) == -1 ) continue; 
    list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      if ( GetTier2ExpectPower(*iterCol, cSystem->vecHeadName.at(chIdx)) > m_ptrMap->GetMaxPower()) {
        penalty += tier2Penalty.at(*iterCol) * 
          GetTier2ExpectPower(*iterCol, cSystem->vecHeadName.at(chIdx)) - m_ptrMap->GetMaxPower();
      }
    }
  }
  return penalty;
}

double
MinPowerSACluster::GetEntropyPenalty(const double entropyPenalty)
{
  double tmpJEntropy = nextSupNum*indEntropy + matrixComputer->computeLog2Det(1.0, cSystem->allSupStru);
  if (fidelityRatio*wholeSystemEntopy - tmpJEntropy < 0) {
    return 0.0;
  }
  else {
    return entropyPenalty * (fidelityRatio*wholeSystemEntopy - tmpJEntropy);
  }
}

// -------------------------------------------------------------------------- //
// @Description: targetHeadIndex, IsolateNodeName
// @Provides: 
// -------------------------------------------------------------------------- //

void MinPowerSACluster::decideIsolate4b(){

    /* minimal rate increasing */
    list<list<int> >::iterator itli1 = cSystem->listCluMember->begin();
    int discardableSize=0;
    for(; itli1!=cSystem->listCluMember->end(); itli1++)if(itli1->size()>1)discardableSize++;
    assert(discardableSize!=0);
    //Randomly choose a cluster
    itli1 = cSystem->listCluMember->begin();
    isolatedHeadIndex = -1;
    int chooseCursor = (int) ((double)rand() / ((double)RAND_MAX + 1) * discardableSize)+1;//+1 Becasue the intial might
    assert(chooseCursor<=maxChNum&&chooseCursor>=0);

    itli1 = cSystem->listCluMember->begin();
    isolatedHeadIndex=0;
    for(int i=0; i<chooseCursor; itli1++,isolatedHeadIndex++)
    {
        if(itli1->size()>1)i++;
        if(i==chooseCursor)break;
    }

    assert(cSystem->vecHeadName[isolatedHeadIndex]!=-1);
    assert(cSystem->vecClusterSize[isolatedHeadIndex]==itli1->size());
    assert(cSystem->vecClusterSize[isolatedHeadIndex]>1);
    double minGain = DBL_MAX;
    double tempInter = 0;
    int    minGainName = -1;
    list<int>::iterator it2=itli1->begin();
    for(; it2!=itli1->end(); it2++) {
        tempInter = 0;
        if(cSystem->vecHeadName[isolatedHeadIndex]==(*it2))continue;
        else {
            if( m_ptrMap->GetGijByPair(cSystem->vecHeadName[isolatedHeadIndex], (*it2) ) < minGain) {
                minGain = m_ptrMap->GetGijByPair(cSystem->vecHeadName[isolatedHeadIndex], (*it2) ) ;
                minGainName = (*it2);
            }
        }
    }
    assert(minGainName!=-1);
    for(int i=0;i<cSystem->vecHeadName.size();i++){
        //cout<<i<<" "<<cSystem->vecHeadName[i]<<endl;
        if(cSystem->vecHeadName[i]==-1){
            targetHeadIndex=i;
            break;
        }
    }

    IsolateNodeName = minGainName;
    assert(isolatedHeadIndex!=-1 && IsolateNodeName!=-1);
}

/*
   Movement Function
   add new member and adjust the "nodes" value

*/
void MinPowerSACluster::addMemberSA(int inputHeadIndex, int inputMemberName)
{
  cSystem->addMemberCs(inputHeadIndex,inputMemberName,iniDone);
  nodes[inputMemberName].ptrHead = cSystem->returnHeadPtr(inputHeadIndex);
}
/*
   Movement Function
   discard the node form the specific head(cluster)
   */
void MinPowerSACluster::discardMemberSA(int inputHeadIndex, int inputMemberName)
{
  cSystem->discardMemberCs(inputHeadIndex,inputMemberName);
  ptrHeadLastDiscard = nodes[inputMemberName].ptrHead;
  powerLastDiscard = nodes[inputMemberName].power;
  nodes[inputMemberName].ptrHead=NULL;
  nodes[inputMemberName].power =0;
}


/*
   Rotate Member SA
   */
void MinPowerSACluster::rotateHeadSA(int inputHeadIndex, int inputMemberName)
{
	rotatedHeadNameLast = cSystem->vecHeadName[inputHeadIndex];
	cSystem->vecHeadName[inputHeadIndex] = inputMemberName;
}

void MinPowerSACluster::isolateHeadSA(int isoName,int IsolateCluI, int targetH){
    //cout<<"doing isolation: iso-"<<isoName<<" From "<<IsolateCluI<<" to "<<targetH<<endl;


    lastIsolateNodeName=isoName;
    lastIsolatedClusterIndex=IsolateCluI;
    //slot= targetHeadIndex
    nodes[isoName].ptrHead=&(cSystem->vecHeadName[targetHeadIndex]);
    cSystem->isolate_FromHead(isoName,IsolateCluI,targetH);
    nextChNum=curChNum+1;
    //cout<<"nexyCh"<<nextChNum<<endl;

}


void MinPowerSACluster::join_fromHeadSA(int JoiningHeadIndex,int targetH){

     lastJoingingMachine.clear();
     lastJoiningHeadIndex=JoiningHeadIndex;
     lastJoiningHead=cSystem->vecHeadName[lastJoiningHeadIndex];
     std::list<std::list <int> >::iterator it_LiInt=cSystem->listCluMember->begin();
     for(int i=0;i<JoiningHeadIndex;i++)it_LiInt++;
     std::list<int>::iterator it_Int=it_LiInt->begin();
     for(int i=0;i<it_LiInt->size();i++,it_Int++){
         lastJoingingMachine.push_back(*it_Int);
         nodes[*it_Int].ptrHead=&(cSystem->vecHeadName[targetH]);
     }

     cSystem->join_FromHead(JoiningHeadIndex,targetH);
     nextChNum=curChNum-1;

}


void MinPowerSACluster::calculateMatrics_minResors(//Calculate next performance matircs
    const vector<double>& sizePenalty, 
    const vector<double>& tier2Penalty, 
    const double& entropyPenalty)
{
  next2nd_Joule = returnTransientJoule();
  next1st_Joule= power1st*next1st_ms/1000;
  //:<<"MaxPowerLoad="<<returnMaxPowerLoad()<<" Max="<<powerMax<<endl;

  if (nextEventFlag==1)
  {
    nextSupNum = curSupNum + 1;//Serve one more node
    nextJEntropy = nextSupNum*indEntropy + matrixComputer->computeLog2Det(1.0, cSystem->allSupStru);
  }
  else if(nextEventFlag==2)
  {
    nextSupNum = curSupNum -1;
    nextJEntropy = nextSupNum * indEntropy + matrixComputer->computeLog2Det(1.0, cSystem->allSupStru);
  }
  else if (nextEventFlag==3)  //Do Nothing
  {
    //cout<<"Original Head="<<rotatedHeadNameLast<<endl;

  }
  else if (nextEventFlag==4)  //Do Nothing
  {
    //cout<<"Original Head="<<rotatedHeadNameLast<<endl;

  }
  else if (nextEventFlag==5)  //Do Nothing
  {
    //cout<<"Original Head="<<rotatedHeadNameLast<<endl;

  }
  else
  {
    cout<<"Error in calculateMatrics_minResors "<<endl;
    assert(0);
  }
  nextPayoff = GetPayOff(sizePenalty, tier2Penalty, entropyPenalty);
}

/*
   after do the neighbor change
   -confirm the neighbor change and decide whether reverse or not
   -confirm3c add reset some metric if structure change(add,discard)
   by"if(nextEventFlag==1||nextEventFlag==2)confirmStructureChange();"
   */
void MinPowerSACluster::ConfirmNeighbor1()
{
  //constraint: the all the machine must be supported
  bool nextAllServe = (nextJEntropy>(fidelityRatio*wholeSystemEntopy)?true:false) && CheckTier2Feasible();
  for (int i = 0; i < cSystem->vecHeadName.size(); ++i) {
    if (cSystem->vecHeadName.at(i) >= 0 && cSystem->vecClusterSize.at(i) > m_tier2NumSlot + 1) {
      nextAllServe = false;
    }
  }
  bool curAllServe = (curJEntropy>(fidelityRatio*wholeSystemEntopy)?true:false) && CheckTier2Feasible() ; 

  for (int i = 0; i < m_prevVecClusterSize.size(); ++i) {
    if (m_prevVecHeadName.at(i) >= 0 && m_prevVecClusterSize.at(i) > m_tier2NumSlot + 1) {
      curAllServe = false;
    }
  }
  nextAllServe = true;
  curAllServe = true;
  /* decision flow */
  if ( !curAllServe && nextAllServe ) {
    passNext2Cur();
    for(int i=0; i<totalNodes; i++) 
      nodes[i].power = nextNodePower[i];
    if( nextEventFlag == 1 || nextEventFlag == 2 ) 
      confirmStructureChange();
  }
  else if ( ( nextPayoff < m_curPayoff )  )
  {
    passNext2Cur();
    for(int i=0; i<totalNodes; i++)
      nodes[i].power = nextNodePower[i];
    if( nextEventFlag == 1 || nextEventFlag == 2 )
      confirmStructureChange();
  }
  else if( 
       ( nextPayoff > m_curPayoff )  )
  {
    double probAnnealing = exp (-1.0*abs(nextPayoff-m_curPayoff)/temparature);
    //cout<<"Show Payoff "<<nextPayoff<<"  "<<m_curPayoff<<endl;
    //cout<<"  Prob Annealing:  "<<probAnnealing<<endl;
    double annealingChoose = (double)rand()/((double)RAND_MAX+1);
    if ( annealingChoose > probAnnealing )
    {
      reverseMoveSA();
    }
    else//accept the move
    {
      passNext2Cur();
      for(int i=0; i<totalNodes; i++)
        nodes[i].power = nextNodePower[i];
      if(nextEventFlag==1||nextEventFlag==2)
        confirmStructureChange();
    }
  }
  targetHeadIndex = -1;
  targetNode = -1;
  nextEventFlag = -1;
}

void
MinPowerSACluster::ConfirmNeighbor2(
    vector<double>& sizePenalty, 
    vector<double>& tier2Penalty, 
    double& entropyPenalty,
    const vector<double>& tmpSizePenalty,
    const vector<double>& tmpTier2Penalty,
    const double& tmpEntropyPenalty)
{
  double tmpPayoff = GetPayOff(tmpSizePenalty, tmpTier2Penalty, tmpEntropyPenalty);

  if (tmpPayoff < m_curPayoff) {
    double probAnnealing = exp (-1.0*abs(tmpPayoff-m_curPayoff)/temparature);
    double annealingChoose = (double)rand()/((double)RAND_MAX+1);
    if ( annealingChoose < probAnnealing ) {//accept the move
      sizePenalty.assign(tmpSizePenalty.begin(), tmpSizePenalty.end());
      tier2Penalty.assign(tmpTier2Penalty.begin(), tmpTier2Penalty.end());
      entropyPenalty = tmpEntropyPenalty;
      m_curPayoff = tmpPayoff;
    }
  }
  else {
      sizePenalty.assign(tmpSizePenalty.begin(), tmpSizePenalty.end());
      tier2Penalty.assign(tmpTier2Penalty.begin(), tmpTier2Penalty.end());
      entropyPenalty = tmpEntropyPenalty;
  }
}

void 
MinPowerSACluster::passNext2Cur() 
{


    curJEntropy = nextJEntropy;
    curSupNum = nextSupNum;

    curChNum=nextChNum;

    m_curPayoff = nextPayoff;

    cur1st_Joule = next1st_Joule;
    cur1st_ms = next1st_ms;

    cur2nd_Joule = next2nd_Joule;
    cur2nd_ms=next2nd_ms;

    switch (nextEventFlag) {
    case 1:
        //cout<<"add "<<targetNode<<" to "<<cSystem->vecHeadName[targetHeadIndex]<<endl;
        tempAddT++;
        break;

    case 2:
        //cout<<"discard "<<targetNode+1<<" from "<<cSystem->vecHeadName[targetHeadIndex]+1<<endl;
        tempDisT++;
        break;
    case 3:
        //cout<<"HR "<<targetNode+1<<" to Replace "<<rotatedHeadNameLast+1<<endl;
        tempHRT++;
        break;
    case 4:
        lastJoinPassAccu=0;
        tempJoinT++;
        break;
    case 5:
        //cout<<"Pass Iso"<<endl;
        lastIsoPassAccu=0;
        tempIsoT++;
        //cout<<"nextChNum"<<nextChNum;
        break;
    }


}



/*
    reverse the last move in Csystem and node
*/
void MinPowerSACluster::reverseMoveSA()
{
    //cout<<"reverse"<<endl;
    if (nextEventFlag == 1)//reverse the add previously did
    {
        cSystem->reverseAdd(targetHeadIndex, targetNode);
        //delete in the ULAGENT
        nodes[targetNode].ptrHead = NULL;
        nodes[targetNode].power = 0;
        nextEventFlag = -1;
    }
    else if (nextEventFlag ==2)//reverse the discard previously did
    {
        cSystem->reverseDiscard(targetHeadIndex, targetNode);
        if(ptrHeadLastDiscard!=NULL)
        {
            nodes[targetNode].ptrHead = ptrHeadLastDiscard;
            nodes[targetNode].power = powerLastDiscard;
        }
        else
        {
            cout<<"Error, the ptrLastDiscard should't be NULL"<<endl;
            assert(1);
        }
        ptrHeadLastDiscard = NULL;
        powerLastDiscard = 0;
        nextEventFlag =-1;
    }
    else if (nextEventFlag == 3)//reverse the HeadRotate previously did
    {
        cSystem->vecHeadName[targetHeadIndex]=rotatedHeadNameLast;
    }
    else if(nextEventFlag == 4){

        cSystem->reverseJoin(lastJoiningHeadIndex,lastJoiningHead,targetHeadIndex,lastJoingingMachine);
        for(unsigned int i;i<lastJoingingMachine.size();i++){
            nodes[lastJoingingMachine[i]].ptrHead=&cSystem->vecHeadName[lastJoiningHeadIndex];
        }
        lastJoingingMachine.clear();
        --nextChNum = curChNum ;

    }
    else if(nextEventFlag == 5){
         cSystem->reverseisolate(lastIsolateNodeName,lastIsolatedClusterIndex,targetHeadIndex);
        nodes[lastIsolateNodeName].ptrHead=&cSystem->vecHeadName[lastIsolatedClusterIndex];
        --nextChNum ;
    }
    else cout<<"Error next flag ="<<nextEventFlag<<"wrong"<<endl;
}

/*
  Metrics update after structure changed
  -reset "aryFlagHRDone"
*/
void MinPowerSACluster::confirmStructureChange()
{
  //cout<<"Structure Change"<<endl;
  for(int i=0; i<maxChNum; i++)aryFlagHRDone[i]=false;
  for (int i = 0; i < totalNodes; ++i) {
    prevAllSupStru[i] = cSystem->allSupStru[i];
  }
  m_prevVecClusterSize.assign(cSystem->vecClusterSize.begin(), cSystem->vecClusterSize.end());
  m_prevVecHeadName.assign(cSystem->vecHeadName.begin(), cSystem->vecHeadName.end());
}



/*
    Check if this structure is the best
*/
bool MinPowerSACluster::checkBestClusterStructure_DataCentric(int inputRound)
{
//    bool curAllServe = (curJEntropy>fidelityRatio*wholeSystemEntopy?true:false);
  /* new constraint */
  bool curAllServe = ((curJEntropy >= fidelityRatio*wholeSystemEntopy?true:false) && CheckTier2Feasible()); 
  bool sizeFeasible = true;
  for (int i = 0; i < m_prevVecClusterSize.size(); ++i) {
    if ( m_prevVecHeadName.at(i) && m_prevVecClusterSize.at(i) > m_tier2NumSlot + 1) {
      curAllServe = false;
      sizeFeasible = false;
    }
  }
//#ifdef DEBUG
  cout << "IterSA: " << inputRound << endl;
  for (int i = 0; i < cSystem->vecClusterSize.size(); ++i) {
    cout << cSystem->vecClusterSize.at(i) << ' ';
  }
  cout << endl;
  cout << "E L S: " << (curJEntropy >= fidelityRatio*wholeSystemEntopy) <<' '<<CheckTier2Feasible() <<' ' <<sizeFeasible << ", ";
  cout << "min payoff: " << bestFeasiblePayoff << ", ";
  cout << "curr Payoff: " << m_curPayoff << endl;
//#endif
  if(curAllServe)bestAllServeFound=true;
  //cout<<"In check Best"<<endl;
  if ((curJEntropy>bestFeasibleJEntropy) && !curAllServe && !bestAllServeFound)
  {
    //cout<<"In check Best"<<endl;
    roundBest = inputRound;
    bestFeasibleJEntropy=curJEntropy;
    bestFeasibleSupNum=curSupNum ;
    //bestFeasiblePayoff=m_curPayoff;
    best1st_Joule=cur1st_Joule;
    best1st_ms=cur1st_ms;
    best2nd_Joule=cur2nd_Joule;
    best2nd_ms=cur2nd_ms;
    keepBestStructure();
    bestChNum=curChNum;

  }
  else if ((m_curPayoff<bestFeasiblePayoff) && curAllServe)
  {
    ///cout<<"Find new best"<<endl;
    roundBest = inputRound;
    bestFeasibleJEntropy=curJEntropy;
    bestFeasibleSupNum=curSupNum ;
    bestFeasiblePayoff=m_curPayoff;
    best1st_Joule=cur1st_Joule;
    best1st_ms=cur1st_ms;
    best2nd_Joule=cur2nd_Joule;
    best2nd_ms=cur2nd_ms;
    bestChNum=curChNum;
    keepBestStructure();
  }
  //cout<<"best FeasiblePayoff="<<bestFeasiblePayoff<<" with headNum="<<bestChNum<<";Info Ratio="<<bestFeasibleJEntropy/wholeSystemEntopy<<endl;
  return false;
}


void MinPowerSACluster::keepBestStructure()
{
    vecHeadNameBest.assign(cSystem->vecHeadName.begin(),cSystem->vecHeadName.end());
    listCluMemBest.assign(cSystem->listCluMember->begin(), cSystem->listCluMember->end());
    //1st tier computation has been done in the calculateMatrics_minResors()
    vecBestClusterBits.assign(vecClusterHeadBits.begin(),vecClusterHeadBits.end());
    vecBestClusterSize.assign(cSystem->vecClusterSize.begin(),cSystem->vecClusterSize.end());
    vecBestClusterHeadMS.assign(vecClusterHeadMS.begin(),vecClusterHeadMS.end());
    vecBestClusterHeadWatt.assign(vecClusterHeadWatt.begin(),vecClusterHeadWatt.end());
    for(int i=0; i<maxChNum; i++)
    {
        for(int j=0; j<totalNodes; j++)
        {
            bestMaClusterStru[i][j] = cSystem->clusterStru[i][j];
        }
    }
    for(int i =0; i<totalNodes; i++)
    {
        powerBest[i] = nextNodePower[i];
        // cout<<"InKEEP "<<i<<" Power"<< powerBest[i]<<endl;
    }
    for(int i=0; i<totalNodes; i++)
        bestAllSupStru[i]=cSystem->allSupStru[i];
}
/*
    Internal tool
    Purpose: return the closet node index from a certain (X,Y)

*/
int MinPowerSACluster::returnClosetNodeIndexInGroup(int tempX,int tempY, std::vector<int> &inputGroup)
{
    float closetNodeDistance =  numeric_limits<float>::max( );
    float tempD = 0;
    int closetNodeIndex = -1;
    for(unsigned int i =0; i
            < inputGroup.size() ; i++)
    {
        tempD = (nodes[inputGroup[i]].locX - tempX )*(nodes[inputGroup[i]].locX - tempX)\
                +(nodes[inputGroup[i]].locY - tempY) * (nodes[inputGroup[i]].locY - tempY);
        if (tempD<closetNodeDistance)
        {
            closetNodeDistance = tempD;
            closetNodeIndex =i;
        }
    }
    return closetNodeIndex;
}

double MinPowerSACluster::returnTransientJoule() {
    std::list <std::list<int> >::iterator itlist1=cSystem->listCluMember->begin();
    double accuJoule=0;
    for(int i =0; itlist1!=cSystem->listCluMember->end(); itlist1++,i++)
    {
        std::list<int>::iterator it1=itlist1->begin();
        double tempSize = static_cast<double> (itlist1->size());
        if (tempSize==1)continue;
        for(; it1!=itlist1->end(); it1++) {
            if (*it1 ==cSystem->vecHeadName[i])continue;
            accuJoule+=nextNodePower[(*it1)]*cur2nd_ms/(tempSize-1);
            //  EYESTEVEN
            /*  cout<<"node "<<setw(4)<<*it1<<" : "<<setw(12)<<nextNodePower[(*it1)]<<"(Watt), "<<setw(12)<<nextNodePower[(*it1)]*cur2nd_ms/ \
            //  (tempSize-1)/1000<<"(Joule) " <<" to Head "<<setw(4)<<cSystem->vecHeadName[i]<<" with "<<setw(4)<<(int)tempSize << " in same cluster" <<endl;
            */
        }
    }
    accuJoule/=1000;
    return (accuJoule);
}

void MinPowerSACluster::resetSA3iSystem() {
    iniDone = false;
    bestFeasibleJEntropy = -1;
    bestFeasibleSupNum = -1;
    bestFeasiblePayoff=DBL_MAX;

    temparature = constantIniTemprature;

    //next series dpn't need to be reseted;
    for(int i=0; i<totalNodes; i++)nodes[i].ptrHead=NULL;
    for(int i=0; i<maxChNum; i++)aryFlagHRDone[i]=false;
    vecClusterHeadBits.clear();
    vecClusterHeadBits.resize(maxChNum);
    vecClusterHeadMS.clear();
    vecClusterHeadMS.resize(maxChNum);
    vecClusterHeadWatt.clear();
    vecClusterHeadWatt.resize(maxChNum);

    vecBestClusterBits.clear();
    vecBestClusterBits.resize(maxChNum);
    vecBestClusterHeadMS.clear();
    vecBestClusterHeadMS.resize(maxChNum);
    vecBestClusterHeadWatt.clear();
    vecBestClusterHeadWatt.resize(maxChNum);
    vecBestClusterSize.clear();
    vecBestClusterSize.resize(maxChNum);
}

