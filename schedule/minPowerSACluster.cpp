#include <iostream>
#include<iomanip>
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

using namespace std;
#include "minPowerSACluster.h"
MinPowerSACluster::MinPowerSACluster() {};

MinPowerSACluster::MinPowerSACluster(FILE *fileReadCursor, int inputTotalNodes, int inputMaxChNum,int inputSAFac,  \
                     int inOutputControl,
                     int isStrucOuput,
                     double inputTemprature, double InputSaAlpha, \
                     double inCorrelationFactor, string ipAddr, Map const * const myPtrMap)
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
    while(inputIndex<totalNodes)
    {
        if( fscanf(fileReadCursor, "%f", &dx)<0 || fscanf(fileReadCursor, "%f", &dy)<0 )
        {
            cout << "File map's content is wrong!!\n";
            assert(0);
        }
        else
        {
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
    for(int i=0; i<totalNodes; i++)
    {
        distanceOf2Nodes[i] = new double [totalNodes];
        Gij[i] = new float [totalNodes];
        Gib[i] = sysComputing->returnChannelGain_BS(nodes[i]);
        for(int j = 0; j<totalNodes; j++)
        {
            if (i==j)
            {
                distanceOf2Nodes[i][j]=0.0;
                Gij[i][j] = 1.0;
            }
            else if (i>j)
            {
                distanceOf2Nodes[i][j]=distanceOf2Nodes[j][i];
                Gij[i][j] = Gij[j][i];
            }
            else
            {
                Gij[i][j] = sysComputing->returnChannelGain_2Nodes(nodes[i],nodes[j]);
                double tempDib = pow(nodes[i].locX-nodes[j].locX,2) +pow(nodes[i].locY - nodes[j].locY,2);
                distanceOf2Nodes[i][j] = tempDib;
            }
        }
    }

    /* for (int i=0;i<totalNodes;i++){
       for(int j=0;j<totalNodes;j++)
         cout<<setw(14)<<Gij[i][j]<<" ";
       cout<<endl;
     }*/
    cSystem = new ULCS1b(inputTotalNodes, inputMaxChNum);
    powerBest = new double [totalNodes];
    nextNodePower = new double [totalNodes];
    for (int i=0; i<maxChNum; i++ )
        for(int j=0; j<totalNodes; j++)
            cSystem->clusterStru[i][j] = false;
    maIndexSortDecGain = new int* [totalNodes];
    for (int i =0; i<totalNodes; i++)
        maIndexSortDecGain[i]= new int [totalNodes];
    sortIndex = -1;//Just a check point
    //Keep the sturcture for best Structure
    aryFlagHRDone = new bool [maxChNum];
    for(int i=0; i<maxChNum; i++)aryFlagHRDone[i]=false;

    //-----Best Performance Index-----//
    listCluMemBest = new list<list <int> >;

    //Interference Matrix Keep the interference index/Name set from othe cluster for each cluster.
    bestMaClusterStru = new bool *[maxChNum];
    for (int i=0; i<maxChNum; i++)
    {
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

    mapNodeName2DistanceRank = new map<int, int> [totalNodes];

    for(int i=0; i<totalNodes; i++)
    {
        for(int j=0; j<totalNodes; j++) maIndexSortDecGain[i][j] = j;
        sortIndex = i;
        sort(maIndexSortDecGain[i], maIndexSortDecGain[i]+totalNodes, compareDis(*this));
    }

    for(int i=0; i<totalNodes; i++)
    {
        for(int j=0; j<totalNodes; j++)
        {
            mapNodeName2DistanceRank[i][maIndexSortDecGain[i][j]]=j;
        }
    }

}
MinPowerSACluster::~MinPowerSACluster()
{
    if(!terminated)releaseMemory();
}

vector<int>
MinPowerSACluster::GetAllSupStru() const {
  vector<int> myAllSupStru(totalNodes,0);
  list<list<int> >::const_iterator iterRow = listCluMemBest->begin();
  for (;iterRow != listCluMemBest->end(); ++iterRow) {
    list<int>::const_iterator iterCol = iterRow->begin();
    for (; iterCol != iterRow->end(); ++iterCol) {
      myAllSupStru.at(*iterCol) = 1;
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
    //release maIndexSortDecGain
    for (int i = 0; i<totalNodes; i++) delete [] maIndexSortDecGain[i];
    delete maIndexSortDecGain;
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
    delete listCluMemBest;
    delete []  mapNodeName2DistanceRank;
    delete [] aryFlagHRDone;
    delete [] sysComputing;
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
    and set the initial interference list for each node
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
    else if (!strcmp(iniFlag, "HeadLimited")) 
      normalFlag = setIniHeadLimited();


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
  vector <vector <int> > tempGroup;
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
      vector <int> tempV;
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


bool MinPowerSACluster::setIniHeadLimited()
{
    double sortGib[totalNodes];
    for(int i=0; i<totalNodes; i++) {
        sortGib[i]=Gib[i];
    }
    sort(sortGib,sortGib+totalNodes);

    headCandidatesNum=maxChNum*3;
    for (int i=0; i<totalNodes; i++)
    {
        double targetGainL=(i==totalNodes-1)?0:sortGib[totalNodes-i-2];
        double targetGainR = sortGib[totalNodes-i-1];
        for(int j=0; j<totalNodes; j++) {
            if(Gib[j]<=targetGainR&&Gib[j]>targetGainL&&i<maxChNum) {
                vecHeadCandidates.push_back(j);
                cSystem->addNewHeadCs(j);
                addMemberSAIni(i, j);
                break;
            }
            else if (Gib[j]<=targetGainR&&Gib[j]>targetGainL&&i<headCandidatesNum) {

                vecHeadCandidates.push_back(j);
                cSystem->listUnSupport->push_back(j);
                break;
            }
            else if (Gib[j]<=targetGainR&&Gib[j]>targetGainL) {

                cSystem->listUnSupport->push_back(j);
                break;
            }
        }
    }

    for (int i=0; i<maxChNum; i++)
    {
        nodes[cSystem->vecHeadName[i]].ptrHead = &(cSystem->vecHeadName[i]);
    }
    return true;
}

bool MinPowerSACluster::setIniStruDistanceKmedoids() 
{
  int retryTimes = 0;
  float tempHeadX [maxChNum];
  float tempHeadY [maxChNum];
  int tempHeadList [maxChNum];
  vector <vector <int> > tempGroup;
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
      vector <int> tempV;
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
#ifdef DEBUG
      //find the k-means coordinate of each cluster
      for (int i = 0; i < tempGroup.size(); i++) {
        cout << "cluster: " << i <<"-th ";
        for (int j = 0; j < tempGroup[i].size(); j++) {
          cout << tempGroup[i][j] << ' ';
        }
        cout << endl;
      }
#endif
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




void MinPowerSACluster::writeStruSingleRound(int round)
{
    FILE *fid;
    char str[50];
    cout<<" writeStruSingleRound "<<endl;
    sprintf(str,"LocalWork/loopStructure/roundStru%d.txt",round);
    cout<<str<<endl;
    fid = fopen(str,"w");
    fprintf(fid,"%d\n",totalNodes);
    fprintf(fid,"%d\n",maxChNum);
    fprintf(fid,"%e\n",powerMax);
    fprintf(fid,"%e\n",C2);
    fprintf(fid,"%d\n", bestFeasibleSupNum-maxChNum);
    fprintf(fid,"%d\n", SAIter);
    fprintf(fid,"%f\n", computingTimes);
    for (int i= 0; i<maxChNum; i++)
    {
        for( int j=0; j<totalNodes; j++)
        {
            fprintf(fid,"%d ",cSystem->clusterStru[i][j]);
        }
        fprintf(fid,"\n");
    }

    for (int i= 0; i<maxChNum; i++)
    {
        fprintf(fid,"%d ",cSystem->vecHeadName[i]+1);
    }
    fprintf(fid,"\n");
    for (int i= 0; i<totalNodes; i++)
    {
        fprintf(fid,"%E ", nodes[i].power);
    }
    fprintf(fid,"\n");
    list <list <int> >::iterator it1=cSystem->listCluMember->begin();
    for (; it1!=cSystem->listCluMember->end(); it1++ )
    {
        list <int>::iterator it2 = it1->begin();
        for (; it2!=it1->end(); it2++)
            fprintf(fid,"%d ",(*it2)+1);
        fprintf(fid,"\n");
    }
    fclose(fid);
}


void MinPowerSACluster::writePayoffEachRound_MinResors(int round)
{
    char str[]="MinPowerSAClusterEachR.txt";
    FILE* fid = fopen(str,"a");
    fprintf(fid,"%f %f %f %f %f %d\n", curJEntropy/wholeSystemEntopy,static_cast<double>(curSupNum)/static_cast<double>(totalNodes), curPayoff, cur1st_ms, cur2nd_ms,  round);
    fclose(fid);

}

void MinPowerSACluster::writePayoffEachRound_MinResors_withHead(int round,int head)
{
    char str[]="MinPowerSAClusterEachR.txt";
    FILE* fid = fopen(str,"a");
    fprintf(fid,"%f %f %f %f %f %d %d\n", curJEntropy/wholeSystemEntopy,static_cast<double>(curSupNum)/static_cast<double>(totalNodes), curPayoff, cur1st_ms, cur2nd_ms,  round, head);
    fclose(fid);

}


void MinPowerSACluster::debug_CheckSizeCorrect() {
    for (unsigned int i=0; i<cSystem->vecHeadName.size(); i++) {
        int tempSize = 0;
        for(int j=0; j<totalNodes; j++) {
            if(cSystem->clusterStru[i][j]==1)tempSize++;
        }
        assert(tempSize == cSystem->vecClusterSize[i]);
    }
}



void MinPowerSACluster::computeBestTRR_DataCentric()
{
    best2ndTierTraffic = 0;
    best1stTierTraffic = 0;
    int bestGroupSize[maxChNum];
    for (int i=0; i<maxChNum; i++)bestGroupSize[i]=0;
    double groupRedundancy=0;

    for(int i=0; i<maxChNum; i++)
    {
        for (int j=0; j<totalNodes; j++)
        {
            if (bestMaClusterStru[i][j]==true)
                bestGroupSize[i]++;
        }

        groupRedundancy = matrixComputer->computeLog2Det(1.0,bestMaClusterStru[i]);
        best2ndTierTraffic = best2ndTierTraffic+ bestGroupSize[i]*indEntropy;
        best1stTierTraffic = best1stTierTraffic + bestGroupSize[i]*indEntropy + groupRedundancy;
    }
    bestTrafficReductionRatio = best1stTierTraffic/best2ndTierTraffic;
}


void MinPowerSACluster::computeUpperResourceNoCodingNoPowerControl()
{
    bestUpperLayerResource=0;
    double tempDataLoad[maxChNum];
    for (int i=0; i<maxChNum; i++)
    {
        int vecClusterSizeLocal=0;
        for(int j=0; j<totalNodes; j++)
        {
            if(bestMaClusterStru[i][j]==1)
            {
                vecClusterSizeLocal++;
            }
        }
        tempDataLoad[i]=vecClusterSizeLocal*dataBits;
        double tempResor=tempDataLoad[i]*virtualCompression/rateibMax[vecHeadNameBest[i]]*1000;//ms
        cout<<i<<"-th Cluster  Head "<<vecHeadNameBest[i]<<" RateMax " <<rateibMax[vecHeadNameBest[i]]<<" Data amount "<<dataBits<<" UpperResource"<<tempResor<<endl;
        bestUpperLayerResource+=tempResor;
    }
    char str3[500];
    sprintf(str3,"HumanDataULSA2iNoCoding%d.txt",maxChNum);
    FILE *fid=fopen(str3,"a");
    fprintf(fid,"------------------\n");
    double totalResors=0;
    fprintf(fid,"SupNum %d\n" , bestFeasibleSupNum);
    for(int i=0; i<maxChNum; i++)
    {
        fprintf(fid,"CR:%f %d-th Head %d, Rate(bps):%e, DataLoad(bits):%f, 2nd tier Resource:%f, 1st tier Resource:%f \n", \
                1-virtualCompression,i ,vecHeadNameBest[i] ,rateibMax[vecHeadNameBest[i]], tempDataLoad[i], cur2nd_ms, tempDataLoad[i]*1000/rateibMax[vecHeadNameBest[i]]);
        totalResors+=(tempDataLoad[i]*1000/rateibMax[vecHeadNameBest[i]]);
    }
    totalResors+=cur2nd_ms;
    fprintf(fid,"totalResors %f\n",totalResors);
    fclose(fid);
}

void MinPowerSACluster::computeBestAvgInterference()
{
    list<list <int> >::iterator itl = listCluMemBest->begin();
    for(int i=0; i<maxChNum; i++,itl++) //For each CH we will see there member 1-by-1
    {
        list<list <int> >:: iterator itl2 = listCluMemBest->begin(); //search other cluster from start for each CH
        for (int j=0; j<maxChNum; j++,itl2++) //Search other cluster for each member in cluster (*vecHeadName)[i]
        {
            int interferentSource = -1;
            float maxInterference = -99;
            if (i==j)//Same headIndex which means this is self head.
            {
                maBestInterfernceIndex[i][j]=-1;
                maBestInterference[i][j]=0;
                continue;

            }
            else if (itl2->size()==1)
            {
                maBestInterfernceIndex[i][j] = -1;
                maBestInterference[i][j] = 0;
                continue;
            }
            //Find the strongest interference from other group
            list<int> ::iterator itl3 = itl2->begin();
            for(; itl3!=itl2->end(); itl3++) //search every node(itl3) under each cluster(itl2) != self cluster
            {
                float tempInterference = powerBest[(*itl3)] * Gij[vecHeadNameBest[i]][(*itl3)];
                if(vecHeadNameBest[j]==(*itl3))continue;//other head won't transmit power
                else if (tempInterference>maxInterference )
                {
                    interferentSource = (*itl3);
                    maxInterference = tempInterference;
                }
            }
            maBestInterfernceIndex[i][j] = interferentSource;
            maBestInterference[i][j] = maxInterference;//Interference form j(HeadIndex) to i(HeadIndex)
        }
    }

    bestAvgInterference=0;
    for (int i=0; i<maxChNum; i++)
        for(int j=0; j<maxChNum; j++)
            bestAvgInterference+=maBestInterference[i][j];
    bestAvgInterference/=maxChNum;

}

double MinPowerSACluster::returnComprRatio()
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
  matrixComputer = new CORRE_MA_OPE(totalNodes, correlationFactor, distanceOf2Nodes, (double)quantizationBits);
  indEntropy = 0.5*log2(2*3.1415*exp(1))+quantizationBits;
  double tmpCompR = matrixComputer->returnNSetCorrelationFactorByCompressionRatio \
                    (compRatio,indEntropy,static_cast<double>(totalNodes));
  tempAddT=0;
  tempDisT=0;
  tempHRT=0;
  bool inClu[totalNodes];
  for(int i=0; i<totalNodes; i++)inClu[i]=true;
  double sysRedundancy =matrixComputer->computeLog2Det(1.0, inClu);
  wholeSystemEntopy = totalNodes*indEntropy+sysRedundancy;

  ULSAOutputToolSet<class MinPowerSACluster> resultShow;
  flagAnsFound =false;
  //-----------------------------//
  //Initialize performance matrix//
  //-----------------------------//
  cur1st_ms = return1stTotalNcal1stResors_HomoPower();
  cur2nd_Joule = returnTransientJoule();
  cur1st_Joule = power1st*cur1st_ms/1000.0;

  curSupNum=cSystem->calSupNodes();
  curChNum=maxChNum;
  nextChNum=curChNum;
  curJEntropy = curSupNum*indEntropy + matrixComputer->computeLog2Det(1.0,cSystem->allSupStru);
  curPayoff=cur1st_ms;

  cout << "curPayoff: " << curPayoff << endl;
  bestAllServeFound=false;

  if ( checkBestClusterStructure_DataCentric(0) ) return true;
  cout << "Compression Ratio " << returnComprRatio() << 
    " indEntropy " << indEntropy << endl;
#ifdef ITERATION 
  fstream iterFile("guideIterObserve.out",ios::out);
#endif
  for(int i=0; i<SAIter; i++)
  {
    coolOnce_minResors();

    if( targetHeadIndex == -1 || targetHeadIndex == -1 ) {
      if ( nextEventFlag == 4 ) {
        passNext2Cur();
        for(int j=0; j<totalNodes; j++) 
          nodes[j].power = nextNodePower[j];
      }
      i++;
      continue;
    }

    calculateMatrics_minResors();
    confirmNeighbor3i();
    if (curJEntropy>(fidelityRatio*wholeSystemEntopy)) {
      flagAnsFound=true;
    }
    assert(curSupNum>=0);
    if(checkBestClusterStructure_DataCentric(i))
    {
      cout<<"Congratulation All nodes are served"<<endl;
      break;
    }

    if (isDetailOutputOn) {
      writePayoffEachRound_MinResors_withHead(i,curChNum);
    }
    int tempche=(signed)cSystem->listUnSupport->size();

    assert(( curSupNum + tempche) == totalNodes);
#ifdef ITERATION
    iterFile << curPayoff << endl;
#endif
  }
#ifdef ITERATION
  iterFile.close();
#endif

  end = clock();
  computingTimes = ((float)(end-begin))/CLOCKS_PER_SEC;
  cout << "best "<< bestFeasibleSupNum << "   Information Ratio:" 
    << bestFeasibleJEntropy / wholeSystemEntopy << endl;
  if(!flagAnsFound)
  {
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
  if(bestFeasibleJEntropy>=(wholeSystemEntopy*fidelityRatio)) {
    cout<<timeBuf<<endl;
#ifdef DEBUG
    if (isDetailOutputOn) {
      sprintf(str,"%s_Best4b2Struc2ndN%d_m%d_FR%.1f_r%.1f.txt",timeBuf,totalNodes,maxChNum,fidelityRatio,radius);
      resultShow.writeBestStru(str,*this);
      sprintf(str2,"%s_Detail4b2DataULSA3i%d_m%d_FR%.1f_r%.1f.txt",timeBuf, totalNodes,maxChNum,fidelityRatio,radius);
      resultShow.summaryNwrite2tiers_MinResors_with2ndPowerControl(str2, *this, best2nd_ms);
    }
    sprintf(str3,"tmpAll/ULSA4b7_All_N%d_BW%.1fPW%.3f_FR%.2f_r%.1f.%s.txt",totalNodes,bandwidthKhz,powerMax,fidelityRatio,radius, strIpAddr.c_str());
    resultShow.writePeformance_MinResors_with2ndPowerControl_4b(str3,*this, best2nd_ms, fidelityRatio,bestChNum);
    sprintf(str3,"tmpAll/ULSA4b7_Cluster_N%d_BW%.1fPW%.3f_FR%.2f_r%.1f.%s.txt",totalNodes,bandwidthKhz,powerMax,fidelityRatio,radius, strIpAddr.c_str());
    resultShow.writeClusterInfo(str3,*this,timeBuf);
#endif

    delete matrixComputer;
    return true;
  }
  else{
#ifdef DEBUG
    if (isDetailOutputOn) {
      sprintf(str,"%s_Best4b2Struc2ndN%d_m%d_FR%.1f_r%.1f.txt",timeBuf,totalNodes,maxChNum,fidelityRatio,radius);
      resultShow.writeBestStru(str,*this);
      sprintf(str2,"%s_Detail4b2DataULSA3i%d_m%d_FR%.1f_r%.1f.txt",timeBuf, totalNodes,maxChNum,fidelityRatio,radius);
      resultShow.summaryNwrite2tiers_MinResors_with2ndPowerControl(str2, *this, best2nd_ms);
    }
#endif
    delete matrixComputer;
    return false;
  }

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
void MinPowerSACluster::coolOnce_minResors()
{
  int probAdd = ((curSupNum<(totalNodes)) ?20000 :0);
  int probDiscard = ((curSupNum<(maxChNum+1)) ?0:30000);
//  int probAdd = 0;
//  int probDiscard = 0;
  int probExchange = 0;
  bool checkRotateble=false;//check if there are rotatableSet;
  for(int i=0; i<maxChNum; i++) {
    //cout<<aryFlagHRDone[i]<<" "<<cSystem->vecClusterSize[i]<<endl;
    if( aryFlagHRDone[i] == false && cSystem->vecClusterSize[i] > 1 )
      checkRotateble=true;
  }
  int probHeadRotate = (((curSupNum>2*maxChNum)&&checkRotateble) ?10000:0);//Don't do head rotate if there are only a few nodes



  int tmpJoinCan=0;
  bool chkLessCluster=cSystem->returnIfClusterSmall(thresholdd,tmpJoinCan);

  //int probJoin = (chkLessCluster&&curJEntropy>(fidelityRatio*wholeSystemEntopy))?tmpJoinCan*50:0;
  int probJoin = (chkLessCluster)?tmpJoinCan*30:0;

  //probJoin=((lastJoinPassAccu>thres2-400)?probJoin:0);

  int probIsoltae=((curChNum<maxChNum)?30:0);
  //probIsoltae=((lastJoinPassAccu>thres2)?probIsoltae:0);


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
      decideAdd3i_DC_HeadDetMemRan();
      if(targetHeadIndex!=-1&&targetNode!=-1){
        //cout<<"add "<<targetNode<<" to "<<cSystem->vecHeadName[targetHeadIndex]<<endl;
        addMemberSA(targetHeadIndex,targetNode);
      }
    }
    nextChNum=curChNum;

  }
  else if (nextEventFlag ==2)
  {
    decideDiscard3b();
    discardMemberSA(targetHeadIndex,targetNode);
    nextChNum=curChNum;

    //cout<<"discard "<<targetNode+1<<" from "<<cSystem->vecHeadName[targetHeadIndex]+1<<endl;
  }
  else if (nextEventFlag ==3)
  {
    decideHeadRotate2i_DC_HeadRanMemDet();
    //cout<<"HR "<<targetNode+1<<" to Replace "<<cSystem->vecHeadName[targetHeadIndex]+1<<endl;
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
MinPowerSACluster::decideAdd3i_DC_HeadDetMemRan() {
    targetHeadIndex=-1;
    targetNode=-1;

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
    assert(targetNode != -1 && targetHeadIndex != -1);
}

void MinPowerSACluster::decideHeadRotate2i_DC_HeadRanMemDet()
{
	//----Uniformly choosed
	int rotateAbleSize=0;
	for (int i=0; i<maxChNum; i++)
		if(cSystem->vecClusterSize[i]>1)rotateAbleSize++;
	//if(cSystem->vecClusterSize[i]>1&&aryFlagHRDone[i]==false)rotateAbleSize++;

	//Notice: no handle of "rotateAbleSize == 0", it handle by coolOnce_minResors
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
	double temp1st_ms=0;
	double temp2nd_ms=0;
	double temp2tiers_ms=0;
	double test2tiers_ms=DBL_MAX;

	int OriginalNode=cSystem->vecHeadName[targetHeadIndex];
	list<int>::iterator it1=itlist1->begin();
	for(; it1!=itlist1->end(); it1++)//find minimum interference received among nodes in the cluster
 	{
		if((*it1)==OriginalNode)continue;
		cSystem->vecHeadName[targetHeadIndex]=(*it1);//==Head Rotate
		temp2nd_ms=0;
		temp1st_ms=return1stTotalNcal1stResors_HomoPower();
		temp2tiers_ms=temp1st_ms+temp2nd_ms;
		if(temp2tiers_ms<test2tiers_ms)
		{
			test2tiers_ms=temp2tiers_ms;
			targetNode=*it1;
		}
	}
	cSystem->vecHeadName[targetHeadIndex]=OriginalNode;
}


void MinPowerSACluster::decideHeadJoining4b(){
    int threshold=thresholdd; // To determine the cluster size is large enough 
                              // to perform join operation 
    JoiningHeadIndex=-1;
    targetHeadIndex=-1;

}

// -------------------------------------------------------------------------- //
// @Description: targetHeadIndex, IsolateNodeName
// @Provides: 
// -------------------------------------------------------------------------- //

void MinPowerSACluster::decideIsolate4b(){
    assert(isolatedHeadIndex!=-1&&IsolateNodeName!=-1);
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
     list<list <int> >::iterator it_LiInt=cSystem->listCluMember->begin();
     for(int i=0;i<JoiningHeadIndex;i++)it_LiInt++;
     list<int>::iterator it_Int=it_LiInt->begin();
     for(int i=0;i<it_LiInt->size();i++,it_Int++){
         lastJoingingMachine.push_back(*it_Int);
         nodes[*it_Int].ptrHead=&(cSystem->vecHeadName[targetH]);
     }

     cSystem->join_FromHead(JoiningHeadIndex,targetH);
     nextChNum=curChNum-1;

}


void MinPowerSACluster::calculateMatrics_minResors()//Calculate next performance matircs
{
  next2nd_ms = 0;
  next1st_ms = return1stTotalNcal1stResors_HomoPower();
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
  nextPayoff =(next1st_ms+next2nd_ms);
}

/*
   after do the neighbor change
   -confirm the neighbor change and decide whether reverse or not
   -confirm3c add reset some metric if structure change(add,discard)
   by"if(nextEventFlag==1||nextEventFlag==2)confirmStructureChange();"
   */
void MinPowerSACluster::confirmNeighbor3i()
{
  //constraint: the fidelity ratio must be satisfied
//  bool nextAllServe = (nextJEntropy>(fidelityRatio*wholeSystemEntopy)?true:false);
//  bool curAllServe = (curJEntropy>(fidelityRatio*wholeSystemEntopy)?true:false);
//
  //constraint: the all the machine must be supported
  bool nextAllServe = (find(cSystem->allSupStru, cSystem->allSupStru + totalNodes,false ) == cSystem->allSupStru + totalNodes ); 
  bool curAllServe = (find(prevAllSupStru, prevAllSupStru + totalNodes, false) == prevAllSupStru + totalNodes);
  if  ( ( nextJEntropy >= curJEntropy )&& !nextAllServe && !curAllServe )
  {
    passNext2Cur();
    for(int i=0; i<totalNodes; i++) 
      nodes[i].power = nextNodePower[i];
    if( nextEventFlag == 1 || nextEventFlag == 2 ) 
      confirmStructureChange();
  }
  else if ( !curAllServe && nextAllServe ) {
    passNext2Cur();
    for(int i=0; i<totalNodes; i++) 
      nodes[i].power = nextNodePower[i];
    if( nextEventFlag == 1 || nextEventFlag == 2 ) 
      confirmStructureChange();
  }
  else if ( ( nextPayoff < curPayoff ) && nextAllServe && curAllServe )
  {
    passNext2Cur();
    for(int i=0; i<totalNodes; i++)
      nodes[i].power = nextNodePower[i];
    if( nextEventFlag == 1 || nextEventFlag == 2 )
      confirmStructureChange();
  }
  else if( ( ( nextJEntropy < curJEntropy ) && !nextAllServe && !curAllServe ) || 
      ( !nextAllServe && curAllServe ) || 
      ( ( nextPayoff > curPayoff ) && nextAllServe && curAllServe ) )
  {
    double probAnnealing = exp (-20*abs(nextPayoff-curPayoff)/temparature);
    //cout<<"Show Payoff "<<nextPayoff<<"  "<<curPayoff<<endl;
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
  temparature*=alpha;
}

void MinPowerSACluster::passNext2Cur() {


    curJEntropy = nextJEntropy;
    curSupNum = nextSupNum;

    curChNum=nextChNum;

    curPayoff = nextPayoff;

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

    }
    else if(nextEventFlag == 5){
         cSystem->reverseisolate(lastIsolateNodeName,lastIsolatedClusterIndex,targetHeadIndex);
        nodes[lastIsolateNodeName].ptrHead=&cSystem->vecHeadName[lastIsolatedClusterIndex];
       /* list <list <int> >:: iterator itt=cSystem->listCluMember->begin();
       for(unsigned ii=0; itt!=cSystem->listCluMember->end(); ii++,itt++)
        { cout<<ii<<": ";
        list <int>::iterator itt2 = itt->begin();
        for(; itt2!=itt->end(); itt2++)
        {
            cout<<(*itt2)<<", ";
        }
        cout<<endl;

        }*/

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
}



/*
    Check if this structure is the best
*/
bool MinPowerSACluster::checkBestClusterStructure_DataCentric(int inputRound)
{
//    bool curAllServe = (curJEntropy>fidelityRatio*wholeSystemEntopy?true:false);
    bool curAllServe = (find(prevAllSupStru, prevAllSupStru + totalNodes, false) == prevAllSupStru + totalNodes);
    if(curAllServe)bestAllServeFound=true;
    //cout<<"In check Best"<<endl;
    if ((curJEntropy>bestFeasibleJEntropy)&&!curAllServe&&!bestAllServeFound)
    {
        roundBest = inputRound;
        bestFeasibleJEntropy=curJEntropy;
        bestFeasibleSupNum=curSupNum ;
        //bestFeasiblePayoff=curPayoff;
        best1st_Joule=cur1st_Joule;
        best1st_ms=cur1st_ms;
        best2nd_Joule=cur2nd_Joule;
        best2nd_ms=cur2nd_ms;
        keepBestStructure();
        bestChNum=curChNum;

    }
    else if ((curPayoff<bestFeasiblePayoff)&&curAllServe)
    {
        //cout<<"Find new best"<<endl;
        roundBest = inputRound;
        bestFeasibleJEntropy=curJEntropy;
        bestFeasibleSupNum=curSupNum ;
        bestFeasiblePayoff=curPayoff;
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
    listCluMemBest->assign(cSystem->listCluMember->begin(), cSystem->listCluMember->end());
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
int MinPowerSACluster::returnClosetNodeIndexInGroup(int tempX,int tempY, vector<int> &inputGroup)
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
    list <list<int> >::iterator itlist1=cSystem->listCluMember->begin();
    double accuJoule=0;
    for(int i =0; itlist1!=cSystem->listCluMember->end(); itlist1++,i++)
    {
        list<int>::iterator it1=itlist1->begin();
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

double MinPowerSACluster::return1stTotalNcal1stResors_HomoPower() {
    power1st=powerMax;
    double T1=0;
    list <list <int> >::iterator it_LiInt=cSystem->listCluMember->begin();
    for(unsigned int i=0; i<cSystem->vecHeadName.size(); i++,it_LiInt++) {
        if ( cSystem->vecHeadName[i] == -1 || it_LiInt->size() == 0 ) continue;
        double information = it_LiInt->size()*indEntropy+ matrixComputer->computeLog2Det(1.0, cSystem->clusterStru[i]);
        //cout<<i<<" "<<it_LiInt->size()<<" "<<information<<endl;
        double temp = (information/(bandwidthKhz*log2(1+power1st*Gib[cSystem->vecHeadName[i]]/realNoise)));
        vecClusterHeadBits[i]=information;
        vecClusterHeadMS[i]= temp;
        vecClusterHeadWatt[i]=power1st;
        T1+=temp;
    }
    return T1;
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



