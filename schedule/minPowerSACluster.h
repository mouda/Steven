/*
  File: MinPowerSACluster.h
  Brief: Simulated anealing Solver. The purpose is to do the SA to solve the clustering problem with input
        "Mapfile position cursor" "Total nodes number" "Max number of CH allowed" "Power Max"
        Deal with Gaussian Random Source:
        1.C2 from float to double
  Author: Steven
  version:1d
  Date: 2012/09/25
  Latest Progress:

  Abbreviation:


*/
#ifndef MinPowerSACluster_H
#define MinPowerSACluster_H

#include <cmath>
#include <vector>
#include <cstdlib>
#include <list>
#include <ctime>
#include <map>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>

#include "ULAGENT.h"
#include "ULCS1b.h"
#include "CORRE_MA_OPE.h"
#include "simSystem.h"
#include "ULSAOutputToolSet.h"
#include "ULSAOutputToolSet.cpp"
#include "TimeStamp.h"
#include "../lib/SA/SABASE.h"
#include "MyNLP.h"
#include "map.h"

template class ULSAOutputToolSet<class MinPowerSACluster>;

class ULCS1b;
class ULAGENT;
class MinPowerSACluster :public SABASE
{
public:
  //Constructor for scenario wise parameter
  MinPowerSACluster(
      FILE *fileReadCursor, 
      int inputTotalNodes, 
      int inputMaxChNum, 
      int inputSAFac,  
      int inOutputControl,
      int isStrucOuput,
      double inputTemprature, 
      double InputSaAlpha, 
      double inCorrelationFactor, 
      string ipAddr, 
      Map const * const myPtrMap,
      double tier1TxTime
      );

  ~MinPowerSACluster();


  //Make a friend with template function
  template <class T>friend class ULSAOutputToolSet;
  //alpha is SA parameter in the paper.
  void releaseMemory();
  bool terminated;

//-------------------------------------------------------------------//
// @Purpose: Public Function
// @Called: by main
//-------------------------------------------------------------------//
  bool                      setSystem(float inPowerMaxWatt, int inQuantizationBits,double inBandwidthKhz, double fidelity);
  bool                      setInitialStucture(char* inputFlag);
  bool                      setIniStruKmeans();//not public but related to setIniStrucKmeans
  bool                      setIniStruDistanceKmedoids();
  bool                      setIniStruFullResourceKmedoids();
  bool                      setIniHeadLimited();

  /* For Integrate */
  const vector<int>&        GetVecHeadName() const { return vecHeadNameBest; }
  const list<list<int> >&   GetListCluMemeber() const { return *listCluMemBest;}
  vector<int>               GetAllSupStru() const;

  int                       OptimalRateControl( vector<double>& vecRate ) const;
  
  double returnComprRatio();

  // Main part
  bool startCool();

//-------------------------------------------------------------------//
// @Purpose: Public Parameters
// @Called: by main
//-------------------------------------------------------------------//

  //@Set by setSystem();
  double bandwidthKhz;
  float powerMax;
  double power1st;
  int quantizationBits;
  double dataBits;
  double virtualCompression;
  int headCandidatesNum;
  vector<int> vecHeadCandidates;


  //We use the struct as class
  int sortIndex;
  struct compareDis
  {
    const MinPowerSACluster& mySA;// data member
    compareDis(const MinPowerSACluster& sa) : mySA(sa) {} // constructor
    bool operator() (const int &i, const int &j)
    { //i, j are index
      return mySA.Gij[mySA.sortIndex][i] > mySA.Gij[mySA.sortIndex][j];
    }
  };

  //---------------------//
  //Module               //
  //---------------------//
  CORRE_MA_OPE * matrixComputer;
  double correlationFactor;
  double compRatio;

  std::vector <ULAGENT> nodes;

  //-----------------------//
  //SA Movement Variable   //
  //-----------------------//
  int  *ptrHeadLastDiscard;
  double powerLastDiscard;
  int lastJoinPassAccu;
  int lastIsoPassAccu;


  int rotatedHeadNameLast;

  bool flagAnsFound;


  int  targetHeadIndex;
  int  targetNode;
  int  JoiningHeadIndex;
  int  IsolateNodeName;
  int  isolatedHeadIndex;
  int  nextEventFlag;
  double* nextNodePower;



  bool iniDone;
  bool *aryFlagHRDone; //If it is true means Head Rotate have been done in this structure, starrt from MinPowerSACluster


 //--------------------------------------------//
  //Current,Next-Need Initialization and Passing//
  //--------------------------------------------//
  int nextChNum;
  int curChNum;
  double m_curPayoff;
  double  curJEntropy;
  int  curSupNum;
  //bool curFeasible;
  //double curInfeasibility;
  double cur1st_ms;
  double cur1st_Joule;
  double cur2nd_ms;
  double cur2nd_Joule;


  double nextPayoff;
  double  nextJEntropy;
  int  nextSupNum;
  //bool nextFeasible;
  //double  nextInfeasibility;
  double next1st_ms;
  double next1st_Joule;
  double next2nd_ms;
  double next2nd_Joule;


  double wholeSystemEntopy;
  double indEntropy;

  //@BEST KPI
  bool   bestAllServeFound;
  double best1st_Joule;
  double best2nd_Joule;
  double best1st_ms;
  double best2nd_ms;
  double bestFeasibleJEntropy;
  double bestFeasiblePayoff;
  int    bestFeasibleSupNum;
  int    bestChNum;

  vector<double> vecBestClusterBits;
  vector<int> vecBestClusterSize;
  vector<double> vecBestClusterHeadMS;
  vector<double> vecBestClusterHeadWatt;

  vector<double> vecBestReceivedInterference;
  vector<double> vecBestSINR_forVerification;
  vector<double> vecBestBpshz_forVerification;
  vector<double> vecChooseIndex;
  double best1stTierTraffic;
  double best2ndTierTraffic;
  double bestTrafficReductionRatio;



  //@Best Parameter
  double **maBestInterference;
  int **maBestInterfernceIndex;
  double bestAvgInterference;
  double bestAvgPowerOFAllNodes;
  double bestUpperLayerResource;
  //@Best Structure
  vector <int> vecHeadNameBest;
  list <list <int> > *listCluMemBest;
  bool** bestMaClusterStru;
  double ** maStrengthInterBest;
  double *powerBest;
  int  roundBest;;
  bool *bestAllSupStru;
  bool *prevAllSupStru;
  clock_t begin, end;
  float computingTimes;

  double radius;


  private:

//-------------------------------------------------------------------//
// @Purpose: Private Function
// @Called: Internal
//-------------------------------------------------------------------//
  //After Start Cool
  void addMemberSAIni(int inputHeadIndex0, int inputMemberName0);
  void coolOnce_minResors();

  void adaptSearchingProbability();
  void updateHeadLocation();

  void decideExchangeNode();
  void decideAdd3i_DC_HeadDetMemRan();
  void decideDiscard3b();
  void decideHeadRotate2i_DC_HeadRanMemDet();

  void decideHeadJoining4b();

  void decideIsolate4b();



  void addMemberSA(int inputHeadIndex, int inputMemberName);
  void discardMemberSA(int inputHeadIndex2, int inputMemberName2);
  void rotateHeadSA(int inputHeadIndex3, int inputMemberName3);

  void isolateHeadSA(int inputMemberName4,int IsolateCluI, int targetH);
  void join_fromHeadSA(int JoiningHeadIndex,int targetH);
  int lastJoiningHead;
  int lastJoiningHeadIndex;

  int lastIsolateNodeName;
  int lastIsolatedClusterIndex;
  void calculateMatrics_minResors();


  vector<int> lastJoingingMachine;
  void confirmNeighbor3i();
  void passNext2Cur();
  void reverseMoveSA();
  void confirmStructureChange();



  // Internal Aid Function
  bool checkBestClusterStructure_DataCentric(int inputRound);
  void keepBestStructure();
  int  returnClosetNodeIndexInGroup(int tempX,int tempY, std::vector <int> &inputGroup);
  double returnTransientJoule();


  //Display Function
  //void showAll
  //-----------------//
  //Judgemnt Database//
  //-----------------//
  double **distanceOf2Nodes;
  float **Gij;// channel gain from node i to node j
  float *Gib; // channel gain from node i to base station
  double *rateibMax; //
  double realNoise;

  void resetSA3iSystem();
  vector<double> vecClusterHeadBits;
  vector<double> vecClusterHeadWatt;
  vector<double> vecClusterHeadMS;
  SimSystem *sysComputing;
  map <int, int>  *mapNodeName2DistanceRank;
  int ** maIndexSortDecGain;
  ULCS1b*  cSystem;   // system cluseter structure
  double fidelityRatio;// Temporary set by here 2013/02/21
  double tempAddT;
  double tempDisT;
  double tempHRT;
  double tempIsoT;
  double tempJoinT;
  static const int thresholdd=5;
  //static const int thres2=800;
  static const int thres2=200;

  // display control
  int isDetailOutputOn;
  string strIpAddr; 

  vector<int>             m_vecHeadName;

  Map const * const       m_ptrMap;
  const double            m_tier1TxTime; 
};
#endif
