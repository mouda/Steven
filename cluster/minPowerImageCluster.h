/*
  File: MinPowerImageCluster.h
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
#ifndef MinPowerImageCluster_H
#define MinPowerImageCluster_H

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
#include "imageSource.h"
#include "simSystem.h"
#include "TimeStamp.h"
#include "../lib/SA/SABASE.h"
#include "myTier1NLP.h"
#include "map.h"


bool pairCompare(const std::pair<int,double>& lPair, const std::pair<int,double>& rPair);

class ULCS1b;
class ULAGENT;
class MinPowerImageCluster :public SABASE
{
public:
  //Constructor for scenario wise parameter
  MinPowerImageCluster(
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
      ImageSource const * const,
      double tier1TxTime,
      int tier2NumSlot,
      bool  logFlag
      );

  ~MinPowerImageCluster();


  //Make a friend with template function
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
  bool                      setIniGraphPartition();
  bool                      setIniBanancedModelCluster();
  double                    GetNodeDistance( const int lName, const int rName);
  bool                      CheckTwoLinkFeasible(const int lChName, const int lName, const int rChName, const int rName);
  bool                      CheckLinkFeasible(const int chName, const int name);
  bool                      CheckAllFeasible();
  double                    GetTier2ExpectPower(const int Name,const int ChName);
  bool                      CheckTier2Feasible();

  double                    GetSizePenalty(const vector<double>& );  
  double                    GetTier2Penalty(const vector<double>& );
  double                    GetEntropyPenalty( const double );
  int                       GetSupportNodes();

  /* For Integrate */
  const std::vector<int>&         GetVecHeadName() const { return vecHeadNameBest; }
  const list<list<int> >&         GetListCluMemeber() const { return listCluMemBest;}
  std::vector<int>                GetAllSupStru() const;
  double                          GetBestPayoff() const { return bestFeasiblePayoff;}

  double                          OptimalRateControl() const;
  
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
  std::vector<int> vecHeadCandidates;


  //We use the struct as class
  int sortIndex;
  struct compareDis
  {
    const MinPowerImageCluster& mySA;// data member
    compareDis(const MinPowerImageCluster& sa) : mySA(sa) {} // constructor
    bool operator() (const int &i, const int &j)
    { //i, j are index
      return mySA.Gij[mySA.sortIndex][i] > mySA.Gij[mySA.sortIndex][j];
    }
  };

  //---------------------//
  //Module               //
  //---------------------//
  ImageSource const * const m_ptrImageSource;
    
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
  bool *aryFlagHRDone; //If it is true means Head Rotate have been done in this structure, starrt from MinPowerImageCluster


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

  double m_cur1st_watt;


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

  std::vector<double> vecBestClusterBits;
  std::vector<int> vecBestClusterSize;
  std::vector<double> vecBestClusterHeadMS;
  std::vector<double> vecBestClusterHeadWatt;

  std::vector<double> vecBestReceivedInterference;
  std::vector<double> vecBestSINR_forVerification;
  std::vector<double> vecBestBpshz_forVerification;
  std::vector<double> vecChooseIndex;
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
  std::vector <int> vecHeadNameBest;
  list <list <int> > listCluMemBest;
  bool** bestMaClusterStru;
  double ** maStrengthInterBest;
  double *powerBest;
  int  roundBest;;
  bool *bestAllSupStru;
  bool *prevAllSupStru;
  std::vector<int> m_prevVecClusterSize;
  std::vector<int> m_prevVecHeadName;
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
  void GetNeighbor1( const int iterSA);
  void GetNeighbor2( const int iterSA, const vector<double>&, const vector<double>&, const double&, 
      vector<double>&, vector<double>&, double& );
  double GetPayOff( const vector<double>&, const vector<double>&, const double& );
  void calculateMatrics_minResors(const vector<double>&, const vector<double>&, const double& );
  void ConfirmNeighbor1();
  void ConfirmNeighbor2(vector<double>&, vector<double>&, double&, const vector<double>&, const vector<double>&, const double& );

  void decideExchangeNode();
  void decideAddClosetAddableNode();
  void decideAdd3i_DC_HeadDetMemRan();
  void decideAddSmallestSize();
  void decideAddRandSelectCluster();
  void decideDiscard3b();
  void decideDiscard3o();
  void decideDiscard3f();
  void decideDiscardMinGain();
  void decideHeadRotate2i_DC_HeadRanMemDet();

  void decideHeadJoining4b();
  double GetClusterEntropy(const int name );
  double GetJoinHeadEntropy(const int lName, const int rName);

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


  std::vector<int> lastJoingingMachine;
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
  std::vector<double> vecClusterHeadBits;
  std::vector<double> vecClusterHeadWatt;
  std::vector<double> vecClusterHeadMS;
  SimSystem *sysComputing;
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

  std::vector<int>             m_vecHeadName;

  Map const * const       m_ptrMap;
  const double            m_tier1TxTime; 
  int                   m_tier2NumSlot;
  fstream               m_logFile; 
  bool                  m_logFlag;

};
#endif
