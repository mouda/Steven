/*
  File: ULSA4b7_NOGUIDE.h
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
#ifndef ULSA4b7_NOGUIDE_H
#define ULSA4b7_NOGUIDE_H

#include<cmath>
#include<vector>
#include<cstdlib>
#include<list>
#include<ctime>
#include<map>

#include "../commonLibrary/ULAGENT.h"
#include "../commonLibrary/ULCS1b.h"
#include "../commonLibrary/ULConstraintSolver.h"
#include "../commonLibrary/CORRE_MA_OPE.h"
#include "../commonLibrary/SimSystem.h"
#include "../commonLibrary/ULSAOutputToolSet.h"
#include "../commonLibrary/ULSAOutputToolSet.cpp"
#include "../commonLibrary/TimeStamp.h"
#include "../lib/SA/SABASE.h"

template class ULSAOutputToolSet<class ULSA4b7_NOGUIDE>;

class ULCS1b;
class ULAGENT;
class ULSA4b7_NOGUIDE :public SABASE
{
public:
  ULSA4b7_NOGUIDE();
  //Constructor for scenario wise parameter
  ULSA4b7_NOGUIDE(FILE *fileReadCursor, int inputTotalNodes, int inputMaxChNum,int inputSAFac,  
               int inOutputControl,
               int isStrucOuput,
               double inputTemprature, double InputSaAlpha, 
               double inCorrelationFactor, string ipAddr);

  ~ULSA4b7_NOGUIDE();


  //Make a friend with template function
  template <class T>friend class ULSAOutputToolSet;
  //alpha is SA parameter in the paper.
  void releaseMemory();
  bool terminated;

//-------------------------------------------------------------------//
// @Purpose: Public Function
// @Called: by main
//-------------------------------------------------------------------//
  bool setSystem(float inPowerMaxWatt, int inQuantizationBits,double inBandwidthKhz, double fidelity);

  bool setInitialStucture(char* inputFlag);
  bool setIniStruKmeans();//not public but related to setIniStrucKmeans
  bool setIniStruDistanceKmedoids();
  bool setIniStruHalfResourceKmedoids();
  bool setIniStruFullResourceKmedoids();
  bool setIniHeadLimited();


  void writeStruSingleRound(int round);
  void writePayoffEachRound_MinResors(int round);
  void writePayoffEachRound_MinResors_withHead(int round,int head);

  //Incomplete!!!

  void do1sttierPowerControlforNext_DataCentric();
  void do1sttierPowerControlforBest_DataCentric();
  void do1sttierPowerControlforTEMP_DataCentric(double &temp1stJoule, double &temp1sMS);
  void do1sttierPowerMaxforBest_DataCentric();
  //

  void debug_CheckSizeCorrect();
  void computeBestTRR_DataCentric();
  void computeUpperResourceNoCodingNoPowerControl();
  void computeBestAvgInterference();
  void computeBestAvgPower();
  double computeRate2Nodes(int i, int j);
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
    const ULSA4b7_NOGUIDE& mySA;// data member
    compareDis(const ULSA4b7_NOGUIDE& sa) : mySA(sa) {} // constructor
    bool operator() (const int &i, const int &j)
    { //i, j are index
      return mySA.Gij[mySA.sortIndex][i] > mySA.Gij[mySA.sortIndex][j];
    }
  };

  //---------------------//
  //Module               //
  //---------------------//
  ULConstraintSolver* consSol;
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
  bool *aryFlagHRDone; //If it is true means Head Rotate have been done in this structure, starrt from ULSA4b7_NOGUIDE


 //--------------------------------------------//
  //Current,Next-Need Initialization and Passing//
  //--------------------------------------------//
  int nextChNum;
  int curChNum;
  double curPayoff;
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
  void do1sttierPowerControlforCur_DataCentric();
  void coolOnce_minResors();

  void adaptSearchingProbability();
  void updateHeadLocation();

  void decideAdd3i_DC_HeadDetMemRan();
  void decideAddRandom();
  void decideDiscard3b();
  void decideDiscardRandom();
  void decideDiscard3i_DC_HeadRanNodeDet_CompressionRatio();
  void decideHeadRotate1();
  void decideHeadRotate2f();
  void decideHeadRotate3c();
  void decideHeadRotate2i_DC_HeadRanMemDet();

  void decideHeadJoining4b();
  void decideIsolation4b();
  double estimateJoin2ndTierCost(int JoiningHeadIndex, int testIndex);
  void computeOriInterference_GivenTarInJoinI(std::vector<double> &oriInterf,std::vector<double> &interfExcept,int JoiningI,int tarHead);
  void updateJoinEstimatedPower(std::vector<double> &newPower, std::vector<int>&newMem,int JoiningHeadIndex, int targetIndex);
  void computeNewInterference_FromNewTarHI(std::vector<double> &newInterf,std::vector<double>&newPower,std::vector<int>&newMem,int JoiningHeadIndex, int targetIndex);

  void decideIsolate4b();
  double estimateIsolate2ndTierGain(int IsolatNodeName,int isoCluIndex);
  void computeOriInterference_GivenIsolate(std::vector<double> &oriInterf,std::vector<double> &interfExcept,int isoCluIndex);
  void updateIsolateEstimatedpower(std::vector<double> &newPower,int IsolateName,int isoCluIndex);
  void computeNewInterference_FromIsoCluster(std::vector<double> &newInterf,std::vector<double>&newPower,int IsolateName,int isoCluIndex);



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
  double returnTransientAveragePower();
  double returnTransientJoule();
  double return1stTotalNcal1stResors_HomoPower();


  //Display Function
  //void showAll
  //-----------------//
  //Judgemnt Database//
  //-----------------//
  float **distanceOf2Nodes;
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
};
#endif
