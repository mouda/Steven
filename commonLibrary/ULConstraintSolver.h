#ifndef ULConstraintSolver_H
#define ULConstraintSolver_H

#include <list>
#include <vector>
#include "ULAGENT.h"

/*
    ULConstraintSovler is designed for 2nd tier constraint and computation

*/
using namespace std;
class ULConstraintSolver
{
    public:
        //first Constructor is for feasibility check for the 2nd-tier for given cluster structure
        ULConstraintSolver(int inChNum,int inTotalNodes,float inPowerMax, float ininBandNoise,float inC2, vector<int>&  \
                           inVecHeadName ,float** inGij,double* inNextNodePower,list<list<int> >*inListCluMem);
        //Constructor for solving minimize T2 find
        ULConstraintSolver(int inChNum,int inTotalNodes,float inPowerMax, double ininBandNoise,double inBandwidthKhz, double inData \
                                       , vector<int>&inVecHeadName ,float** inGij,double* inNextNodePower,list<list<int> >*inListCluMem);

        ~ULConstraintSolver();

        //Debug Tool
        void printPowerRAll();
        void printPowerAll();
        void printPowerDifAll();
        //Display Tool
        void showVerificationResult(vector <double> &vecBestReceivedInterference,vector<double> &vecBestSINR_forVerification, \
                                    vector<double> &vecBestBpsHz_forVerification);

        int returnCriticalNode();
        int returnCriticalHeadIndex();
        //------
        //PublicIntergae//
        bool solve();
        int statusFlag;

        double solve_withT2Adj_BinerySearch(double iniT2);
        double solve_withT2Adj_BinerySearch_2(double iniT2);

        // solve with power control mechanism ad return T2: unfinished;
        double solve_withT2Adj(double iniT2);

        void returnMaxNextNodePower(int &nodeIndex, double &maxNextPower);
        void returnMaxPowerGradient();
        void updateT2_DataHomo(int &nodeIndex, double &maxNextPower, bool feasiblility);

        //solve with homogeneous power
        void solve_withPowerHomo(double &totalT2);



        //----
        double need1stResorsIndividualMs; // This value is also use for Headrotate

        void Do1stTierPowerMax_ForBest(double &returnJoule, double &return1sttierMS, \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS, vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float *Gib,vector<int>&vecBestHeadName );

        void Do1stTierPowerControl(double &returnJoule, double &return1sttierMS, \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS, vector<double> &vecClusterHeadBits, double *rateib_MaxPower,float  *Gib );

        void Do1stTierPowerControl_ForBest(double &returnJoule, double &return1sttierMS, \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS, vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float *Gib,vector<int>&vecBestHeadName );

        void Do1stTierPowerControl_ForBest_PibAdjust(double &returnJoule, double &return1sttierMS, \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS,vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float  *Gib,vector<int>&vecBestHeadName,double pibAdjust_Scale );


        bool Do1stTierPowerControl_1stLimited_EqDiv(double &returnJoule,  \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS, vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float  *Gib , double timeMs1st);


        bool Do1stTierPowerControl_ForBest_1stLimited_EqDiv(double &returnJoule,  \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS, vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float *Gib,vector<int>&vecBestHeadName,  double timeMs1st );




         void setinBand(double ininBand);

        int returnHeadIndex_ByNodeName(int NodeName);
        void updateInterference();
        void changeAllMemberPower();
        bool checkConverged();
        bool checkDifference();
        double returnInfeasibility();
        double returnAverageInterference();
        double singleNodeInformation_Bits;

        //private:
        //From Caller Function
        int chNum;
        int totalNodes;
        float powerMax;
        float inBandNoise;
        float C2;

        float ** Gij;
        int **maIndexInterference;
        double **maStrengthInterference;
        vector <int> *vecHeadName;
        double *nextNodePower;
        list<list<int> > *listCluMember;
        vector <ULAGENT> *nodes;
        int criticalNode;
        int criticalHeadIndex;
    private:
        //For Contructor 2-nd---
        double bandwidthKhz;
        double best2nd_ms;
        double data_homo;

        //----------------------
        //Internal Variable
        double* powerDifference;
        double* powerDifferenceRatio;
        double avgRatio;
        bool exceedPc;
        static const double scale = 1; //This scale is for avoiding computation error.
};
#endif
//Not Applied
//void changeSingleClusterMemberPower(int headIndex,list <list<int> >::iterator it1);
//void updateSingleClusterInterference(int headIndex2,list <list<int> >::iterator it2);
