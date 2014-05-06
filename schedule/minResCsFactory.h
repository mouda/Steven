#ifndef _MINRESCSFACTORY_
#define _MINRESCSFACTORY_

#include "csFactory.h"
#include "ULAGENT.h"
#include "csPowerUpdater.h"

class MinResCsFactory: public CsFactory
{
  public:
    MinResCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~MinResCsFactory();
    ClusterStructure*   CreateClusterStructure();
  private:

    bool                SASearch();
    double              Return1stTotalNcal1stResors_HomoPower(); 
    double              returnTransientJoule( const vector<double>& vecPower );
    bool                checkBestClusterStructure_DataCentric(int inputRound, CSPowerUpdater& myPowerUpdater,  vector<double>& vecPower);
    void                keepBestStructure( const vector<double>& vecPower);
    void                coolOnce_minResors( CSPowerUpdater& myPowerUpdater, vector<double>& vecPower);
    bool                Kmedoid( vector<int>& vecHeadNames, list<list<int> >& listCluMembers );

    void                decideAdd3i_DC_HeadDetMemRan();
    void                addMemberSA(int inputHeadIndex, int inputMemberName);
    void                decideDiscard3b();
    void                discardMemberSA(int inputHeadIndex2, int inputMemberName2);
    void                decideHeadRotate2i_DC_HeadRanMemDet( CSPowerUpdater& myPowerUpdater, vector<double>& vecPower );
    void                rotateHeadSA(int inputHeadIndex3, int inputMemberName3);


    vector<vector<int> >m_matBestCluStru;
    list<list<int> >    m_listCluMemBest;
    vector<double>      m_powerBest;
    vector<int>         m_vecBestAllSupStru;
    vector<double>      m_vecClusterHeadBits;
    vector<double>      m_vecClusterHeadWatt;
    vector<double>      m_vecClusterHeadMS;
    vector<double>      m_vecBestClusterBits;
    vector<int>         m_vecBestClusterSize;
    vector<double>      m_vecBestClusterHeadMS;
    vector<double>      m_vecBestClusterHeadWatt;
    vector<int>         m_vecHeadNameBest;

    vector<double>      m_vecBestReceivedInterference;
    vector<double>      m_vecBestSINR_forVerification;
    vector<double>      m_vecBestBpshz_forVerification;
    vector<double>      m_vecChooseIndex;
    vector<ULAGENT>     m_nodes;

    int                 m_roundBest;
    double              m_fidelityRatio;// Temporary set by here 2013/02/21
    double              m_power1st;

    int                 m_nextChNum;
    int                 m_curChNum;
    double              m_curPayoff;
    double              m_curJEntropy;
    int                 m_curSupNum;
    double              m_cur1st_ms;
    double              m_cur1st_Joule;
    double              m_cur2nd_ms;
    double              m_cur2nd_Joule;


    double              m_nextPayoff;
    double              m_nextJEntropy;
    int                 m_nextSupNum;
    double              m_next1st_ms;
    double              m_next1st_Joule;
    double              m_next2nd_ms;
    double              m_next2nd_Joule;


    double              m_wholeSystemEntropy;
    double              m_indEntropy;

    //@BEST KPI
    bool                m_bestAllServeFound;
    double              m_best1st_Joule;
    double              m_best2nd_Joule;
    double              m_best1st_ms;
    double              m_best2nd_ms;
    double              m_bestFeasibleJEntropy;
    double              m_bestFeasiblePayoff;
    int                 m_bestFeasibleSupNum;
    int                 m_bestChNum;

    int                 m_targetHeadIndex;
    int                 m_targetNode;
    int                 m_JoiningHeadIndex;
    int                 m_IsolateNodeName;
    int                 m_isolatedHeadIndex;
    int                 m_nextEventFlag;
    bool                m_iniDone;
    static const int    m_thresholdd = 5;
    int*                m_ptrHeadLastDiscard;
    double              m_powerLastDiscard;
    int                 m_rotatedHeadNameLast;
    vector<bool>        m_aryFlagHRDone; //If it is true means Head Rotate have been done in this structure, starrt from ULSA4b7_DC
};

#endif
