#ifndef _MINRESCSFACTORY_
#define _MINRESCSFACTORY_

#include "csFactory.h"
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
    bool                Kmedoid( vector<int>& vecHeadNames, list<list<int> >& listCluMembers );


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
};

#endif
