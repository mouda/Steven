#ifndef _MINRESCSFACTORY_
#define _MINRESCSFACTORY_

#include "csFactory.h"
#include "csPowerUpdater.h"

class MinResCsFactory: public CsFactory
{
  public:
    MinResCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~MinResCsFactory();
    ClusterStructure * CreateClusterStructure();
  private:

    bool    SASearch();
    double  Return1stTotalNcal1stResors_HomoPower(); 
    bool    Kmedoid( vector<int>& vecHeadNames, list<list<int> >& listCluMembers );

    double  m_wholeSystemEntropy;


    vector<double>  m_vecClusterHeadBits;
    vector<double>  m_vecClusterHeadWatt;
    vector<double>  m_vecClusterHeadMS;
};

#endif
