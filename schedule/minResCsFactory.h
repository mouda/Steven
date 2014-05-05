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

    bool SASearch();
    bool Kmedoid( vector<int>& vecHeadNames, list<list<int> >& listCluMembers );

    double m_wholeSystemEntropy;
};

#endif
