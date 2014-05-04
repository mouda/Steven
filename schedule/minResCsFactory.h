#ifndef _MINRESCSFACTORY_
#define _MINRESCSFACTORY_

#include "csFactory.h"

class MinResCsFactory: public CsFactory
{
  public:
    MinResCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~MinResCsFactory();
    ClusterStructure * CreateClusterStructure();
  private:

    bool SASearch(vector<int>& vecHeadNames, list<list<int> >& listCluMembers );
};

#endif
