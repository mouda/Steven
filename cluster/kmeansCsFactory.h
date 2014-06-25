#ifndef _KMEANSCSFACTORY_
#define _KMEANSCSFACTORY_

#include "csFactory.h"

class KmeansCsFactory: public CsFactory
{
  public:
    KmeansCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~KmeansCsFactory();
    ClusterStructure * CreateClusterStructure();
  private:
    bool Kmedoid( vector<int>&, list<list<int> >& );
};

#endif
