#ifndef _MINRESCSFACTORY_
#define _MINRESCSFACTORY_

#include "csFactory.h"
#include "ULSA4b7_DC.h"

class MinResCsFactory: public CsFactory
{
  public:
    MinResCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~MinResCsFactory();
    ClusterStructure * CreateClusterStructure();
  private:
    bool Kmedoid( vector<int>&, list<list<int> >& );
};

#endif
