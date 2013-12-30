#ifndef _CSFACTORY_
#define _CSFACTORY_
#include "clusterStructure.h"

class CsFactory
{
  public:
    ClusterStructure * CreateClusterStructure();
  private:
    ClusterStructure * m_ptrCS;

};
#endif
