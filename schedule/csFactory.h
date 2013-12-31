#ifndef _CSFACTORY_
#define _CSFACTORY_
#include "clusterStructure.h"
#include "CORRE_MA_OPE.h"
#include "map.h"

class CsFactory
{
  public:
    CsFactory();
    CsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~CsFactory();
    ClusterStructure * CreateClusterStructure();
  private:
    void Kmedoid();

    int                 m_numNodes;
    int                 m_numMaxHeads;
    ClusterStructure*   m_ptrCS;
    Map const * const   m_ptrMap;
    CORRE_MA_OPE const * const  m_ptrMatComputer;

};
#endif
