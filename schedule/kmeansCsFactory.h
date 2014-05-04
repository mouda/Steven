#ifndef _KMEANSCSFACTORY_
#define _KMEANSCSFACTORY_

#include "csFactory.h"

class KmeansCsFactory: public CsFactory
{
  public:
    KmeansCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~KmeansCsFactory();
  private:
    ClusterStructure*   m_ptrCS;
    Map const * const   m_ptrMap;
    CORRE_MA_OPE const * const  m_ptrMatComputer;
};

#endif
