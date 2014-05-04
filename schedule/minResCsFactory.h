#ifndef _MINRESCSFACTORY_
#define _MINRESCSFACTORY_

#include "csFactory.h"

class MinResCsFactory: public CsFactory
{
  public:
    MinResCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~MinResCsFactory();
  private:
    ClusterStructure*   m_ptrCS;
    Map const * const   m_ptrMap;
    CORRE_MA_OPE const * const  m_ptrMatComputer;
};

#endif
