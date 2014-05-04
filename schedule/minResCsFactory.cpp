#include "minResCsFactory.h"

MinResCsFactory::MinResCsFactory( Map const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  m_ptrMap(myMap), m_ptrMatComputer(myMatComputer),
  m_ptrCS(0),
  CsFactory(myMap, myMatComputer) 
{

}

MinResCsFactory::~MinResCsFactory()
{

}
