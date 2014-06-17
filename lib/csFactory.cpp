#include "csFactory.h"


CsFactory::CsFactory( Map const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  m_ptrMap(myMap), m_ptrMatComputer(myMatComputer),
  m_ptrCS(0)
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads = m_ptrMap->GetNumInitHeads();
}

CsFactory::~CsFactory()
{
  if (m_ptrCS != NULL) {
    delete m_ptrCS;
  }
}


