#include "imageCsFactory.h"


ImageCsFactory::ImageCsFactory( ImageMap const * const myMap, 
    CORRE_MA_OPE const * const myMatComputer):
  m_ptrMap(myMap), m_ptrMatComputer(myMatComputer),
  m_ptrCS(0)
{
  m_numNodes = m_ptrMap->GetNumNodes();
  m_numMaxHeads = m_ptrMap->GetNumInitHeads();
}

ImageCsFactory::~ImageCsFactory()
{
  if (m_ptrCS != NULL) {
    delete m_ptrCS;
  }
}


