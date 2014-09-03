#ifndef _CSFACTORY_
#define _CSFACTORY_
#include <limits>
#include <armadillo>
#include <vector>
#include <list>

#include "clusterStructure.h"
#include "CORRE_MA_OPE.h"
#include "imageSource.h"
#include "imageMap.h"

using std::vector;
using std::list;
using std::pair;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::fstream;
using std::ios;
using std::stringstream;
using std::make_pair;

class ImageCsFactory
{
  public:
    ImageCsFactory(ImageMap const * const, CORRE_MA_OPE const * const );
    virtual ~ImageCsFactory();
    virtual ClusterStructure * CreateClusterStructure() = 0;
  protected:

    int                 m_numNodes;
    int                 m_numMaxHeads;
    ClusterStructure*   m_ptrCS;
    ImageMap const * const   m_ptrMap;
    CORRE_MA_OPE const * const  m_ptrMatComputer;
};
#endif
