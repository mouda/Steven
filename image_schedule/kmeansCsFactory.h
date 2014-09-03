#ifndef _KMEANSCSFACTORY_
#define _KMEANSCSFACTORY_

#include "imageCsFactory.h"
#include "imageSource.h"

class KmeansCsFactory: public ImageCsFactory
{
  public:
    KmeansCsFactory(ImageMap const * const, CORRE_MA_OPE const * const );
    KmeansCsFactory(ImageMap const * const, ImageSource const * const );
    ~KmeansCsFactory();
    ClusterStructure * CreateClusterStructure();
  private:
    bool Kmedoid( vector<int>&, list<list<int> >& );
};

#endif
