#ifndef _FILECSFACTORY_
#define _FILECSFACTORY_

#include "imageCsFactory.h"

class FileCSFactory: public ImageCsFactory
{
  public:
    FileCSFactory(ImageMap const * const, CORRE_MA_OPE const * const );
    FileCSFactory(ImageMap const * const, ImageSource const * const );
    ~FileCSFactory();
    ClusterStructure * CreateClusterStructure();
    void                SetFileName( const std::string& FName){ m_CSFName = FName; }
  private:
    bool ReadCSFromFile(std::vector<int>& , std::list<list<int> >& );

    std::string m_CSFName;
};

#endif
