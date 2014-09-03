#ifndef _MINPOWERIMAGECSFACTORY_
#define _MINPOWERIMAGECSFACTORY_

#include <iostream>
#include "imageCsFactory.h"
#include "imageSource.h"
#include "minPowerImageCluster.h"

class MinPowerImageCsFactory: public ImageCsFactory
{
  public:
    MinPowerImageCsFactory(ImageMap const * const, CORRE_MA_OPE const * const );
    MinPowerImageCsFactory(ImageMap const * const , ImageSource const * const );
    ~MinPowerImageCsFactory();
    void                    SetMapFileName(const std::string& fileName){ m_mapFileName = fileName;}
    void                    SetCompressionRatio( const double compressionRatio){ m_compressionRatio = compressionRatio;}
    void                    SetFidelityRatio( const double fidelityRatio){ m_fidelityRatio = fidelityRatio; }
    void                    SetTier2NumSlot(const int tier2NumSlot){ m_tier2NumSlot = tier2NumSlot; }
    void                    SetTier1TxTime( const double tier1NumSlot ){ m_tier1TxTime = tier1NumSlot;}
    void                    SetIterationLog( const bool flag ){ m_logFlag = flag; }
    ClusterStructure*       CreateClusterStructure();
  private:
    bool Kmedoid( std::vector<int>&, std::list<std::list<int> >& );

    MinPowerImageCluster*    m_ptrToolSA; 
    std::string           m_mapFileName;
    FILE*                 m_fid;
    double                m_compressionRatio;
    double                m_fidelityRatio;
    int                   m_tier2NumSlot;
    double                m_tier1TxTime;
    bool                  m_logFlag;
    ImageSource const * const   m_ptrImageSource;
};

#endif
