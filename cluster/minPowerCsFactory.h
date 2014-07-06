#ifndef _MINPOWERCSFACTORY_
#define _MINPOWERCSFACTORY_

#include <iostream>
#include "csFactory.h"
#include "minPowerSACluster.h"

class MinPowerCsFactory: public CsFactory
{
  public:
    MinPowerCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~MinPowerCsFactory();
    void                    SetMapFileName(const std::string& fileName){ m_mapFileName = fileName;}
    void                    SetCompressionRatio( const double compressionRatio){ m_compressionRatio = compressionRatio;}
    void                    SetFidelityRatio( const double fidelityRatio){ m_fidelityRatio = fidelityRatio; }
    void                    SetTier2NumSlot(const int tier2NumSlot){ m_tier2NumSlot = tier2NumSlot; }
    void                    SetTier1TxTime( const double tier1NumSlot ){ m_tier1TxTime = tier1NumSlot;}
    void                    SetIterationLog( const bool flag ){ m_logFlag = flag; }
    ClusterStructure*       CreateClusterStructure();
  private:
    bool Kmedoid( std::vector<int>&, std::list<std::list<int> >& );

    MinPowerSACluster*    m_ptrToolSA; 
    std::string           m_mapFileName;
    FILE*                 m_fid;
    double                m_compressionRatio;
    double                m_fidelityRatio;
    int                   m_tier2NumSlot;
    double                m_tier1TxTime;
    bool                  m_logFlag;
};

#endif
