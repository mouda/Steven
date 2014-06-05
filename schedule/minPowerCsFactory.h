#ifndef _MINPOWERCSFACTORY_
#define _MINPOWERCSFACTORY_

#include "csFactory.h"
#include "minPowerSACluster.h"

class MinPowerCsFactory: public CsFactory
{
  public:
    MinPowerCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~MinPowerCsFactory();
    void                    SetMapFileName(const string& fileName){ m_mapFileName = fileName;}
    void                    SetCompressionRatio( const double compressionRatio){ m_compressionRatio = compressionRatio;}
    void                    SetFidelityRatio( const double fidelityRatio){ m_fidelityRatio = fidelityRatio; }
    ClusterStructure*       CreateClusterStructure();
  private:
    bool Kmedoid( vector<int>&, list<list<int> >& );

    MinPowerSACluster*    m_ptrToolSA; 
    std::string           m_mapFileName;
    FILE*                 m_fid;
    double                m_compressionRatio;
    double                m_fidelityRatio;
};

#endif
