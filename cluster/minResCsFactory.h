#ifndef _MINRESCSFACTORY_
#define _MINRESCSFACTORY_

#include "csFactory.h"
#include "ULSA4b7_DC.h"

class MinResCsFactory: public CsFactory
{
  public:
    MinResCsFactory(Map const * const, CORRE_MA_OPE const * const );
    ~MinResCsFactory();
    void                    SetMapFileName(const string& fileName){ m_mapFileName = fileName;}
    void                    SetCompressionRatio( const double compressionRatio){ m_compressionRatio = compressionRatio;}
    void                    SetFidelityRatio( const double fidelityRatio){ m_fidelityRatio = fidelityRatio; }
    void                    SetTier2NumSlot(const int tier2NumSlot){ m_tier2NumSlot = tier2NumSlot; }
    void                    SetTier1TxTime( const double tier1NumSlot ){ m_tier1TxTime = tier1NumSlot;}
    ClusterStructure*       CreateClusterStructure();
  private:
    bool Kmedoid( vector<int>&, list<list<int> >& );

    ULSA4b7_DC*   m_ptrToolSA; 
    std::string   m_mapFileName;
    FILE*         m_fid;
    double        m_compressionRatio;
    double        m_fidelityRatio;
    int                   m_tier2NumSlot;
    double                m_tier1TxTime;
};

#endif
