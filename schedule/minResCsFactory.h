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
    ClusterStructure*       CreateClusterStructure();
  private:
    bool Kmedoid( vector<int>&, list<list<int> >& );

    ULSA4b7_DC*   m_ptrToolSA; 
    std::string   m_mapFileName;
    FILE*         m_fid;
    double        m_compressionRatio;
    double        m_fidelityRatio;
};

#endif
