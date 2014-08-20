#ifndef _CSPOWERUPDATER_
#define _CSPOWERUPDATER_
#include "map.h"
#include "clusterStructure.h"
#include <vector>
#include <iostream>

class CSPowerUpdater
{
  public:
    CSPowerUpdater(
        Map const * const
        ); 
    ~CSPowerUpdater();
    double PowerUpdate( std::vector<double>&, const std::vector<int>&, ClusterStructure const * const );
    double Solve_withT2Adj_BinerySearch_2( double inIniT2, std::vector<double>&, const std::vector<int>&, ClusterStructure const * const ); 
    void  showVerificationResult(std::vector<double>&, const std::vector<int>&, ClusterStructure const * const, std::vector <double>& , std::vector<double>& v, std::vector<double>& );
    void    decideDiscard3b();

  private:

    void Init();
    void UpdateInterference( std::vector<double>&, const std::vector<int>&, ClusterStructure const * const );
    void ChangeAllMemberPower( std::vector<double>&, std::vector<double>&, std::vector<double>&, ClusterStructure const * const);
    bool CheckDifference( const std::vector<double>&, ClusterStructure const * const ) const;
    bool CheckConverged( const std::vector<double>&, ClusterStructure const * const ) ;
    void returnMaxNextNodePower(int &nodeIndex, double &nodePower, std::vector<double>&, ClusterStructure const * const);
    int returnHeadIndex_ByNodeName(int nodeName, ClusterStructure const * const);
  
    int                 m_criticalNode;
    int                 m_criticalHeadIndex;
    int                 m_statusFlag;
    const double        m_threshold;
    double              m_inBandNoise;
    double              m_avgRatio;
    double              m_C2;
    double              m_idtEntropy;
    double              m_best2nd_ms;
    bool                m_exceedPc;
    static const double m_scale = 1; //This scale is for avoiding computation error.
    int**               m_maIndexInterference;
    double**            m_maStrengthInterference;

    Map const * const               m_ptrMap;
};

#endif
