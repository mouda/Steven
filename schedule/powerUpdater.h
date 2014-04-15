#ifndef _POWERUPDATER_
#define _POWERUPDATER_

#include "map.h"
#include "clusterStructure.h"
#include <vector>

class PowerUpdater
{
  public:
    PowerUpdater(
        Map const * const,
        ClusterStructure const * const
        ); 
    ~PowerUpdater();

    double Solve( const std::vector<int>& );

  private:

    void Init();
    void UpdateInterference( std::vector<double>& );
    void ChangeAllMemberPower( std::vector<double>& ) const;
    bool CheckDifference( const std::vector<double>& ) const;
    bool CheckConverged( const std::vector<double>& ) ;

    const double        m_threshold;
    double              m_avgRatio;
    static const double m_scale = 1; //This scale is for avoiding computation error.

    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;


};
#endif
