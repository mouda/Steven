#ifndef _NEWIMAGEPOWERUPDATER_
#define _NEWIMAGEPOWERUPDATER_

#include "imageMap.h"
#include "clusterStructure.h"
#include <vector>
#include <iostream>

class NewImagePowerUpdater
{
  public:
    NewImagePowerUpdater(
        ImageMap const * const,
        ClusterStructure const * const,
        const double txTimeSlot,
        const double txNumSlot
        ); 
    ~NewImagePowerUpdater();

    double Solve( std::vector<double>&, const std::vector<int>& );

  private:

    void Init();
    void UpdateInterference( std::vector<double>&, const std::vector<int>& );
    void ChangeAllMemberPower( std::vector<double>&, std::vector<double>&, std::vector<double>& );
    bool CheckDifference( const std::vector<double>& ) const;
    bool CheckConverged( const std::vector<double>& ) ;

    const double        m_threshold;
    double              m_inBandNoise;
    double              m_avgRatio;
    double              m_C2;
    std::vector<double> m_vecC2;
    double              m_txTimePerSlot;
    bool                m_exceedPc;
    static const double m_scale = 1; //This scale is for avoiding computation error.
    int                 **m_maIndexInterference;
    double              **m_maStrengthInterference;

    ImageMap const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;


};
#endif
