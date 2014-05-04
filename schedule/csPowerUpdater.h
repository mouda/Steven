#ifndef _CSPOWERUPDATER_
#define _CSPOWERUPDATER_

class CSPowerUpdater
{
  public:
    CSPowerUpdater(
        Map const * const,
        ClusterStructure const * const,
        const double txTimeSlot
        ); 
    ~CSPowerUpdater();
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
    double              m_idtEntropy;
    double              m_txTimePerSlot;
    bool                m_exceedPc;
    static const double m_scale = 1; //This scale is for avoiding computation error.
    int                 **m_maIndexInterference;
    double              **m_maStrengthInterference;

    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
};

#endif
