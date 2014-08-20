#ifndef _MINPOWERIMAGESCHEDULER_
#define _MINPOWERIMAGESCHEDULER_ 
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>
#include <coin/CoinPragma.hpp>
#include <coin/CoinTime.hpp>
#include <coin/CoinError.hpp>

#include <coin/BonOsiTMINLPInterface.hpp>
#include <coin/BonIpoptSolver.hpp>
#include <coin/BonCbc.hpp>
#include <coin/BonBonminSetup.hpp>

#include <coin/BonOACutGenerator2.hpp>
#include <coin/BonEcpCuts.hpp>
#include <coin/BonOaNlpOptim.hpp>
#include "scheduler.h"
#include "imageSource.h"
#include <vector>

#include "minPowerImageMILP.hpp"

class MinPowerImageScheduler: public Scheduler
{
  public:
    MinPowerImageScheduler(
        const double txTime, 
        const int tier2NumSlot,
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE* , 
        ClusterStructure const * const);
    ~MinPowerImageScheduler();

    void SetGaussianField(CORRE_MA_OPE* myGField) { m_ptrMatComputer = myGField;}
    double ScheduleOneSlot( std::vector<int>& );
    double ScheduleOneSlot( std::vector<int>& , std::vector<double>&, const std::vector<double>& );
    string PrintSelf(){ return m_type; }
    double GetTxTimePerSlot() const { return m_txTimePerSlot; }
  private:

    void InitSolution(); /* solve the problem */

    int                               m_numNodes;
    int                               m_numHeads;
    const double                      m_txTimePerSlot;
    const double                      m_bandwidthKhz;
    double                            m_maxPower;
    const int                         m_tier2NumSlot;
    Map const * const                 m_ptrMap;
    ClusterStructure const * const    m_ptrCS;
    CORRE_MA_OPE*                     m_ptrMatComputer;
    std::vector<double>               m_vecNodePower;
    const string                      m_type;

    Ipopt::SmartPtr<Bonmin::TMINLP>   m_MILP;

    std::vector<double>               m_vecSolution;
    int                               m_slotCounter;

};

#endif
