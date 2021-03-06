#ifndef _BBSECHEDULER_
#define _BBSECHEDULER_ 
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>
//#include <coin/IpIpoptApplication.hpp>
//#include <coin/IpSolveStatistics.hpp>
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

#include <vector>
#include <cmath>
#include "scheduler.h"
#include "MyTMINLP.hpp"



class BranchBoundScheduler: public Scheduler
{
  public:
    BranchBoundScheduler(const double txTime, 
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE* , 
        ClusterStructure const * const);
    ~BranchBoundScheduler();

    void    SetEpsilon(const double epsilon) {m_epsilon = epsilon;}
    void    SetGaussianField(CORRE_MA_OPE* myGField) { m_ptrMatComputer = myGField;}
    double  ScheduleOneSlot( std::vector<int>& );
    double  ScheduleOneSlot( std::vector<int>& , std::vector<double>&, const std::vector<double>& );
    string  PrintSelf(){ return m_type; }
    double GetTxTimePerSlot() const { return m_txTimePerSlot; }
  private:
    double OmegaValue( const int nodeName );
    int  SolverHook(std::vector<int>&, Eigen::MatrixXd&  );
    int                             m_numNodes;
    int                             m_numMaxHeads;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    double                          m_maxPower;
    double                          m_epsilon;
    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE*                   m_ptrMatComputer;
    std::vector<double>             m_vecNodePower;
    std::vector<int>                m_vecSched;
    const string                    m_type;
    Ipopt::SmartPtr<Bonmin::TMINLP> m_nlp;
    Eigen::MatrixXd                 m_A;
    Eigen::MatrixXd                 m_B;
    Eigen::MatrixXd                 m_C;
    Eigen::MatrixXd                 m_X;
    Eigen::MatrixXd                 m_Signma;

};
#endif
