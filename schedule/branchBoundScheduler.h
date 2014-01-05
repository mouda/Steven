#ifndef _BBSECHEDULER_
#define _BBSECHEDULER_ 
#include "scheduler.h"
#include <eigen3/Eigen/Core>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include "MyNLP.h"


using namespace Ipopt;

class BranchBoundScheduler: public Scheduler
{
  public:
    BranchBoundScheduler(const double txTime, 
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE const * const, 
        ClusterStructure const * const);
    ~BranchBoundScheduler();
    double ScheduleOneSlot( vector<int>& );
    string PrintSelf(){ return m_type; }
  private:
    double OmegaValue( const int nodeName );
    int                             m_numNodes;
    int                             m_numMaxHeads;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    double                          m_maxPower;
    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE const * const      m_ptrMatComputer;
    vector<double>                  m_vecNodePower;
    const string                    m_type;
    SmartPtr<TNLP>                  m_nlp;
    SmartPtr<IpoptApplication>      m_nlpApp;
    Eigen::MatrixXd                 m_A;
    Eigen::MatrixXd                 m_B;
    Eigen::MatrixXd                 m_C;
    Eigen::MatrixXd                 m_X;

};
#endif
