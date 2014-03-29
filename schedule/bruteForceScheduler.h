#ifndef _BRUTEFORCESCHEDULER_
#define _BRUTEFORCESCHEDULER_

#include "scheduler.h"
class BruteForceScheduler: public Scheduler
{
  public:
    BruteForceScheduler(const double txTime, 
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE* , 
        ClusterStructure const * const);
    ~BruteForceScheduler();

    void SetGaussianField(CORRE_MA_OPE* myGField) { m_ptrMatComputer = myGField;}
    double ScheduleOneSlot( std::vector<int>& );
    bool ScheduleOneSlot( std::vector<int>& , const std::vector<double>& );
    string PrintSelf(){ return m_type; }
    double GetTxTimePerSlot() const { return m_txTimePerSlot; }

  private:

    void Init();
    double  OmegaValue( const int nodeName );
    void    Perm(
        Eigen::MatrixXd& matX, 
        std::vector<int>& vecSupport,
        const int& ChIdx, 
        double& maxValue, 
        std::vector<int>& vecSolution,
        const std::vector<double>& vecVariance);

    bool EigenMatrixIsSmaller(const Eigen::MatrixXd&, const Eigen::MatrixXd& );
    
    int                             m_numNodes;
    int                             m_numMaxHeads;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    double                          m_maxPower;
    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE*                   m_ptrMatComputer;
    std::vector<double>             m_vecNodePower;
    std::vector<int>                m_vecSched; /* record the scheduled node */
    const string                    m_type;
    
    Eigen::MatrixXd                 m_A;
    Eigen::MatrixXd                 m_B;
    Eigen::MatrixXd                 m_C;
    Eigen::MatrixXd                 m_constraints;
    Eigen::MatrixXd                 m_X;
    Eigen::MatrixXd                 m_Signma;

}; 
#endif
